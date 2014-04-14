
#include <string>
#include <vector>
#include <map>

using namespace std;

#include "state_reconstructor_simple.h"
#include "tree.h"
#include "rate_model.h"
#include "node.h"
#include "utils.h"
#include "vector_node_object.h"
#include "sequence.h"

#define verbose false

#define MINBL 0.000000001

StateReconstructorSimple::StateReconstructorSimple(RateModel & _rm):tree(NULL),nstates(_rm.nstates),rm(_rm),dc("dist_conditionals"),store_p_matrices(false),use_stored_matrices(false),sp_alphas("sp_alphas"),alphas("alphas"){}
    /*
     * initialize each node with segments
     */
void StateReconstructorSimple::set_tree(Tree * tr){
    tree = tr;
    if(verbose)
        cout << "initializing nodes..." << endl;
    for(int i=0;i<tree->getNodeCount();i++){
        if(tree->getNode(i)->getBL()<MINBL){
            tree->getNode(i)->setBL(MINBL * 100);
	}
	VectorNodeObject<double> * dcs = new VectorNodeObject<double>(nstates);
	tree->getNode(i)->assocObject(dc,*dcs);
	delete dcs;
    }
}

bool StateReconstructorSimple::set_tip_conditionals(vector<Sequence> & distrib_data){
    bool allsame = true;
    string testsame = distrib_data[0].get_sequence();
    for(unsigned int i=0;i<distrib_data.size();i++){
	Sequence seq = distrib_data[i];
	Node * nd = tree->getExternalNode(seq.get_id());
	if(verbose)
	    cout << nd->getName() << " ";
	for(int j=0;j<nstates;j++){
	    if(seq.get_sequence().at(j) == '1')
		(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 1.0;
	    else
		(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 0.0;
	    if(verbose)
		cout << seq.get_sequence().at(j);
	}
	if(verbose)
	    cout << endl;
	if (testsame != seq.get_sequence())
	    allsame = false;
    }
    if (allsame == true && verbose == true){
	cerr << "all the tips have the same characters" << endl;
    }
    return allsame;
}

VectorNodeObject<double> StateReconstructorSimple::conditionals(Node & node){
    VectorNodeObject<double> distconds = *((VectorNodeObject<double>*) node.getObject(dc));
    VectorNodeObject<double> * v = new VectorNodeObject<double> (nstates, 0);
    cx_mat p;
    if(use_stored_matrices == false){
	p= rm.setup_P(node.getBL(),store_p_matrices);
    }else{
	p = rm.stored_p_matrices[node.getBL()];
    }
    for( int j=0;j<nstates;j++){
	for( int k=0;k<nstates;k++){
	    v->at(j) += (distconds.at(k)*real(p(j,k)));
	}
    }
    for(unsigned int j=0;j<distconds.size();j++){
	distconds[j] = v->at(j);
    }
    if(store_p_matrices == true){
	node.assocObject(sp_alphas,distconds);
	node.assocObject(alphas,distconds);
    }
    delete v;
    return distconds;
}

void StateReconstructorSimple::ancdist_conditional_lh(Node & node){
    VectorNodeObject<double> distconds(nstates, 0);
    if (node.isExternal()==false){//is not a tip
	Node * c1 = node.getChild(0);
	Node * c2 = node.getChild(1);
	ancdist_conditional_lh(*c1);
	ancdist_conditional_lh(*c2);
	VectorNodeObject<double> v1;
	VectorNodeObject<double> v2;
	v1 =conditionals(*c1);
	v2 =conditionals(*c2);
	for ( int i=0;i<nstates;i++){
	    distconds.at(i)= v1[i] * v2[i];
	}
	if(node.isRoot()){
	    //with equal freq!
	    for(int i=0;i<nstates;i++){
		distconds.at(i) = distconds.at(i) * (1./nstates);
	    }
	}
    }else{
	distconds = *((VectorNodeObject<double>*)node.getObject(dc));
    }
    for(unsigned int i=0;i<distconds.size();i++){
	((VectorNodeObject<double>*)node.getObject(dc))->at(i) = distconds.at(i);
    }
}

double StateReconstructorSimple::eval_likelihood(){
    ancdist_conditional_lh(*tree->getRoot());
    return (-log(calculate_vector_double_sum(*
          (VectorNodeObject<double>*) tree->getRoot()->getObject(dc))));
    //return double(-(calculate_vector_Superdouble_sum(*(VectorNodeObject<double>*) tree->getRoot()->getObject(dc))).getLn());
}


void StateReconstructorSimple::set_store_p_matrices(bool i){
    store_p_matrices = i;
}

void StateReconstructorSimple::set_use_stored_matrices(bool i){
    use_stored_matrices = i;
}

StateReconstructorSimple::~StateReconstructorSimple(){

}

