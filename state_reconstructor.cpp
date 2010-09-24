
#include <string>
#include <vector>

using namespace std;

#include "state_reconstructor.h"
#include "tree.h"
#include "rate_model.h"
#include "node.h"
#include "utils.h"
#include "vector_node_object.h"

StateReconstructor::StateReconstructor(RateModel & _rm):tree(NULL),nstates(_rm.nstates),rm(_rm),dc("dist_conditionals"),
		store_p_matrices(false),use_stored_matrices(false),revB("revB"),
		rev(false),rev_exp_number("rev_exp_number"),rev_exp_time("rev_exp_time"),
		stochastic(false),stored_EN_matrices(map<double, mat >()),
		stored_ER_matrices(map<double, mat >()),sp_alphas("sp_alphas"),alphas("alphas"){

}

void StateReconstructor::set_tip_conditionals(map<string,vector<int> > distrib_data){
	int numofleaves = tree->getExternalNodeCount();
	for(int i=0;i<numofleaves;i++){
		for(int j=0;j<nstates;j++){
			if(distrib_data[tree->getExternalNode(i)->getName()][j] == 1)
				(((VectorNodeObject<double>*) tree->getExternalNode(i)->getObject(dc)))->at(j) = 1.0;
		}
	}
}

VectorNodeObject<double> StateReconstructor::conditionals(Node & node){
	VectorNodeObject<double> distconds;
	distconds = *((VectorNodeObject<double>*) node.getObject(dc));
	VectorNodeObject<double> * v = new VectorNodeObject<double> (nstates, 0);
	cx_mat p;
	if(use_stored_matrices == false){
		p= rm.setup_P(node.getBL(),store_p_matrices);
	}else{
		p = rm.stored_p_matrices[node.getBL()];
	}
	for(unsigned int j=0;j<nstates;j++){
		for(unsigned int k=0;k<nstates;k++){
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

void StateReconstructor::ancdist_conditional_lh(Node & node){
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
			for (unsigned int i=0;i<nstates;i++){
				double lh_part1 = 0.0;
				double lh_part2 = 0.0;
				for (unsigned int j=0;j<nstates;j++){
					lh_part1 += v1.at(j);
					lh_part2 += v2.at(j);
				}
				distconds.at(i)= lh_part1 * lh_part2;
			}
		}else{
			distconds = *((VectorNodeObject<double>*)node.getObject(dc));
		}
		for(unsigned int i=0;i<distconds.size();i++){
			((VectorNodeObject<double>*)node.getObject(dc))->at(i) = distconds.at(i);
		}
	}

double StateReconstructor::eval_likelihood(){
	ancdist_conditional_lh(*tree->getRoot());
	return (-log(calculate_vector_double_sum(*
			(VectorNodeObject<double>*) tree->getRoot()->getObject(dc))));
}

void StateReconstructor::prepare_ancstate_reverse(){

}

void StateReconstructor::reverse(Node &){

}

vector<double> StateReconstructor::calculate_ancstate_reverse(Node & node){

}

void StateReconstructor::prepare_stochmap_reverse_all_nodes(int, int){

}

vector<double> StateReconstructor::calculate_reverse_stochmap(Node &, bool){

}
