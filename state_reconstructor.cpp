
#include <string>
#include <vector>

using namespace std;

#include "state_reconstructor.h"
#include "tree.h"
#include "rate_model.h"
#include "node.h"
#include "utils.h"
#include "vector_node_object.h"
#include "sequence.h"

StateReconstructor::StateReconstructor(RateModel & _rm):tree(NULL),nstates(_rm.nstates),rm(_rm),dc("dist_conditionals"),
		store_p_matrices(false),use_stored_matrices(false),revB("revB"),
		rev(false),rev_exp_number("rev_exp_number"),rev_exp_time("rev_exp_time"),
		stochastic(false),stored_EN_matrices(map<double, mat >()),
		stored_ER_matrices(map<double, mat >()),sp_alphas("sp_alphas"),alphas("alphas"){}
	/*
	 * initialize each node with segments
	 */
void StateReconstructor::set_tree(Tree * tr){
	tree = tr;
	cout << "initializing nodes..." << endl;
	for(int i=0;i<tree->getNodeCount();i++){
		if(tree->getNode(i)->getBL()<0.000001)
			tree->getNode(i)->setBL(0.000001);
		VectorNodeObject<double> * dcs = new VectorNodeObject<double>(nstates);
		tree->getNode(i)->assocObject(dc,*dcs);
		delete dcs;
	}
}

void StateReconstructor::set_tip_conditionals(vector<Sequence> & distrib_data){
	for(unsigned int i=0;i<distrib_data.size();i++){
		Sequence seq = distrib_data[i];
		Node * nd = tree->getExternalNode(seq.get_id());
		cout << nd->getName() << " ";
		for(int j=0;j<nstates;j++){
			if(seq.get_sequence().at(j) == '1')
				(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 1.0;
			else
				(((VectorNodeObject<double>*) nd->getObject(dc)))->at(j) = 0.0;
			cout << seq.get_sequence().at(j);
		}
		cout << endl;
	}
}

VectorNodeObject<double> StateReconstructor::conditionals(Node & node){
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
		for ( int i=0;i<nstates;i++){
			distconds.at(i)= v1[i] * v2[i];
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
	reverse(tree->getRoot());
}

void StateReconstructor::reverse(Node * node){
	rev = true;
	VectorNodeObject<double> * revconds = new VectorNodeObject<double> (nstates, 0);//need to delete this at some point
	if (node == tree->getRoot()) {
		for(int i=0;i<nstates;i++){
			revconds->at(i) = 1.0;//prior
		}
		node->assocObject(revB,*revconds);
		delete revconds;
		for(int i = 0;i<node->getChildCount();i++){
			reverse(node->getChild(i));
		}
	}else{
	//else if(node.isExternal() == false){
		//calculate A i
		//sum over all alpha k of sister node of the parent times the priors of the speciations
		//(weights) times B of parent j
		VectorNodeObject<double> * parrev = ((VectorNodeObject<double>*)node->getParent()->getObject(revB));
		VectorNodeObject<double> sisdistconds;
		if(node->getParent()->getChild(0) != node){
			VectorNodeObject<double>* talph = ((VectorNodeObject<double>*) node->getParent()->getChild(0)->getObject(alphas));
			sisdistconds = *talph;
		}else{
			VectorNodeObject<double>* talph = ((VectorNodeObject<double>*) node->getParent()->getChild(1)->getObject(alphas));
			sisdistconds = *talph;
		}

		VectorNodeObject<double> tempA (nstates,0);
		//needs to be the same as ancdist_cond_lh
		for ( int i = 0; i < nstates; i++) {
			//root has i, curnode has left, sister of cur has right
			//for ( int j = 0; j < nstates; j++) {
				tempA[i] += (sisdistconds.at(i)*parrev->at(i));
			//}
		}
		//now calculate node B
		//VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		vector<double> tempmoveA(tempA);
		for(int j=0;j<nstates;j++){revconds->at(j) = 0;}
		//RateModel * rm = tsegs->at(ts).getModel();
		cx_mat * p = &rm.stored_p_matrices[node->getBL()];
		mat * EN = NULL;
		mat * ER = NULL;
		VectorNodeObject<double> tempmoveAer(tempA);
		VectorNodeObject<double> tempmoveAen(tempA);
		if(stochastic == true){
			//initialize the segment B's
			for( int j=0;j<nstates;j++){tempmoveAer[j] = 0;}
			for( int j=0;j<nstates;j++){tempmoveAen[j] = 0;}
			EN = &stored_EN_matrices[node->getBL()];
			ER = &stored_ER_matrices[node->getBL()];
		}
		for( int j=0;j < nstates;j++){
			for ( int i = 0; i < nstates; i++) {
				revconds->at(j) += tempmoveA[i]*real((*p)(i,j));//tempA needs to change each time
				if(stochastic == true){
					tempmoveAer[j] += tempmoveA[i]*(((*ER)(i,j)));
					tempmoveAen[j] += tempmoveA[i]*(((*EN)(i,j)));
				}
			}
		}
		for( int j=0;j<nstates;j++){tempmoveA[j] = revconds->at(j);}
		if(stochastic == true){
			node->seg_sp_stoch_map_revB_time = tempmoveAer;
			node->seg_sp_stoch_map_revB_number = tempmoveAen;
		}

		node->assocObject(revB,*revconds);
		delete revconds;
		for(int i = 0;i<node->getChildCount();i++){
			reverse(node->getChild(i));
		}
	}
}

vector<double> StateReconstructor::calculate_ancstate_reverse(Node & node){
	if (node.isExternal()==false){//is not a tip
		VectorNodeObject<double> * Bs = (VectorNodeObject<double> *) node.getObject(revB);
		Node * c1 = node.getChild(0);
		Node * c2 = node.getChild(1);
		VectorNodeObject<double>* v1  = ((VectorNodeObject<double>*) c1->getObject(alphas));
		VectorNodeObject<double>* v2 = ((VectorNodeObject<double>*) c2->getObject(alphas));
		vector<double> LHOODS (nstates,0);
		for ( int i = 0; i < nstates; i++) {
			//for ( int j=0;j<nstates;j++){
			//	LHOODS[i] += (v1->at(i)*v2->at(j));//*weight);
			//}
			LHOODS[i] = (v1->at(i)*v2->at(i)) * Bs->at(i);
		}
		return LHOODS;
	}
}

void StateReconstructor::prepare_stochmap_reverse_all_nodes(int from, int to){
	stochastic = true;
	//calculate and store local expectation matrix for each branch length
	for(int k = 0; k < tree->getNodeCount(); k++){
		double dur =  tree->getNode(k)->getBL();
		cx_mat eigvec(nstates,nstates);eigvec.fill(0);
		cx_mat eigval(nstates,nstates);eigval.fill(0);
		bool isImag = rm.get_eigenvec_eigenval_from_Q(&eigval, &eigvec);
		mat Ql(nstates,nstates);Ql.fill(0);Ql(from,to) = rm.get_Q()(from,to);
		mat W(nstates,nstates);W.fill(0);W(from,from) = 1;
		cx_mat summed(nstates,nstates);summed.fill(0);
		cx_mat summedR(nstates,nstates);summedR.fill(0);
		for(int i=0;i<nstates;i++){
			mat Ei(nstates,nstates);Ei.fill(0);Ei(i,i)=1;
			cx_mat Si(nstates,nstates);
			Si = eigvec * Ei * inv(eigvec);
			for(int j=0;j<nstates;j++){
				cx_double dij = (eigval(i,i)-eigval(j,j)) * dur;
				mat Ej(nstates,nstates);Ej.fill(0);Ej(j,j)=1;
				cx_mat Sj(nstates,nstates);
				Sj = eigvec * Ej * inv(eigvec);
				cx_double Iijt = 0;
				if (abs(dij) > 10){
					Iijt = (exp(eigval(i,i)*dur)-exp(eigval(j,j)*dur))/(eigval(i,i)-eigval(j,j));
				}else if(abs(dij) < 10e-20){
					Iijt = dur*exp(eigval(j,j)*dur)*(1.+dij/2.+pow(dij,2.)/6.+pow(dij,3.)/24.);
				}else{
					if(eigval(i,i) == eigval(j,j)){
						//WAS Iijt = dur*exp(eigval(j,j)*dur)*expm1(dij)/dij;
						if (isImag)
							Iijt = dur*exp(eigval(j,j)*dur)*(exp(dij)-1.)/dij;
						else
							Iijt = dur*exp(eigval(j,j)*dur)*(expm1(real(dij)))/dij;
					}else{
						//WAS Iijt = -dur*exp(eigval(i,i)*dur)*expm1(-dij)/dij;
						if (isImag)
							Iijt = -dur*exp(eigval(i,i)*dur)*(exp(-dij)-1.)/dij;
						else
							Iijt = -dur*exp(eigval(i,i)*dur)*(expm1(real(-dij)))/dij;
					}
				}
				summed += (Si  * Ql * Sj * Iijt);
				summedR += (Si * W * Sj * Iijt);
			}
		}
		stored_EN_matrices[dur] = (real(summed));
		stored_ER_matrices[dur] = (real(summedR));
	}
}

vector<double> StateReconstructor::calculate_reverse_stochmap(Node & node, bool tm){
	if (node.isExternal()==false){//is not a tip
		vector<double> totalExp (nstates,0);
		vector<double> Bs;
		if(tm)
			Bs = node.seg_sp_stoch_map_revB_time;
		else
			Bs =  node.seg_sp_stoch_map_revB_number;
		Node * c1 = node.getChild(0);
		Node * c2 = node.getChild(1);
		VectorNodeObject<double> * v1  = ((VectorNodeObject<double>*) c1->getObject(alphas));
		VectorNodeObject<double> * v2  = ((VectorNodeObject<double>*) c2->getObject(alphas));
		VectorNodeObject<double> LHOODS (nstates,0);
		for ( int i = 0; i < nstates; i++) {
			//for (int j=0;j<nstates;j++){
				//int ind1 = leftdists[j];
				//int ind2 = rightdists[j];
				//LHOODS[i] += (v1.at(ind1)*v2.at(ind2)*weight);
			//}
			LHOODS[i] = v1->at(i) * v2->at(i) * Bs.at(i);
			//cout << v1->at(i) << " " <<  v2->at(i)<< " " << Bs.at(i) << endl;
		}
		for(int i=0;i<nstates;i++){
			totalExp[i] = LHOODS[i];
		}
		//not sure if this should return a double or not when doing a bigtree
		return totalExp;
	}else{
		vector<double> totalExp (nstates,0);
		vector<double> Bs;
		if(tm)
			Bs = node.seg_sp_stoch_map_revB_time;
		else
			Bs =  node.seg_sp_stoch_map_revB_number;
		VectorNodeObject<double> LHOODS (nstates,0);
		VectorNodeObject<double>* distconds = ((VectorNodeObject<double>*) node.getObject(dc));
		for (int i = 0; i < nstates; i++) {
			LHOODS[i] = Bs.at(i) * (distconds->at(i) );
		}
		for(int i=0;i<nstates;i++){
			totalExp[i] = LHOODS[i];
		}
		return totalExp;
	}
}

void StateReconstructor::set_store_p_matrices(bool i){
	store_p_matrices = i;
}

void StateReconstructor::set_use_stored_matrices(bool i){
	use_stored_matrices = i;
}

