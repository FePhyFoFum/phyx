
#ifndef TREE_RECONSTRUCTOR_H_
#define TREE_RECONSTRUCTOR_H_

#include "tree.h"
#include "node.h"
#include "rate_model.h"
#include "vector_node_object.h"
#include "sequence.h"

class StateReconstructor{

private:
	Tree * tree;
	int nstates;
	RateModel & rm;
	string dc;
	bool store_p_matrices;
	bool use_stored_matrices;

	//reverse bits
	string revB;
	bool rev;
	//end reverse bits

	//stochastic mapping bits
	string rev_exp_number;
	string rev_exp_time;
	bool stochastic;
	//map of period int and then branch length double
	map<double, mat > stored_EN_matrices;
	map<double, mat > stored_ER_matrices;
	//end mapping bits
	string sp_alphas;
	string alphas;

	VectorNodeObject<double> conditionals(Node & node);
	void ancdist_conditional_lh(Node & node);

public:
	StateReconstructor(RateModel &);
	void set_tree(Tree *);
	double eval_likelihood();
	void set_tip_conditionals(vector<Sequence> & distrib_data);
	void prepare_ancstate_reverse();
	void reverse(Node &);
	vector<double> calculate_ancstate_reverse(Node & node);
	void prepare_stochmap_reverse_all_nodes(int, int);
	vector<double> calculate_reverse_stochmap(Node &, bool);
	void set_store_p_matrices(bool i);
	void set_use_stored_matrices(bool i);
};
#endif
