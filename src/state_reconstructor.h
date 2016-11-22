
#ifndef _STATE_RECONSTRUCTOR_H_
#define _STATE_RECONSTRUCTOR_H_

#include "tree.h"
#include "node.h"
#include "rate_model.h"
#include "vector_node_object.h"
#include "sequence.h"
#include "superdouble.h"

using namespace std;

class StateReconstructor{

private:
    Tree * tree;
    vector<double> periods;
    bool use_periods;
    
    int nstates;
    RateModel & rm;
    vector<RateModel> & rm_periods;
    string dc;
    string andc;
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
    map<Superdouble, mat > stored_EN_matrices;
    map<Superdouble, mat > stored_ER_matrices;
    //end mapping bits
    string sp_alphas;
    string alphas;
    
    VectorNodeObject<Superdouble> conditionals(Node & node);
    VectorNodeObject<Superdouble> conditionals_periods(Node & node);
    void ancdist_conditional_lh(Node & node);
    
public:
    StateReconstructor(RateModel &, vector<RateModel> & _vrm);
    void set_periods(vector<double> & ps, vector<RateModel> & rms);
    void set_tree(Tree *);
    double eval_likelihood();
    void set_periods_model();
    bool set_tip_conditionals(vector<Sequence> & distrib_data);
    bool set_tip_conditionals_already_given(vector<Sequence> & distrib_data);
    void prepare_ancstate_reverse();
    void reverse(Node *);
    vector<Superdouble> calculate_ancstate_reverse_sd(Node & node);
    vector<double> calculate_ancstate_reverse(Node & node);
    void prepare_stochmap_reverse_all_nodes(int, int);
    void prepare_stochmap_reverse_all_nodes_all_matrices();
    vector<double> calculate_reverse_stochmap(Node &, bool);
    void set_store_p_matrices(bool i);
    void set_use_stored_matrices(bool i);
    ~StateReconstructor();
};
#endif /* _STATE_RECONSTRUCTOR_H_ */
