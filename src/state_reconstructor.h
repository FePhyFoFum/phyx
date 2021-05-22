#ifndef PX__STATE_RECONSTRUCTOR_H
#define PX__STATE_RECONSTRUCTOR_H

#include <vector>
#include <string>
#include <map>

#include "rate_model.h"
#include "superdouble.h"
#include "vector_node_object.h"

class Sequence; // forward declaration
class Node; // forward declaration
class Tree; // forward declaration


class StateReconstructor {
private:
    Tree * tree;
    std::vector<double> periods;
    bool use_periods;
    
    int nstates;
    RateModel& rm;
    std::vector<RateModel>& rm_periods;
    std::string dc;
    std::string andc;
    bool store_p_matrices;
    bool use_stored_matrices;
    
    //reverse bits
    std::string revB;
    bool rev;
    //end reverse bits
    
    //stochastic mapping bits
    std::string rev_exp_number;
    std::string rev_exp_time;
    bool stochastic;
    //map of period int and then branch length double
    std::map<Superdouble, mat > stored_EN_matrices;
    std::map<Superdouble, mat > stored_ER_matrices;
    //end mapping bits
    std::string sp_alphas;
    std::string alphas;
    
    VectorNodeObject<Superdouble> conditionals (Node& node);
    VectorNodeObject<Superdouble> conditionals_periods (Node& node);
    void ancdist_conditional_lh (Node& node);
    
public:
    StateReconstructor (RateModel&, std::vector<RateModel>& _vrm);
    void set_periods (std::vector<double>& ps, std::vector<RateModel>& rms);
    void set_tree (Tree *);
    double eval_likelihood ();
    void set_periods_model ();
    bool set_tip_conditionals (std::vector<Sequence>& distrib_data);
    bool set_tip_conditionals_already_given (std::vector<Sequence>& distrib_data);
    void prepare_ancstate_reverse ();
    void reverse (Node *);
    std::vector<Superdouble> calculate_ancstate_reverse_sd (Node& node);
    std::vector<double> calculate_ancstate_reverse (Node& node);
    void prepare_stochmap_reverse_all_nodes (int, int);
    void prepare_stochmap_reverse_all_nodes_all_matrices ();
    std::vector<double> calculate_reverse_stochmap (Node&, bool);
    void set_store_p_matrices (bool i);
    void set_use_stored_matrices (bool i);
    ~StateReconstructor ();
};

#endif /* PX__STATE_RECONSTRUCTOR_H */
