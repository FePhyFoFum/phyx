
#ifndef TREE_RECONSTRUCTOR_SIMPLE_H_
#define TREE_RECONSTRUCTOR_SIMPLE_H_

#include <map>
#include <string>
#include <vector>

using namespace std;

#include "tree.h"
#include "node.h"
#include "rate_model.h"
#include "vector_node_object.h"
#include "sequence.h"

class StateReconstructorSimple{

private:
    Tree * tree;
    int nstates;
    int nsites;
    RateModel & rm;
    string dc;
    bool store_p_matrices;
    bool use_stored_matrices;
    
    string sp_alphas;
    string alphas;
    mat p;
    map<Node *,vector<double> > tip_conditionals;
    map<Node *,vector<double> > internal_conditionals;
    VectorNodeObject<double> v_storage;  
    VectorNodeObject<double> conditionals(Node & node);
    void ancdist_conditional_lh(Node & node);
    
public:
    StateReconstructorSimple(RateModel &, int);
    void set_tree(Tree *);
    double eval_likelihood();
    bool set_tip_conditionals(vector<Sequence> & distrib_data,int );
    void set_store_p_matrices(bool i);
    void set_use_stored_matrices(bool i);
    ~StateReconstructorSimple();
};
#endif
