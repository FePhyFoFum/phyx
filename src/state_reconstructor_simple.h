
#ifndef _STATE_RECONSTRUCTOR_SIMPLE_H_
#define _STATE_RECONSTRUCTOR_SIMPLE_H_

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
    
    map<double, mat> map_ps;
    map<double, mat> map_ps0;
    map<double, mat> map_ps1;
    map<double, mat> map_ps2;
    map<Node *,vector<vector<double> > > conditionals_map;
    vector<double> v_storage;  //just junk storage
    vector<double> v1;
    vector<double> v2;
    
    //VectorNodeObject<double> 
    void conditionals(vector<double> * v, Node & node,int site);
    void conditionals2(vector<double> * v, Node & node,int site);
    void ancdist_conditional_lh(Node & node,int site);
    
public:
    StateReconstructorSimple(RateModel &, int);
    void set_tree(Tree *);
    double eval_likelihood(int site);
    bool set_tip_conditionals(vector<Sequence> & distrib_data,int );
    void clear_map_ps();
    double pp0;
    double pp1;
    double pp2;
    ~StateReconstructorSimple();
};
#endif /* _STATE_RECONSTRUCTOR_SIMPLE_H_ */
