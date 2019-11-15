
#ifndef _STATE_RECONSTRUCTOR_SIMPLE_H_
#define _STATE_RECONSTRUCTOR_SIMPLE_H_

#include <string>
#include <vector>
#include <map>

#include "tree.h"
#include "node.h"
#include "rate_model.h"
#include "vector_node_object.h"
#include "sequence.h"

class StateReconstructorSimple {

private:
    Tree * tree;
    int nstates;
    int nsites;
    RateModel& rm;
    std::string dc;
    
    std::map<double, mat> map_ps;
    std::map<double, mat> map_ps0;
    std::map<double, mat> map_ps1;
    std::map<double, mat> map_ps2;
    std::map<Node *, std::vector< std::vector<double> > > conditionals_map;
    std::vector<double> v_storage;  //just junk storage
    std::vector<double> v1;
    std::vector<double> v2;
    
    //VectorNodeObject<double> 
    void conditionals(std::vector<double> * v, Node& node, int site);
    void conditionals2(std::vector<double> * v, Node& node, int site);
    void ancdist_conditional_lh(Node& node, int site);
    
public:
    StateReconstructorSimple(RateModel&, int);
    void set_tree(Tree *);
    double eval_likelihood(int site);
    bool set_tip_conditionals(std::vector<Sequence>& distrib_data, int);
    void clear_map_ps();
    double pp0;
    double pp1;
    double pp2;
    ~StateReconstructorSimple();
};
#endif /* _STATE_RECONSTRUCTOR_SIMPLE_H_ */
