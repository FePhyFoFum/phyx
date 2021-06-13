#ifndef PX_CONT_MODELS_H
#define PX_CONT_MODELS_H

#include <vector>
#include <map>
#include <cmath>

#include "tree.h"

class Node; // forward declaration

#include <armadillo>

using namespace arma;

double norm_pdf_multivariate (rowvec& x, rowvec& mu, mat& sigma);
double norm_log_pdf_multivariate (rowvec& x, rowvec& mu, mat& sigma);
Node * getMRCA_forVCV (Node * curn1,Node * curn2);
Node * getMRCAFromPath_forVCV (std::vector<Node *> * path1,Node * curn2);
void calc_vcv(Tree * tr, mat& vcv);
void calc_square_change_anc_states (Tree * tree, int index);
void calc_postorder_square_change (Node * node, std::map<Node *, int>& nodenum,
    mat * fullMcp, mat * fullVcp, int index);
double calc_bm_node_postorder (Node * node, int nch, double sigma);
double calc_bm_prune (Tree * tr, double sigma);

#endif /* PX_CONT_MODELS_H */
