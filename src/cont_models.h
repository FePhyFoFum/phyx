#ifndef _CONT_MODELS_H_
#define _CONT_MODELS_H_

#include <string>
#include <vector>
#include <math.h>


#include "node.h"
#include "tree.h"
#include "tree_utils.h"
#include <armadillo>

using namespace arma;


double norm_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma);
double norm_log_pdf_multivariate(rowvec & x, rowvec & mu, mat & sigma);
Node * getMRCA_forVCV(Node * curn1,Node * curn2);
Node * getMRCAFromPath_forVCV(vector<Node *> * path1,Node * curn2);
void calc_vcv(Tree * tr, mat & vcv);
void calc_square_change_anc_states(Tree * tree, int index);
void calc_postorder_square_change(Node * node,map<Node *,int> & nodenum,
    mat * fullMcp, mat * fullVcp, int index);
double calc_bm_node_postorder(Node * node, int nch, double sigma);
double calc_bm_prune(Tree * tr, double sigma);
#endif /* _CONT_MODELS_H_ */
