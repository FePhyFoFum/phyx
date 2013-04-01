#ifndef CONTMODELS_H_
#define CONTMODELS_H_

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

#endif
