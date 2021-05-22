#ifndef PX__OPTIMIZE_CONT_MODELS_NLOPT_H
#define PX__OPTIMIZE_CONT_MODELS_NLOPT_H

#include <vector>
#include <nlopt.h>

#include "cont_models.h"

class Tree; // forward declaration

#include <armadillo>
using namespace arma;

double nlopt_bm_sr(unsigned n, const double *x, double *grad, void *data);
double nlopt_bm_sr_log(unsigned n, const double *x, double *grad, void *data);
double nlopt_ou_sr_log(unsigned n, const double *x, double *grad, void *data);
double nlopt_bm_bl(unsigned n, const double *x, double *grad, void *data);
std::vector<double> optimize_single_rate_bm_nlopt(rowvec& _x, mat& _vcv, bool log);
std::vector<double> optimize_single_rate_bm_ou_nlopt(rowvec& _x, mat& _vcv);
std::vector<double> optimize_single_rate_bm_bl(Tree * tr);

#endif /* PX__OPTIMIZE_CONT_MODELS_NLOPT_H */
