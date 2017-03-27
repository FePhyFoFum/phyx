/*
 * optimize_tnc.h
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#ifndef _OPTIMIZE_CONT_MODELS_NLOPT_H_
#define _OPTIMIZE_CONT_MODELS_NLOPT_H_

#include <nlopt.h>

#include "cont_models.h"

#include <armadillo>
using namespace arma;

double nlopt_bm_sr(unsigned n, const double *x, double *grad, void *data);
double nlopt_bm_sr_log(unsigned n, const double *x, double *grad, void *data);
double nlopt_ou_sr_log(unsigned n, const double *x, double *grad, void *data);
double nlopt_bm_bl(unsigned n, const double *x, double *grad, void *data);
vector<double> optimize_single_rate_bm_nlopt(rowvec & _x, mat & _vcv,bool log);
vector<double> optimize_single_rate_bm_ou_nlopt(rowvec & _x, mat & _vcv);
vector<double> optimize_single_rate_bm_bl(Tree * tr);

#endif /* _OPTIMIZE_CONT_MODELS_NLOPT_H_ */
