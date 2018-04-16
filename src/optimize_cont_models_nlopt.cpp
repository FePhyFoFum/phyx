/*
 * optimize_tnc.cpp
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#include <iostream>
#include <stdio.h>
#include <nlopt.hpp>
#include <math.h>

#include "cont_models.h"

#include <armadillo>
using namespace arma;

#define LARGE 1000000000

typedef struct {
    rowvec x;
    mat ovcv;
} analysis_data;

typedef struct {
    Tree * tree;
} analysis_data_tree;

double nlopt_bm_sr(unsigned n, const double *x, double *grad, void *data) {
    if (x[1] <= 0) {
        return LARGE;
    }
    //cout << x[0] << " " << x[1] << endl;
    analysis_data * d = (analysis_data *) data;
    mat tvcv = (d->ovcv) * x[1];
    rowvec m = rowvec(d->x.n_cols); m.fill(x[0]);
    double like = norm_pdf_multivariate(d->x, m, tvcv);
    return -like;
}

double nlopt_bm_sr_log(unsigned n, const double *x, double *grad, void *data) {
    if (x[1] <= 0) {
        return LARGE;
    }
    //cout << x[0] << " " << x[1] << endl;
    analysis_data * d = (analysis_data *) data;
    mat tvcv = (d->ovcv) * x[1];
    rowvec m = rowvec(d->x.n_cols); m.fill(x[0]);
    double like = norm_log_pdf_multivariate(d->x, m, tvcv);
    return -like;
}

/*
 * single alpha ou
 */
double nlopt_ou_sr_log(unsigned n, const double *x, double *grad, void *data) {
    if (x[1] <= 0 || x[2] <= 0) {
        return LARGE;
    }
    double alpha = x[2];
    analysis_data * d = (analysis_data *) data;
    mat vcvDiag(d->ovcv.n_cols,d->ovcv.n_cols); vcvDiag.zeros(); 
    vcvDiag.diag() = (d->ovcv).diag();
    mat tm(vcvDiag.n_cols,vcvDiag.n_cols);
    tm.ones();
    mat diagi = trans(vcvDiag * tm);
    mat diagj = vcvDiag * tm;
    mat Tij = diagi + diagj - (2 * d->ovcv);
    mat ouvcv = (1. / (2. * alpha)) * exp(-alpha * Tij) % (1. - exp(-2. * alpha * d->ovcv));
    ouvcv = ouvcv * x[1];
    rowvec m = rowvec(d->x.n_cols); m.fill(x[0]);
    double like = norm_log_pdf_multivariate(d->x,m,ouvcv);
    return -like;
}


double nlopt_bm_bl(unsigned n, const double *x, double *grad, void *data){
    for (unsigned int i=0;i<n;i++){
        if (x[i] <= 0){
            return LARGE;
        }   
    }
    double sigma =1;// x[0];//1;
    analysis_data_tree * d = (analysis_data_tree *) data;
    Tree * tr = d->tree;
    for (int i=0;i<tr->getNodeCount();i++){
        if (tr->getNode(i) != tr->getRoot()){
            tr->getNode(i)->setBL(x[i+1]);
        }
    }
    double like = calc_bm_prune(tr,sigma);
    //cout << like <<" " << sigma << endl;
    return -like;
}

vector<double> optimize_single_rate_bm_nlopt(rowvec & _x, mat & _vcv, bool log) {
    analysis_data a;
    a.x = _x;
    a.ovcv = _vcv;

    //nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);
    //nlopt::opt opt(nlopt::LN_BOBYQA,2);
    //BOBYQA is better but the other finishes more
    nlopt::opt opt(nlopt::LN_SBPLX, 2);
    //nlopt::opt opt(nlopt::LN_PRAXIS,2);

    opt.set_lower_bounds(0.000000001);
    opt.set_upper_bounds(100000);
    opt.set_ftol_abs(0.000001);
    if (log) {
        opt.set_min_objective(nlopt_bm_sr_log, &a);
    } else {
        opt.set_min_objective(nlopt_bm_sr, &a);
    }
    opt.set_xtol_rel(0.000001);
    opt.set_maxeval(5000);

    double minf;
    //2 parameters, 1 anc, 2 rate
    vector<double> x(2,1);
    nlopt::result result = opt.optimize(x, minf);
    vector<double> results;
    results.push_back(x[0]);
    results.push_back(x[1]);
    results.push_back(minf);
    return results;
}

vector<double> optimize_single_rate_bm_ou_nlopt(rowvec & _x, mat & _vcv) {
    analysis_data a;
    a.x = _x;
    a.ovcv = _vcv;

    //nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
    //BOBYQA is better but the other finishes more
    //nlopt::opt opt(nlopt::LN_BOBYQA,3);
    nlopt::opt opt(nlopt::LN_SBPLX, 3);
    //nlopt::opt opt(nlopt::LN_PRAXIS,3);
    opt.set_min_objective(nlopt_ou_sr_log, &a);
    opt.set_lower_bounds(0.000000001);
    opt.set_upper_bounds(100000);
    opt.set_xtol_rel(0.000001);
    opt.set_ftol_rel(0.00001);
    opt.set_maxeval(5000);
    double minf;
    //2 parameters, 1 anc, 2 rate, 3 alpha
    vector<double> x(3, 1);
    nlopt::result result = opt.optimize(x, minf);
//    cout << result << endl;
    vector<double> results;
    results.push_back(x[0]); results.push_back(x[1]); results.push_back(x[2]);
    results.push_back(minf);
    return results;
}

vector<double> optimize_single_rate_bm_bl(Tree * tr) {
    analysis_data_tree a;
    a.tree = tr;
    int n = 1+tr->getNodeCount() - 1;
    //nlopt::opt opt(nlopt::LN_NELDERMEAD, n);
    //BOBYQA is better but the other finishes more
    //nlopt::opt opt(nlopt::LN_BOBYQA,n);
    nlopt::opt opt(nlopt::LN_SBPLX,n);
    //nlopt::opt opt(nlopt::LN_COBYLA,n);
    //nlopt::opt opt(nlopt::LN_PRAXIS,n);
    opt.set_min_objective(nlopt_bm_bl, &a);
    opt.set_lower_bounds(0.0001);
    opt.set_upper_bounds(100000);
    opt.set_xtol_rel(0.000001);
    opt.set_ftol_rel(0.000001);
    opt.set_maxeval(100000);
    double minf;
    vector<double> x(n, 1);
    nlopt::result result = opt.optimize(x, minf);
    cout << result << endl;
    vector<double> results;
    for (int i=0;i<tr->getNodeCount();i++){
        if (tr->getNode(i) != tr->getRoot()){
            tr->getNode(i)->setBL(x[i+1]);
        }
    }
    results.push_back(x[0]); 
    results.push_back(minf);
    return results;
}


