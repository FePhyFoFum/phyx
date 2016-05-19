/*
 * optimize_tnc.cpp
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#include "optimize_state_reconstructor_periods_nlopt.h"
#include <iostream>
#include <stdio.h>
#include <nlopt.hpp>
#include <math.h>
#include <vector>

#include "state_reconstructor.h"
#include "rate_model.h"

#include <armadillo>
using namespace arma;


StateReconstructor * nloptsr_periods;
vector<RateModel> * nloptrm_periods;
vector<mat> * nloptfree_variables_periods;

//double nlopt_sr(int n, const double *x, void *state) {
//double nlopt_sr(unsigned n, const double *x, double *grad, void *my_func_data);
double nlopt_sr_periods(unsigned n, const double *x, double *grad, void *my_func_data) {
    for (unsigned int k=0; k < nloptfree_variables_periods->size(); k++) {
        for (unsigned int i=0; i < (*nloptfree_variables_periods)[k].n_rows; i++) {
            for (unsigned int j=0; j < (*nloptfree_variables_periods)[k].n_cols; j++) {
                if (i != j) {
                    (*nloptrm_periods)[k].set_Q_cell(i,j,x[int(nloptfree_variables_periods->at(k)(i,j))]);
                    if ((*nloptrm_periods)[k].get_Q()(i,j) < 0 || (*nloptrm_periods)[k].get_Q()(i,j) >= 1000) {
                        return 1000000000000;
                    }
                }
            }
        }
    }
    double like;
    for (unsigned int i=0; i < nloptrm_periods->size(); i++) {
        nloptrm_periods->at(i).set_Q_diag();
    }
    like = nloptsr_periods->eval_likelihood();
    for (unsigned int i=0; i < nloptrm_periods->size(); i++) {
        if (nloptrm_periods->at(i).neg_p == true) {
            like = 10000000000000;
            break;
        }
    }
//    cout << like << endl;
    if (like < 0 || like == std::numeric_limits<double>::infinity()) {
        like = 10000000000000;
    }
    return like;
}

void optimize_sr_periods_nlopt(vector<RateModel> * _rm,StateReconstructor * _sr, vector<mat> * _free_mask, int _nfree) {
    nloptsr_periods = _sr;
    nloptrm_periods = _rm;
    nloptfree_variables_periods = _free_mask;
    
    nlopt::opt opt(nlopt::LN_NELDERMEAD, _nfree);
    //nlopt::opt opt(nlopt::LN_BOBYQA, _nfree);
    //nlopt::opt opt(nlopt::LN_PRAXIS, _nfree);
    //nlopt::opt opt(nlopt::LN_SBPLX, _nfree);
    //nlopt::opt opt(nlopt::LN_COBYLA, _nfree);
    //nlopt::opt opt(nlopt::LN_NEWUOA, _nfree);
    
    opt.set_lower_bounds(0.0000);
    opt.set_upper_bounds(100000);
    opt.set_min_objective(nlopt_sr_periods, NULL);
    opt.set_xtol_rel(0.001);
    opt.set_maxeval(10000);
    
    vector<double> x(_nfree,0);
    for (unsigned int k=0; k < _rm->size(); k++) {
        for (unsigned int i=0; i < _rm->at(k).get_Q().n_rows; i++) {
            for (unsigned int j=0; j < _rm->at(k).get_Q().n_cols; j++) {
                if (i != j) {
                    x[int((*_free_mask)[k](i, j))] = _rm->at(k).get_Q()(i, j);
                    //cout << x[int((*_free_mask)[k](i,j))] << " ";
                }
            }
            cout << endl;
        }
    }
    //double minf;
    vector<double> result = opt.optimize(x);
    for (unsigned int k=0; k < _rm->size(); k++) {
        for (unsigned int i=0; i < _rm->at(k).get_Q().n_rows; i++) {
            for (unsigned int j=0; j < _rm->at(k).get_Q().n_cols; j++) {
                if (i != j) {
                    (*_free_mask)[k](i, j) = result[int((*_free_mask)[k](i, j))];
                    //cout << x[int((*(*_free_mask)[k])(i,j))] << " ";
                }
            }
        //cout << endl;
        }
    }
    //if ((&minf) < 0) {
    //   printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at %0.10g\n", x[0], x[1], minf);
    //}
}



