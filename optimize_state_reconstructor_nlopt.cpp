/*
 * optimize_tnc.cpp
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#include "optimize_nlopt.h"
#include <iostream>
#include <stdio.h>
#include <nlopt.hpp>
#include <math.h>

#include "state_reconstructor.h"
#include "rate_model.h"

#include <armadillo>
using namespace arma;


StateReconstructor * nloptsr;
RateModel * nloptrm;
mat * nloptfree_variables;

//double nlopt_sr(int n, const double *x, void *state){
//double nlopt_sr(unsigned n, const double *x, double *grad, void *my_func_data);
double nlopt_sr(unsigned n, const double *x, double *grad, void *my_func_data){
	for (unsigned int i=0;i<nloptfree_variables->n_rows;i++){
		for (unsigned int j=0;j<nloptfree_variables->n_cols;j++){
			if(i!=j){
				nloptrm->set_Q_cell(i,j,x[int((*nloptfree_variables)(i,j))]);
				if(nloptrm->get_Q()(i,j) < 0 || nloptrm->get_Q()(i,j) >= 1000){
					return 100000000;
				}
			}
		}
	}
	double like;
	nloptrm->set_Q_diag();
	like = nloptsr->eval_likelihood();
	//cout << like << endl;
	if(like < 0 || like == std::numeric_limits<double>::infinity())
		like = 100000000;
	return like;
}

void optimize_sr_nlopt(RateModel * _rm,StateReconstructor * _sr, mat * _free_mask, int _nfree){
	nloptsr = _sr;
	nloptrm = _rm;
	nloptfree_variables = _free_mask;

	//nlopt_opt opt;
	nlopt::opt opt(nlopt::LN_SBPLX, _nfree);


	opt.set_lower_bounds(0.0);
	opt.set_upper_bounds(1000);
	opt.set_min_objective(nlopt_sr,NULL);
	opt.set_xtol_rel(0.001);

	vector<double> x(_nfree,0);
	for(unsigned int i=0;i<_rm->get_Q().n_rows;i++){
		for(unsigned int j=0;j<_rm->get_Q().n_cols;j++){
			if(i!=j){
				x[int((*_free_mask)(i,j))] = _rm->get_Q()(i,j);
				//cout << x[int((*_free_mask)(i,j))] << " ";
			}
		}
		//cout << endl;
	}

	//double minf;
	vector<double> result = opt.optimize(x);
	for(unsigned int i=0;i<_rm->get_Q().n_rows;i++){
		for(unsigned int j=0;j<_rm->get_Q().n_cols;j++){
			if(i!=j){
				(*_free_mask)(i,j) = result[int((*_free_mask)(i,j))];
				//cout << x[int((*_free_mask)(i,j))] << " ";
			}
		}
		//cout << endl;
	}
	//if ((&minf) < 0) {
	 //   printf("nlopt failed!\n");
	//}
	//else {
	//    printf("found minimum at %0.10g\n", x[0], x[1], minf);
	//}
}



