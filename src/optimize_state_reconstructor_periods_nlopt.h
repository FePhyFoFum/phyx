/*
 * optimize_tnc.h
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#ifndef OPTIMIZE_NLOPT_PERIODS_H_
#define OPTIMIZE_NLOPT_PERIODS_H_

#include <nlopt.h>

#include "state_reconstructor.h"
#include "rate_model.h"

#include <armadillo>
using namespace arma;

void optimize_sr_periods_nlopt(vector<RateModel> * _rm,StateReconstructor * _sr, vector<mat> * _free_mask, int _nfree);

#endif
