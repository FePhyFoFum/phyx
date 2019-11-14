/*
 * optimize_tnc.h
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#ifndef _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_
#define _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_

#include <vector>
#include <nlopt.h>

#include "state_reconstructor.h"
#include "rate_model.h"

#include <armadillo>
using namespace arma;

void optimize_sr_periods_nlopt(std::vector<RateModel> * _rm,StateReconstructor * _sr,
    std::vector<mat> * _free_mask, int _nfree);

#endif /* _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_ */
