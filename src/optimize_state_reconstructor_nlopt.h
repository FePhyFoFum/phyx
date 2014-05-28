/*
 * optimize_tnc.h
 *
 *  Created on: Feb 9, 2010
 *      Author: smitty
 */

#ifndef _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_
#define _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_

#include <nlopt.h>

#include "state_reconstructor.h"
#include "rate_model.h"

#include <armadillo>
using namespace arma;

void optimize_sr_nlopt(RateModel * _rm,StateReconstructor * _sr, mat * _free_mask, int _nfree);

#endif /* _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_ */
