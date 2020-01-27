#ifndef _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_
#define _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_

#include <vector>
#include <nlopt.h>

class RateModel; // forward declaration
class StateReconstructor; // forward declaration

#include <armadillo>
using namespace arma;

void optimize_sr_periods_nlopt(std::vector<RateModel> * _rm,StateReconstructor * _sr,
    std::vector<mat> * _free_mask, int _nfree);

#endif /* _OPTIMIZE_STATE_RECONSTRUCTOR_PERIODS_NLOPT_H_ */
