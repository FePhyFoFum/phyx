#ifndef _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_
#define _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_

#include <nlopt.h>

class RateModel; // forward declaration
class StateReconstructor; // forward declaration

#include <armadillo>
using namespace arma;

void optimize_sr_nlopt(RateModel * _rm,StateReconstructor * _sr, mat * _free_mask, int _nfree);

#endif /* _OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H_ */
