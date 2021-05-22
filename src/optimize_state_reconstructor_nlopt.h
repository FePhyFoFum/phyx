#ifndef PX__OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H
#define PX__OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H

#include <nlopt.h>

class RateModel; // forward declaration
class StateReconstructor; // forward declaration

#include <armadillo>
using namespace arma;

void optimize_sr_nlopt(RateModel * _rm,StateReconstructor * _sr, mat * _free_mask, int _nfree);

#endif /* PX__OPTIMIZE_STATE_RECONSTRUCTOR_NLOPT_H */
