#ifndef _OPTIMIZE_STATE_RECONSTRUCTOR_GSL_H_
#define _OPTIMIZE_STATE_RECONSTRUCTOR_GSL_H_

class RateModel; // forward declaration
class StateReconstructor; // forward declaration

#include <armadillo>
using namespace arma;

#include <gsl/gsl_vector.h>

class OptimizeStateReconstructor {
private:
    RateModel * rm;
    StateReconstructor * sr;
    mat * free_variables;
    int nfree_variables;
    int maxiterations;
    double stoppingprecision;

    double GetLikelihoodWithOptimized(const gsl_vector * variables);
    static double GetLikelihoodWithOptimized_gsl(const gsl_vector * variables, void *obj);

public:
    OptimizeStateReconstructor(RateModel *,StateReconstructor *, mat *, int);
    mat optimize();
};

#endif /* _OPTIMIZE_STATE_RECONSTRUCTOR_GSL_H_ */
