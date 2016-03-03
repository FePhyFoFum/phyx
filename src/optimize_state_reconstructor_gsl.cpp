
#include <math.h>
#include <vector>
#include <limits>
using namespace std;

#include <armadillo>
using namespace arma;

#include "optimize_state_reconstructor.h"
#include "rate_model.h"
#include "state_reconstructor.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

OptimizeStateReconstructor::OptimizeStateReconstructor(RateModel * _rm,StateReconstructor * _sr, mat * _free_mask, int _nfree):
    rm(_rm),sr(_sr),free_variables(_free_mask),nfree_variables(_nfree),maxiterations(10000),stoppingprecision(0.001) {}

double OptimizeStateReconstructor::GetLikelihoodWithOptimized(const gsl_vector * variables) {
    for (unsigned int i=0; i < free_variables->n_rows; i++) {
        for (unsigned int j=0; j < free_variables->n_cols; j++) {
            if (i != j) {
                rm->set_Q_cell(i,j,gsl_vector_get(variables,(*free_variables)(i, j)));
                if (rm->get_Q()(i, j) < 0 || rm->get_Q()(i, j) >= 1000) {
                    return 100000000;
                }
            }
        }
    }
    double like;
    rm->set_Q_diag();
    like = sr->eval_likelihood();
    //cout << like << endl;
    if (like < 0 || like == std::numeric_limits<double>::infinity()) {
        like = 100000000;
    }
    return like;
}


double OptimizeStateReconstructor::GetLikelihoodWithOptimized_gsl(const gsl_vector * variables, void *obj) {
    double temp;
    temp= ((OptimizeStateReconstructor*)obj)->GetLikelihoodWithOptimized(variables);
    return temp;
}

/*
 * USES THE SIMPLEX ALGORITHM
 *
 */
mat OptimizeStateReconstructor::optimize() {
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    size_t np = nfree_variables;
    size_t iter = 0, i;
    int status;
    double size;
    /* Initial vertex size vector */
    ss = gsl_vector_alloc (np);
    /* Set all step sizes to .01 */ //Note that it was originally 1
    gsl_vector_set_all (ss, .2);
    /* Starting point */
    //cout<<"Now in OPtimizaRateWithGivenTipVariance in OptimizationFn"<<endl;
    x = gsl_vector_alloc (np);
    for (unsigned int i=0; i < np; i++) {
        gsl_vector_set (x, i, 0.1);
    }
    OptimizeStateReconstructor *pt;
    pt=(this);
    double (*F)(const gsl_vector *, void *);
    F = &OptimizeStateReconstructor::GetLikelihoodWithOptimized_gsl;
    /* Initialize method and iterate */
    gsl_multimin_function minex_func;
    minex_func.f =* F;
    minex_func.params = pt;
    minex_func.n = np;
    s = gsl_multimin_fminimizer_alloc (T, np);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        //cout<<"Now on iteration "<<iter<<endl;
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status != 0) { //0 Means it's a success
        //    printf ("error: %s\n", gsl_strerror (status));
            break;
        }
        size = gsl_multimin_fminimizer_size (s);
        //status = gsl_multimin_test_size (size, 1e-2);
        status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
        if (status == GSL_SUCCESS) {
        //    printf ("converged to minimum at\n");
        }
        //printf ("%5d ", iter);
        for (i = 0; i < np; i++) {
        //    printf ("%10.3e ", gsl_vector_get (s->x, i));
        }
        //printf ("f() = %7.3f size = %.3f\n", s->fval, size);
    }
    while (status == GSL_CONTINUE && iter < maxiterations);
    mat results (free_variables->n_rows, free_variables->n_cols); results.fill(0);
    for (unsigned int i=0; i < results.n_rows; i++) {
        for (unsigned int j=0; j < results.n_cols; j++) {
            if (i != j) {
                results(i, j) = (gsl_vector_get(s->x, (*free_variables)(i, j)));
            }
        }
    }
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return results;
}
