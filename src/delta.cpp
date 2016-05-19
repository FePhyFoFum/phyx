
#include <math.h>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <iostream>

using namespace std;

#include "gmpfrxx/gmpfrxx.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>

#include "delta.h"

#define LARGE 100000000
#define PMAX 1000

Delta::Delta() {}

double Delta::shift(double p, double s, int c, int l, int r) {
    if (p < 0 || s < 0 || p > PMAX || s > PMAX)
        return LARGE;
    double enp = exp(-p);
    double ens = exp(-s);
    double num = (enp*pow((1.0 - enp),(l - 1.0)))*(ens*pow((1.0 - ens),(r - 1.0)));
    double den = 0;
    for (int i=0; i < (l+r-1); i++) {
        int il = i + 1;
        den += (enp*pow((1 - enp),(il - 1)))*(ens*pow((1 - ens),((l + r) - il - 1)));
    }
    double bigleft = 0;
    try {
        bigleft = log(num/den);
    } catch( char * str ) {
        return LARGE;
    }
    double enc = exp(-c);
    num = (enc*pow((1.0 - enc),(l - 1.0)))*(enc*pow((1.0 - enc),(r - 1.0)));
    den = 0;
    for (int i=0; i < (l+r-1); i++) {
        int il = i + 1;
        den = den +(enc*pow((1 - enc),(il - 1)))*(enc*pow((1 - enc),((l + r) - il - 1)));
    }

    double bigright = 0;
    try {
        bigright = log(num/den);
    } catch( char * str ) {
        return LARGE;
    }
    //cout << "bl:" << bigleft << "\tbr:"<< bigright << endl;
    return -(bigleft-bigright);
}

/*
 * use for over/underflow problems
 */

double Delta::bigshift(double p, double s, int c, int l, int r) {
    if (p < 0 || s < 0 || p > PMAX || s > PMAX)
            return LARGE;
    mpfr_class bp,bs,bc,bl,br;
    bp = p;bs = s;bc = c;bl = l;br = r;
    mpfr_class expnbp = exp(-bp);
    mpfr_class expnbs = exp(-bs);
    mpfr_class num = (expnbp*pow((1.0 - expnbp),(bl - 1.0)))*(expnbs*pow((1.0 - expnbs),(br - 1.0)));
    mpfr_class den = 0;
    for (mpfr_class i=0; i < (l+r-1); i++) {
        mpfr_class il=i+1;
        den = den +(expnbp*pow((1 - expnbp),(il - 1)))*(expnbs*pow((1 - expnbs),((bl + br) - il - 1)));
    }
    mpfr_class bigleft = 0;
    try{
        bigleft = log(num/den);
    }catch( char * str ) {
        return LARGE;
    }
    mpfr_class expnbc = exp(-bc);
    num = (expnbc*pow((1.0 - expnbc),(bl - 1.0)))*(expnbc*pow((1.0 - expnbc),(br - 1.0)));
    den = 0;
    for (mpfr_class i=0; i < (l+r-1); i++) {
        mpfr_class il=i+1;
        den = den +(expnbc*pow((1 - expnbc),(il - 1)))*(expnbc*pow((1 - expnbc),((bl + br) - il - 1)));
    }

    mpfr_class bigright = 0;
    try {
        bigright = log(num/den);
    } catch( char * str ) {
        return LARGE;
    }
    //cout << "bl:" << bigleft << "\tbr:"<< bigright << endl;
    mpfr_class f = -(bigleft-bigright);
    double x = f.get_d();
    return x;
}

double Delta::cdf(double del) {
    double pvalue = 1-gsl_cdf_gaussian_P(del, 1.3);
    return pvalue;
}

vector<double> Delta::delta(int l, int r,int o) {
    OptimizeShift os(this);
    os.setCLR(1,l+r,o);
    vector<double> resout = os.optimize_shift();
    double out;
    if (l > 100 || r > 100 || o > 100) {
        out =  bigshift(resout[0],resout[1],1,l+r,o);
    } else {
        out =  shift(resout[0],resout[1],1,l+r,o);
    }
    os.setCLR(1,l,r);
    resout = os.optimize_shift();
    double in;
    if (l > 100 || r > 100 || o > 100) {
        in = bigshift(resout[0],resout[1],1,l,r);
    } else {
        in = shift(resout[0],resout[1],1,l,r);
    }
    vector<double> ret;
    ret.push_back(-(out-in));
    ret.push_back(cdf(-(out-in)));
    return ret;
}


OptimizeShift::OptimizeShift(Delta * delt):delta(delt),maxiterations(10000),stoppingprecision(0.001),c(0),l(0),r(0) {}

void OptimizeShift::setCLR(int ce, int le, int ri) {
    c = ce; l = le; r = ri;
}

double OptimizeShift::GetShift(const gsl_vector * variables) {
    double p=gsl_vector_get(variables,0);
    double s=gsl_vector_get(variables,1);
    double shift;
    if (c > 100 || l > 100 || r > 100) {
        shift = delta->bigshift(p,s,c,l,r);
    } else {
        shift = delta->shift(p,s,c,l,r);
    }
    if (shift == std::numeric_limits<double>::infinity() || isnan(shift)) {
        shift = 100000000;
    }
    return shift;
}


double OptimizeShift::GetShift_gsl(const gsl_vector * variables, void *obj) {
    double temp;
    temp = ((OptimizeShift*)obj)->GetShift(variables);
    return temp;
}

vector<double> OptimizeShift::optimize_shift() {
    //need to check the performance on this
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    size_t np = 2;
    size_t iter = 0;
    int status;
    double size;
    /* Initial vertex size vector */
    ss = gsl_vector_alloc (np);
    /* Set all step sizes to .01 */ //Note that it was originally 1
    gsl_vector_set_all (ss, 0.1);
    /* Starting point */
    //cout<<"Now in OPtimizaRateWithGivenTipVariance in OptimizationFn"<<endl;
    x = gsl_vector_alloc (np);
    gsl_vector_set (x,0,0.5);
    gsl_vector_set (x,1,0.5);
    OptimizeShift *pt;
    pt=(this);
    double (*F)(const gsl_vector *, void *);
    F = &OptimizeShift::GetShift_gsl;
    /* Initialize method and iterate */
    gsl_multimin_function minex_func;
    minex_func.f=*F;
    minex_func.params=pt;
    minex_func.n = np;
    s = gsl_multimin_fminimizer_alloc (T, np);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status != 0) { //0 Means it's a success
            printf ("error: %s\n", gsl_strerror (status));
            break;
        }
        size = gsl_multimin_fminimizer_size (s);
        //status = gsl_multimin_test_size (size, 1e-2);
        status = gsl_multimin_test_size (size, stoppingprecision); //since we want more precision
    }
    while (status == GSL_CONTINUE && iter < maxiterations);
    //cout << "iterations:" << iter << endl;
    vector<double> results;
    results.push_back(gsl_vector_get(s->x,0));
    results.push_back(gsl_vector_get(s->x,1));
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return results;
}
