// a few handy routines that aren't in the mpfr distribution

#include "mpfr.h"

#ifdef __cplusplus
extern "C" {
#endif

int mpfr_mul_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t);
int mpfr_div_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t);
int mpfr_d_div(mpfr_ptr z, double a, mpfr_srcptr x, mpfr_rnd_t);
int mpfr_add_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t);
int mpfr_sub_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t);
int mpfr_d_sub(mpfr_ptr z, double a, mpfr_srcptr x, mpfr_rnd_t);
void mpfr_mul_d_clear();

#ifdef __cplusplus
} // end extern "C" block
#endif
