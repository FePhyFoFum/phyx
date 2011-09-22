#include "mpfr_mul_d.h"

static mpfr_t y;
static int initialized = 0;

int mpfr_mul_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_mul(z, x, y, r);
}

int mpfr_div_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_div(z, x, y, r);
}

int mpfr_d_div(mpfr_ptr z, double a, mpfr_srcptr x, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_div(z, y, x, r);
}

int mpfr_add_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_add(z, x, y, r);
}

int mpfr_sub_d(mpfr_ptr z, mpfr_srcptr x, double a, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_sub(z, x, y, r);
}

int mpfr_d_sub(mpfr_ptr z, double a, mpfr_srcptr x, mpfr_rnd_t r)
{
    if (!initialized) {
	mpfr_init2(y, 53);
	initialized = 1;
    }
    mpfr_set_d(y, a, GMP_RNDN);
    return mpfr_sub(z, y, x, r);
}

void mpfr_mul_d_clear()
{
    if (initialized) {
	initialized = 0;
	mpfr_clear(y);
    }
}
