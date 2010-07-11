// this file implements the routines required by gmpfrxx.h
// it is modeled on mpfr-2.2.1/out_str.c
//                   gmp-4.2.1/mpz/set_f.c
//                   gmp-4.2.1/mpq/set_f.c

// (written by Jon Wilkening on Aug 8, 2007)


#include "gmpfrxx.h"
#include <cstdlib>
#include <climits>

using namespace std;


// MpFrC is a dummy class holding the static members of mpfr_class
mpfr_rnd_t  MpFrC::rnd = GMP_RNDN;
int         MpFrC::base = 10;


istream & operator>>(istream &s, mpfr_ptr a) {
    string tmp;
    s >> tmp;
    mpfr_set_str(a, tmp.c_str(),
		 mpfr_class::get_base(),
		 mpfr_class::get_rnd());
    // a = tmp.c_str();
    return s;
}

// modeled on mpfr-2.2.1/out_str.c
ostream & operator<<(ostream &os, mpfr_srcptr a)
{
    char *s, *t, *s0;
    mp_exp_t  e;

    // for debugging:
    // mpfr_out_str(stdout, 10, 0, a, RND); printf("\n");

    if (mpfr_nan_p(a)) {
	os << "@NaN@";
	return os;
    }

    if (mpfr_inf_p(a)) {
	if (MPFR_SIGN(a) > 0)
	    os << "@Inf@";
	else
	    os << "-@Inf@";
	return os;
    }

    if (mpfr_zero_p(a)) {
	if (MPFR_SIGN(a) > 0)
	    os << "0";
	else
	    os << "-0";
	return os;
    }

    s = mpfr_get_str (NULL, &e,
		      mpfr_class::get_base(),
		      0, a,
		      mpfr_class::get_rnd());

    t = mpfr_get_str (NULL, &e,
		      mpfr_class::get_base(),
		      os.precision(), a,
		      mpfr_class::get_rnd());

    if (strlen(s)<strlen(t))
	mpfr_free_str(t);
    else {
	mpfr_free_str(s);
	s = t;
    }

    s0 = s;
    /* for a=3.1416 we have s = "31416" and e = 1 */

    if (*s == '-')
	os.put(*s++);

    /* outputs mantissa */
    os.put(*s++); e--; /* leading digit */
    os.put('.');
    while (*s != '\0')
	os.put(*s++);       /* rest of mantissa */

    mpfr_free_str(s0);

    /* outputs exponent */
    if (e) {
	os << (mpfr_class::get_base() <= 10 ? 'e' : '@') << (long) e;
    }

    return os;
}

// fixme: it would be more efficient to work directly
// with the limbs, but I don't want to deal with the
// internal representations of both mpfr and gmp
void mpz_set_mpfr(mpz_ptr w, mpfr_srcptr u)
{
    char *s, *t;
    long k, sz;
    mp_exp_t  e;
  
    // abs(u)<1 truncates to zero
    if ( mpfr_get_exp(u)<=0 || mpfr_zero_p(u)
	 || mpfr_nan_p(u) || mpfr_inf_p(u) ) {
	mpz_set_ui(w, 0);
	return;
    }

    // note: this is done in hex to represent u exactly
    // example: for u=3.1416 we have s = "31416" and e = 1
    s = mpfr_get_str (NULL, &e, 16, 0, u, GMP_RNDN);
    sz = strlen(s);

    if (*s == '-')
	e++;

    if (e<=sz) {
	s[e] = '\0'; // truncate
	mpz_set_str(w, s, 16);
	mpfr_free_str(s);
	return;
    }

    t = (char *) calloc(e+1, sizeof(char));
    
    for (k=0; k<sz; k++)
	t[k] = s[k];
    for (k=sz; k<e; k++)
	t[k] = '0';
    t[e] = '\0';

    mpz_set_str(w, t, 16);
    free(t);
    mpfr_free_str(s);
}


void mpq_set_mpfr(mpq_ptr w, mpfr_srcptr u)
{
    char *s, *t;
    long k, sz;
    mp_exp_t  e, bits;
    
    if ( mpfr_zero_p(u) || mpfr_nan_p(u) || mpfr_inf_p(u) ) {
	mpq_set_ui(w, 0, 1);
	return;
    }

    // note: this is done in binary to represent u exactly
    // example:   u=10.1001  -->  s = "101001", e = 2
    s = mpfr_get_str (NULL, &e, 2, 0, u, GMP_RNDN);
    sz = strlen(s);

    t = s+sz;
    while (*(--t) == '0')
	*t = '\0'; // trim trailing zeros

    sz = strlen(s);
    bits = (*s == '-') ? sz-1 : sz; // bits>0 since u!=0

    if (e>=bits) {
	e = sz+e-bits;
	// no denominator needed
	t = (char *) malloc((e+1)*sizeof(char));
	strcpy(t,s);
	for (k=sz; k<e; k++)
	    t[k] = '0';
	t[e] = '\0';
	mpq_set_str(w, t, 2);
	free(t);
	mpfr_free_str(s);
	return;
    }

    e = sz+bits-e+2;
    t = (char *) malloc((e+1)*sizeof(char));
    strcpy(t,s);
    t[sz++] = '/';
    t[sz++] = '1';
    for (k=sz; k<e; k++)
	t[k] = '0';
    t[e] = '\0';
    mpq_set_str(w, t, 2);
    free(t);
    mpfr_free_str(s);
}
