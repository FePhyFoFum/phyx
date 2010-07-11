// example file to illustrate gmpfrxx.h

#include "gmpfrxx.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <complex>

using namespace std;


void example1() {
    mpfr_class f;
    mpq_class  q;

    f.set_prec(128); // precision in bits
    f = "1.7e12";
    
    q = "-5/3";

    f = 12.5*(f*q + mpfr_class("123512351235")/7);

    cout.precision(60); // max significant digits to print

    cout << "example 1\n";
    cout << q << endl;
    cout << f << endl;

    f.prec_round(16);
    cout << f << endl;
    
    f.set_base(16);  // output in base 16
    cout << f << endl;
    mpfr_class::set_base(10); // change it back to base 10

    ostringstream s;
    s << f;
    cout << s.str() << endl;
}


void example2() {
    mpz_class  z("10000000000000000005");
    mpfr_class f(z,64);

    f = -f/7000; // roundoff error occurs
    mpq_class  q(f); // get f as an exact rational number
    mpfr_class g(q,64); // convert q back to floating point

    cout << "\nexample 2\n";
    cout << f << endl;
    cout << q << endl;
    cout << f - g << endl;
}


void example3() {
    mpfr_class f(-4.75, 5);  // 5 bit rep of 4.75
    mpfr_class  g(0.0, 32);  // 32 bits

    cout << "\nexample 3\n";
    // f.set_base(2);
    cout << "f = " << f << endl;

    mpfr_class::set_rnd(GMP_RNDD);
    g = rint(f);
    cout << "GMP_RNDD    : " << g << endl;
    mpfr_class::set_rnd(GMP_RNDU);
    g = rint(f);
    cout << "GMP_RNDU    : " << g << endl;
    mpfr_class::set_rnd(GMP_RNDZ);
    g = rint(f);
    cout << "GMP_RNDZ    : " << g << endl;
    mpfr_class::set_rnd(GMP_RND_MAX);
    g = rint(f);
    cout << "GMP_RND_MAX : " << g << endl;
    mpfr_class::set_rnd();
    g = rint(f);
    cout << "GMP_RNDN    : " << g << endl;

    mpz_class z(f); // truncates regardless of rounding mode
    cout << "z = " << z << endl;
}


void example4() {
    int n = 10;
    gmp_randclass  w(gmp_randinit_default);
    vector<mpfr_class> g( n, mpfr_class(0.0,4) );

    w.seed(1);

    cout << "\nexample 4\n";
    cout << n << " random numbers, each 4 bits\n"
	 << "uniformly distributed between 0<=f<1\n"
	 << "(sorted after generating)\n";

    for (int i=0; i<n; i++)
	g[i] = w.get_f();

    sort(g.begin(), g.end());

    for (int i=0; i<n; i++)
	cout << g[i] << endl;
}


void example5() {
    mpfr_class t(1.0, 10); // 10 bits, t=1
    mpfr_class s, c, S, C; // 53 bits
    s.set_prec(100); // 100 bits now
    S.set_prec(100);

    // binary left shift multiplies t by 2^p (p unsigned long)
    t = t<<3; // t=8 now

    cout << "\nexample 5\n";
    cout << "t = " << t << endl;

    // direct call to mpfr routine.
    // (a bit more efficient if sin, cos are both needed)
    mpfr_sin_cos(s.get_mpfr_t(), c.get_mpfr_t(), t.get_mpfr_t(), GMP_RNDN);
    cout << "sin(t), cos(t) = " << s << "  " << c << endl;

    // template version of the same
    // note that target still determines precision of result
    S = sin(t);
    C = cos(t);
    cout << "sin(t), cos(t) = " << S << "  " << C << endl;

    // even here the target determines the precision of all
    // intermediate calculations
    S = sin(t)*sin(t) + cos(t)*cos(t);
    cout << S << endl;

    S = sin(const_pi());
    cout << S << endl;
}

void example6() {
    mpfr_class x(1.0,53);
    mpfr_class y(1.0,106);

    y = const_pi();
    complex<mpfr_class> z(x,y);

    cout << "\nexample 6\n";
    cout.precision(0);
    cout << z << endl;
    // cout << exp(z) << endl;
}

int main() {
    example1();
    example2();
    example3();
    example4();
    example5();
    example6();

    return 0;
}
