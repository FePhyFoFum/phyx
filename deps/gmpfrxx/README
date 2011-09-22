Copyright 1991, 1996, 1999, 2000, 2007 Free Software Foundation, Inc.

gmpfrxx is an extension the GNU MP Library and the MPFR Library.

gmpfrxx is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

gmpfrxx is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with gmpfrxx.  If not, see http://www.gnu.org/licenses/.


----  about gmpfrxx  ----

This code allows one to freely combine integer, rational and
floating point objects in arbitrarily complex algebraic expressions.
When possible, it avoids creating temporary variables and directly
calls the appropriate gmp or mpfr function.

To learn to use gmpfrxx.h, start by looking over the example file
and try compiling/running it.  Then consult the gmp documentation
(i.e. gmp-man-4.2.1.pdf, from the gmp webpage): everywhere mpf_class
is mentioned, the mpfr_class will behave similarly, calling the
appropriate mpfr function instead (see mpfr.pdf from the mpfr webpage).
Note: you cannot include both gmpfrxx.h and gmpxx.h in the same file.

Rounding is controlled via a static member function in mpfr_class.
rounding modes: {GMP_RNDN, GMP_RNDZ, GMP_RNDU, GMP_RNDD, GMP_RND_MAX}
static member functions controlling rounding:
mpfr_rnd_t  mpfr_class::get_rnd()
void        mpfr_class::set_rnd(mpfr_rnd_t r=GMP_RNDN)

We include a static member function to set the base for
subsequent calls of cin, cout, or any stream or string i/o.
int  mpfr_class::get_base()
void mpfr_class::set_base(base=10), 2<=base<=36

We include a static member to set the default rounding precision:
mpfr_prec_t  mpfr_class::get_dprec()
void         mpfr_class::set_dprec(mpfr_prec_t  p=53)

We added the following mpfr routines not in gmpxx:
  prec_round, rint, log, log2, log10, exp, exp2, exp10,
  cos, sin, tan, sec, csc, cot, acos, asin, atan, atan2,
  cosh, sinh, tanh, sech, csch, coth, acosh, asinh, atanh,
  fac_ui, log1p, expm1, eint, gamma, lngamma, zeta, erf,
  erfc, agm, const_log2, const_pi, const_euler, const_catalan

on 11/16/08 we added
  sqr, cbrt, root, pow, pow_ui, fac_ui, max, min, dim, cmpabs,
  lgamma, zeta_ui, j0, j1, jn, y0, y1, yn, frac, remainder

please e-mail Jon Wilkening <wilken at math berkeley edu>
to report bugs or make suggestions (especially if the gmpxx
behavior is different than what my code does)

last modified: 11/16/2008
tested with  mpfr-2.3.1, gmp-4.2.2
(should work with other versions as well.
 no changes were needed for mpfr-2.3.2)

---

to compile the example file:

(1) install or update mpfr and gmp.  see notes below.
I have tested the code with  mpfr-2.3.1,  gmp-4.2.2


(2) edit the Makefile and add paths to the libraries (if they
aren't in the usual library search path)

(3) make
(4) ./example

--
installing mpfr and gmp:

I found the installation procedure for both codes to be straightforward.
The notes below were for mpfr-2.2.1, which was current at the time of
the first release of gmpfrxx.h

---
for gmp, I unpacked the source files and installed via

tar xvzf gmp-4.2.1.tar.gz
cd gmp-4.2.1
./configure --prefix=/home/wilken --enable-cxx
make
make check
make install

this creates header files in /home/wilken/include and libraries
in /home/wilken/lib

---
for mpfr, I unpacked the source and installed via

tar xvzf mpfr-2.2.1.tar.gz
cd mpfr-2.2.1
wget http://www.mpfr.org/mpfr-2.2.1/patches
patch -N -Z -p1 < patches
./configure --prefix=/home/wilken
make
make check
make install

this creates header files in /home/wilken/include and libraries
in /home/wilken/lib
