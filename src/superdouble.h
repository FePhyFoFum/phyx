/*
 * superdouble, a class with long double precision but less subject to overflow or underflow
 * superdouble X=mantissa * 10^exponent
 *
 * Copyright Brian C. O'Meara, Oct. 2, 2008
 * http://www.brianomeara.info
 *
 * Additional developers wanted!!!
 *
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 more development by Stephen A. Smith
 modernizing by Joseph W. Brown (2021)
 */

#ifndef PX_SUPERDOUBLE_H
#define PX_SUPERDOUBLE_H

#include <cmath>
#include <iostream>

class Superdouble {
private:
    long double mantissa;
    int exponent;
    bool stilldouble;
    double upperlimit;
    double lowerlimit;
    void adjustDecimal ();
    friend std::ostream& operator << (std::ostream& os, const Superdouble& x);
    
public:
    Superdouble(long double m=1.0l, int e=0);
    ~Superdouble();
    Superdouble operator* ( Superdouble x);
    Superdouble operator* ( double x);
    Superdouble operator/ ( Superdouble x);
    Superdouble operator+ ( Superdouble x);
    Superdouble operator- ( Superdouble x);
    void operator++ ();
    void operator -- ();
    void operator*= (Superdouble x);
    void operator/= (Superdouble x);
    void operator+= (Superdouble x);
    void operator-= (Superdouble x);
    bool operator < (const Superdouble& x) const;
    bool operator > (const Superdouble& x) const;
    bool operator >= (const Superdouble& x) const;
    bool operator <= (const Superdouble& x) const;
    int getExponent () const;
    double getMantissa ();
    Superdouble getLn ();
    Superdouble abs ();
    void switch_sign ();
    
    operator double() {
        return static_cast<double> (mantissa) * pow(10., exponent);
    }
};
#endif /* PX_SUPERDOUBLE_H */
