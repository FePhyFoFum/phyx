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
 */

#ifndef SUPERDOUBLE_H
#define SUPERDOUBLE_H
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;


class Superdouble  {
private:
	long double mantissa;
	int exponent;
	bool stilldouble;
	double upperlimit;
	double lowerlimit;
	void adjustDecimal();
	friend ostream& operator<<(ostream& os, const Superdouble& x);
	
public:
	Superdouble(long double mantissa=1.0, int exponent=0);
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
	bool operator < (const Superdouble &x)const ;
	bool operator > (const Superdouble &x)const ;
	bool operator >= (const Superdouble &x)const ;
	bool operator <= (const Superdouble &x)const ;
	int getExponent();
	double getMantissa();
	Superdouble getLn();
	Superdouble abs();
	void switch_sign();
	
	operator double() {return mantissa*pow(10.,exponent);};
	
};
#endif
