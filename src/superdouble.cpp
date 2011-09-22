
#include "superdouble.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;


Superdouble::Superdouble(long double m, int e):stilldouble(false),upperlimit(1e+100),lowerlimit(1e-100){
	mantissa=m;
	exponent=e;
	if (exponent == 0){
		stilldouble = true;
	}else{
		adjustDecimal();
	}
}

Superdouble::~Superdouble() {}

int Superdouble::getExponent(){
	return exponent;
}

double Superdouble::getMantissa(){
	return mantissa;
}

void Superdouble::adjustDecimal() {
	stilldouble = false;
	if (mantissa==0 || isinf(mantissa) || isnan(mantissa)) {
		exponent=0;
		stilldouble = true;
	}
	else {
		while (fabs(mantissa)>=10) {
			mantissa*=0.1;
			exponent+=1;
		}
		while (fabs(mantissa)<1) {
			mantissa*=10.0;
			exponent+=-1;
		}
	}
}

ostream& operator<<(ostream& os, const Superdouble& x)
{
	os<<x.mantissa <<"e"<<x.exponent;
	return os;
}

Superdouble Superdouble::operator * ( Superdouble  x){
	Superdouble result(mantissa*x.mantissa,exponent+x.exponent);
	if (result.stilldouble == true){
		if (fabs(result.mantissa) > upperlimit || fabs(result.mantissa) < lowerlimit){
			result.adjustDecimal();
		}
	}else{
		result.adjustDecimal();
	}
	return result;
}

Superdouble Superdouble::operator * ( double  x){
	Superdouble result(mantissa*x,exponent);
	if (result.stilldouble == true){
		if (fabs(result.mantissa) > upperlimit || fabs(result.mantissa) < lowerlimit){
			result.adjustDecimal();
		}
	}else{
		result.adjustDecimal();
	}
	return result;
}


//add stilldouble
Superdouble Superdouble::operator / ( Superdouble  x){
	Superdouble result(mantissa/x.mantissa,exponent-x.exponent);
	result.adjustDecimal();
	return result;
}

//add stilldouble
Superdouble Superdouble::operator + ( Superdouble  x){

	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		Superdouble result(mantissa+(x.mantissa*(pow(10,exponentdif))),exponent);
		result.adjustDecimal();
		return result;
	}
	else {
		Superdouble result(x.mantissa,x.exponent);
		result.adjustDecimal();
		return result;
	}
}

//add stilldouble
Superdouble Superdouble::operator - ( Superdouble  x){
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		Superdouble result(mantissa-(x.mantissa*(pow(10,exponentdif))),exponent);
		result.adjustDecimal();
		return result;
	}
	else {
		Superdouble result(-1.0*x.mantissa,x.exponent);
		result.adjustDecimal();
		return result;
	}
}

//add stilldouble
void Superdouble::operator ++ (){
	mantissa++;
	adjustDecimal();
}

//add stilldouble
void Superdouble::operator -- (){
	mantissa--;
	adjustDecimal();
}

//add stilldouble
void Superdouble::operator *= (Superdouble x){
	mantissa*=x.mantissa;
	exponent+=x.exponent;
	adjustDecimal();
}

//add stilldouble
void Superdouble::operator /= (Superdouble x){
	mantissa/=x.mantissa;
	exponent-=x.exponent;
	adjustDecimal();
}

void Superdouble::operator += ( Superdouble  x){
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		if (stilldouble == true && x.stilldouble == true){
			mantissa += x.mantissa;
			if (fabs(mantissa) > upperlimit || fabs(mantissa) < lowerlimit){
				adjustDecimal();
			}
		}else{
			int exponentdif=x.exponent-exponent;
			mantissa=mantissa+(x.mantissa*(pow(10,exponentdif)));
			adjustDecimal();
		}
	}else {
		if (stilldouble == true && x.stilldouble == true){
			mantissa = x.mantissa;
			exponent=x.exponent;
		}else{
			mantissa=x.mantissa;
			exponent=x.exponent;
			adjustDecimal();
		}
	}
}

//add stilldouble
void Superdouble::operator -= ( Superdouble  x){
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		mantissa=mantissa-(x.mantissa*(pow(10,exponentdif)));
		adjustDecimal();
	}
	else {
		mantissa=-1.0*x.mantissa;
		exponent=x.exponent;
		adjustDecimal();
	}
}

bool Superdouble::operator > (const Superdouble & x)const{
	if (exponent > x.exponent)
		return true;
	else if(exponent == x.exponent && mantissa > x.mantissa)
		return true;
	else
		return false;
}

bool Superdouble::operator >= (const Superdouble & x)const{
	if (exponent > x.exponent)
		return true;
	else if(exponent == x.exponent && mantissa >= x.mantissa)
		return true;
	else
		return false;
}

bool Superdouble::operator < (const Superdouble & x)const{
	if (exponent < x.exponent)
		return true;
	else if(exponent == x.exponent && mantissa < x.mantissa)
		return true;
	else
		return false;
}

bool Superdouble::operator <= (const Superdouble & x)const{
	if (exponent < x.exponent)
		return true;
	else if(exponent == x.exponent && mantissa <= x.mantissa)
		return true;
	else
		return false;
}

//this just switches the sign of the superdouble
void Superdouble::switch_sign(){
	mantissa = -1*mantissa;
}

/*bool Superdouble::operator > (double x){
	if (double() > x)
		return true;
	else
		return false;
}*/

Superdouble Superdouble::getLn(){
	//ln(a * 10^b) = ln(a) + ln(10^b) = ln(a) + log10 (10^b) / log10 (e^1) = ln(a) + b/log10(e^1)
	Superdouble result(log(mantissa)+(1.0*(exponent))/log10(exp(1)),0);
	result.adjustDecimal();
	return result;
}

Superdouble Superdouble::abs(){
	if (mantissa < 0){
		Superdouble result(-mantissa,exponent);
		return result;
	}else{
		Superdouble result(mantissa,exponent);
		return result;
	}

}

