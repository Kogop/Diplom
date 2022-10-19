#pragma once
#include <iostream>
#include <math.h>

using namespace std;
int real, imag;
const double pi = 4.0 * atan(1.0);

class complex1
{
public:
	double q;
	double real, imag;
	int n;
	complex1 Sum(complex1 z1, complex1 z2);
	complex1 Dif(complex1 z1, complex1 z2);
	complex1(double r1, double i1) {
		real = r1;
		imag = i1;
	}
	

	complex1() {
		imag = 0;
	}
	
	double abs(complex1 z1) {
		return sqrt(z1.real * z1.real + z1.imag * z1.imag);
	}
	double arg() {
		double f = 0;
		if (real > 0) { f = atan((double)imag / real); }
		if ((real < 0) && (imag >= 0)) { f = 3.14 + atan((double)imag / real); }
		if ((real < 0) && (imag < 0)) { f = -3.14 + atan((double)imag / real); }
		if ((real == 0) && (imag > 0)) { f = (double)3.14 / 2; }
		if ((real == 0) && (imag < 0)) { f = (double)-3.14 / 2; }
		return f;
	}
	friend void printcomplex1(complex1 z);
	complex1 Umnoj(complex1 z1, complex1 z2);
	complex1 Delen(complex1 z1, complex1 z2);
	complex1 Deg(complex1 z1, complex1 z2);
	complex1 Exp(complex1 z);
	complex1 Ln(complex1 z);
	complex1 Sin(complex1 z);
	complex1 Cos(complex1 z);

	complex1 Deg1(complex1 z1, complex1 z2);
	complex1 Root(complex1 z);

	friend complex1 operator*(const complex1& z1, double q);
	friend complex1 operator*(double q ,const complex1& z1);
	friend complex1 operator*(const complex1& z1, const complex1& z2);
	friend complex1 operator+(const complex1& z1, double q);
	friend complex1 operator+(double q, const complex1& z1);
	friend complex1 operator+(const complex1& z1, const complex1& z2);
	friend complex1 operator/(const complex1& z1, double q);
	friend complex1 operator/(double q, const complex1& z1);
	friend complex1 operator/(const complex1& z1, const complex1& z2);
	friend complex1 operator-(const complex1& z1, double q);
	friend complex1 operator-(double q, const complex1& z1);
	friend complex1 operator-(const complex1& z1, const complex1& z2);
};

const complex1 _i = complex1 (0, 1);



void printcomplex1(complex1 z) {
	cout << z.real << ",  " << z.imag << endl;
}

complex1 complex1::Sum(complex1 z1, complex1 z2) {

	return complex1(z1.real + z2.real, z1.imag + z2.imag);
}

complex1 complex1::Dif(complex1 z1, complex1 z2) {

	return complex1(z1.real - z2.real, z1.imag - z2.imag);
}
complex1 complex1::Umnoj(complex1 z1, complex1 z2) {

	return complex1(z1.real * z2.real - z1.imag * z2.imag, z1.real * z2.imag + z2.real * z1.imag);	
}

complex1 complex1::Delen(complex1 z1, complex1 z2) {

	return complex1(((z1.real * z2.real + z1.imag * z2.imag) / (z2.real * z2.real + z2.imag * z2.imag)) + _i*((z2.real * z1.imag - z1.real * z2.imag) / (z2.real * z2.real + z2.imag * z2.imag)));
}

complex1 complex1::Deg(complex1 z1, complex1 z2) {

	
	return complex1 (Exp(z2*Ln(z1)));
}

complex1 complex1::Deg1(complex1 z1, complex1 z2) {
	double fi;
	fi = 2 * atan((z1.imag) / (z1.real + abs(z1)));

	return complex1 (abs(z1) * abs(z1)*(cos(2*fi) + _i * sin(2*fi)));
}


complex1 complex1::Ln(complex1 z) {

return complex1(log(z.abs(z)), z.arg());
}
complex1 complex1::Exp(complex1 z) {
	
	return complex1 (exp(z.real) * cos(z.imag), exp(z.real) * sin(z.imag));
}
complex1 complex1::Sin(complex1 z) {

	complex1 c(0,-1),d, e, f, g;
	d = Exp(z * _i); e = Exp(z * c);
	
	return complex1 ((e - d)/ 2*_i);
}

complex1 complex1::Cos(complex1 z) {

	complex1 c(0, -1), d, e, f, g;
	d = Exp(z * _i); e = Exp(z * c);
	
	return complex1((d + e) / 2);
}

complex1 complex1::Root(complex1 z) {
	double fi;
	complex1 a1, a2;
	fi = 2 * atan((z.imag) / (z.real + abs(z)));
	a1 = (sqrt(abs(z)) * (cos((fi) / 2) + _i * sin((fi) / 2)));
	return a1;
	
}
/*
class privet {
public:
	void sayPrivet()
	{
		cout << "Privet";
	}
}; 
*/


inline complex1 operator *(const complex1& z1, double q) {
	return complex1(z1.real * q, z1.imag * q);
}

inline complex1 operator *(double q, const complex1& z1) {
	return complex1(q * z1.real, q * z1.imag);
}

inline complex1 operator *(const complex1& z1, const complex1& z2) {
	return complex1(z1.real * z2.real - z1.imag * z2.imag, z1.real * z2.imag + z2.real * z1.imag);
}

inline complex1 operator +(const complex1& z1, double q) {
	return complex1(z1.real + q, z1.imag);
}

inline complex1 operator +(double q, const complex1& z1) {
	return complex1(q + z1.real, z1.imag);
}

inline complex1 operator +(const complex1& z1, const complex1& z2) {
	return complex1(z1.real + z2.real, z1.imag + z2.imag);
}

inline complex1 operator /(const complex1& z1, double q) {
	
	return complex1((z1.real * q) / (q * q), (q * z1.imag) / (q * q)); 
}

inline complex1 operator /(double q, const complex1& z1) {

	return complex1((z1.real * q) / (q * q), (q * z1.imag) / (q * q));
}

inline complex1 operator /(const complex1& z1, const complex1& z2) {
	return complex1((z1.real * z2.real + z1.imag * z2.imag) / (z2.real * z2.real + z2.imag*z2.imag), (z2.real * z1.imag - z1.real * z2.imag) / (z2.real * z2.real + z2.imag * z2.imag));
}

inline complex1 operator -(const complex1& z1, double q) {
	return complex1(z1.real - q, z1.imag);
}

inline complex1 operator -(double q, const complex1& z1) {
	return complex1(q - z1.real, z1.imag);
}

inline complex1 operator -(const complex1& z1, const complex1& z2) {
	return complex1(z1.real - z2.real, z1.imag - z2.imag);
}
