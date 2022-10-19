#pragma once
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>


double f(double x);
double f1(double x, double y);

class IntegraL
{
public:
	double Rng1();
	double Trapize(double h, double x1, double x2);
	double Right(double h, double x1, double x2);
	double Mid(double h, double x1, double x2);
	double Left(double h, double x1, double x2);
	double Simp(double h, double x1, double x2);
	double Katess2(double h, double x1, double x2);
};

double IntegraL::Trapize(double h, double x1, double x2) {
	double S = 0.0;
	double i = 0.0;

	for (i = x1; i < x2; i = i + h)
	{
		S += (f(i) + f(i - h));
	}
	return  (S * (h / 2.0));
}

double IntegraL::Right(double h, double x1, double x2) {
	double S = 0.0;
	for (double i = x1; i < x2; i = i + h)
	{
		S += (f(i + h));
	}
	return  (S * h);
}

double IntegraL::Mid(double h, double x1, double x2) {
	double S = 0.0;
	for (double i = x1; i < x2; i = i + h)
	{
		S += (f(i + (h / 2.0)));
	}
	return (S * h);
}

double IntegraL::Rng1() {
	double S = 0.0;
	double x1, x2;
	int k = 2;
	x1 = (double)(rand() % 100) / 25 - 2;
	x2 = (double)(rand() % 100) / 25 - 2;


	if (((-pow(x1, 3) + pow(x2, 3)) < 1) && (-x1 + x2 < 1) && (-2 < x1 < 2) && (-2 < x2 < 2)) {

		++k;
	};


	for (double i = 0.0; i < 10000; i++)
	{
		S += f1(x1, x2);
		S = S / (k * k);
	}
	return (S);
};

double IntegraL::Left(double h, double x1, double x2) {
	double S = 0.0;
	
	for (double i = x1; i < x2; i = i + h)
	{
		S += (f(i));
	}
	return (S * h);
}

double IntegraL::Simp(double h, double x1, double x2) {
	double I=0.0;
	double odd = 0.0;
	double even = 0.0;
	for (double i = x1 + h; i < x2; i = i + h)
	{
		odd += f(i - h);
		even += f(i);
	}
	I = x1 + x2 + 4*odd + 2*even;
	return I * h / 6;
}

double IntegraL::Katess2(double h, double x1, double x2) {
	double S = 0.0;
	for (double i = x1; i < x2; i = i + h)
	{
		S += 4 * f(i) + 2 * f(i);
	}
	return (S * h / 6);
}