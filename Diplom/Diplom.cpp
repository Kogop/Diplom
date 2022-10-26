// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
//#include <complex.h>
#include "Complex.h"
//#include "IntegraL.h"
#include <iostream>
#include <cmath>


complex1 t;
double k = 1;
const int n = 50;
const int lymda = 2;
const double a = 0;
const double b = 1;
complex1 A[n][n + 1];


complex1 Ker(double x, double y) {
    return(t.Exp(_i * k * abs(x * y)));
}


//nado peredelat skoree vsego, 4et tut ne to voobshe. na dvumernuu matrix nado
double phi(double xi, int i) {
    double x[n + 1], h, s;
    int j;

    h = (b - a) / n;

    for (j = 0; j < n + 1; j++) {

        x[j] = a + j * h;

    }
    if ((xi >= x[i]) && (xi <= x[i + 1])) {
        s = 1;
    }
    else {
        s = 0;
    }
    return(s);
}

complex1 middlepryam2(double a, double b, double a1, double b1) {
    double nn = 1000, h, h1, x, x1; complex1 in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 <= b1 - (h1 / 2)) {
        while (x <= b - (h / 2)) {
            in = in + (Ker(x1, x) * h);
            x = x + h;


        }
        x1 = x1 + h1;
    }
    return in;
}

void Gauss(int k, complex1 Matrix[n][n + 1]) {
    if (Matrix[k][k] != complex1(1.0, 0.0)) {
        complex1 T = Matrix[k][k];
        for (int j = k; j < n + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n; i++) {
        if ((Matrix[i][k] != complex1(0.0, 0.0)) && (i != k)) {
            complex1 T = Matrix[i][k];
            Matrix[i][k] = complex1(0.0, 0.0);
            for (int j = k + 1; j < n + 1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n - 1) {
        Gauss(k + 1, Matrix);
    }
}



int main()
{
     printcomplex1(Ker(1, 2));
}

