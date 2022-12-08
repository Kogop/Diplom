﻿// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
//#include <complex.h>
#include "Complex.h"
//#include "IntegraL.h"
//#include <ccomplex>
//#include "Complex_v2.h"
//#include "Complex_v1.h"
#include <iostream>
#include <cmath>
using namespace std;

complex1 t;
double k = 1;
const int n = 10;
const int lymda = 0.1;
const double a = 0;
const double b = 1;
complex1 A[n][n + 1];


complex1 Ker(double x, double y) {
    return(t.Exp(_i * k * abs(x - y)));
}
complex1 U0(double x) {

    return 1 - lymda * t.Exp(_i * x) * _i * (t.Exp(_i) - 1);
}

//pravilno?
double phi(double xi, int i) { //poka odnomernoe potom peredelat nado na 2 mernoe
    double x[n + 1], h, s;
    int j;

    h = (b - a) / n;

    for (j = 0; j < n + 1; j++) {

        x[j] = a + j * h;

    }
    /*if ((xi >= x[i]) && (xi <= x[i + 1])) {
        s = 1;
    }
    else {
        s = 0;
    }*/
    return((xi >= x[i]) && (xi <= x[i + 1]));
}

complex1 del(double xi, int I, int J) { //vmesto etogo integral ot * itoi and jtoi basisnoy function
    complex1 S(0.0, 0.0);
    double h = (b - a) / n;
    for (double i = 0.0; i < 1.0; i += h)
    {
        S = S + phi(xi, I) * phi(xi, J);
    }
    return S;
}

complex1 Mid(double a, double b) {   //gospodi, Ya je zabil daje chto ya doljen pisatb !!!
                                    // nado bilo (b - a) / n
    complex1 S(0.0, 0.0);
    double h = (b - a) / n;
    for (double i = a; i < b; i = i + h)
    {
        double a = i + (h / 2.0);
        S = S + (Ker(a, 1.0));
    }
    return S * h;
}



complex1 middlepryam2(double a, double b, double a1, double b1) { //noviy integral sdelat vmesto etogo
    double nn = 10, h, h1, x, x1; complex1 in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 <= b1 - (h1 / 2))
    {
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
complex1 un(double xi, complex1 c[n]) {
    int i;
    complex1 s(0.0, 0.0);
    for (i = 0; i < n; i++) {
        s = s + c[i] * phi(xi, i);

    }
    return(s);
}



int main() {
    //cout << " AAAAAAAAAAAAAA";

    double h, x[n + 1], xi[n]; complex1 c[n];
    int i, j;

    h = (b - a) / n;

    for (i = 0; i < n + 1; i++) {

        x[i] = a + i * h;

         //cout<< x[i]<<endl; 
    }
    for (i = 0; i < n; i++) {

        xi[i] = x[i] + (h / 2);
        // cout<< xi[i]<<endl; 
    }
    for (i = 0; i < n; i++) {  //peredelat na galerkina
        for (j = 0; j < n; j++) {
            //cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
            //A[i][j] = вот тут другое - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);


            //A[i][j] = del(xi[i], i, j) * h - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);


            // 
            //cout << real(del(i, j) * h) << " <- JOPA "<< imag(del(i, j) * h) << "<-JOPA pomenwe " << real(lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1])) << "<__ POPA " << imag(lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1])) << " POPOPOPO " << endl;
            A[i][j] = del(xi[i], i, j) * h - lymda * Mid(a, b);
        }
        //double Temp = (xi[i] * xi[i]) - lymda * ((xi[i] / 3) - 0.25);
       // complex Temp1(Temp, 0.0);
        A[i][n] = U0(xi[i]);
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

          /*  cout << real(A[i][j]) << " " << imag(A[i][j])<< "i";
            cout << " ";*/
            printcomplex1(A[i][j]);
        }
        cout << endl;
    }
    cout << endl;
    Gauss(0, A);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

           /* cout << real(A[i][j]) << " " << imag(A[i][j]) << "i";
            cout << " ";*/
            printcomplex1(A[i][j]);

        }
        c[i] = A[i][n];
        cout << endl;
    }

    for (i = 0; i < n; i++) {
        printcomplex1(un(xi[i], c));
        cout << "  " << pow(xi[i], 2) << endl;

    }
    return 0;

}

