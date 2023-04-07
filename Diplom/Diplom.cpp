// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// REWENIE SKALYARNOY ZADA4i DIFRAKCII NA SYSTEME AKUSTI4ESKY MYAGKIH EKRANOV METODOM GALERKINA


//#include <complex.h>
#include "Complex (2).h"
//#include "IntegraL.h"
//#include <ccomplex>
//#include "Complex_v2.h"
//#include "Complex_v1.h"
#include <iostream>
#include <cmath>
using namespace std;

complex t;
double k = 1;
const int n = 10;
const double lymda = 0.5;
const double a = 0.0;
const double b = 1.0;
const double c = 0.0;
const double d = 1.0;
const double h1 = (b - a) / n, h2 = (d - c) / n;
double x1[n + 1], x2[n + 1];
complex A[n * n][n * n + 1];

void printcomplex(complex z) {
    printf("(%5.10f, %5.10f)", real(z), imag(z));
    //cout <<"(" << z.real << "+" << z.imag<< ""<< ") ";
}

//ядро
complex Ker(double x1, double y1, double x2, double y2) {
    return(_i * (x1 - y2));
}
// правая часть
complex U0(double x1, double x2) {

    //return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
    return 1 - lymda * _i * (x1 - 0.5);
}
complex Ux(double x1, double x2) {

    complex ux(x1 * x2, 0);

    return ux;
}
complex Uy(double y1, double y2) {

    complex uy(y1 * y2, 0);

    return uy;
}


//pravilno?
double phi(double xi, int i) { //poka odnomernoe potom peredelat nado na 2 mernoe
    return((xi >= x1[i]) && (xi < x1[i + 1]));
}

//сюда приходит i и j как дабл, но здесь определены как инт, не будет ли проблемы? выглядит очень неправильно
double phi2(double xi1, double xi2, int i, int j) {
    return((xi1 >= x1[i]) && (xi1 < x1[i + 1]) && (xi2 >= x2[j]) && (xi2 < x2[j + 1]));
}

// это интеграл от правой части, но надо ли её интегрировать вообще...
complex middlepryam1(double a1, double b1, double a2, double b2) {
    double nn = 20.0, h1, h2, t1, t2;
    complex in(0.0, 0.0);
    h1 = (b1 - a1) / nn;
    h2 = (b2 - a2) / nn;

    for (int i1 = 0; i1 < nn; i1++) {
        for (int i2 = 0; i2 < nn; i2++) {
            t1 = a1 + (i1 + 0.5) * h1;
            t2 = a2 + (i2 + 0.5) * h2;
            in = in + U0(t1, t2);
        }
    }
    // in = in * h1 * h2;
    return in * h1 * h2;
}


complex del(/*double xi,*/ int I, int J) { //vmesto etogo integral ot * itoi and jtoi basisnoy function
    return complex(I == J, 0.0);
}

//правильно?
complex del2(int I1, int J1, int I2, int J2) {

    return complex(I1 == I2 && J1 == J2, 0.0);
}


// запасная копия если испорчу рабочее
complex middlepryam2_save_copy(double a, double b, double a1, double b1) { //noviy integral sdelat vmesto etogo
    double nn = 10, h, h1, x, x1; complex in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 < b1)
    {
        while (x < b) {
            in = in + (Ker(x1, x, x1, x));
            x = x + h;
        }
        x1 = x1 + h1;
    }
    return in * h * h1;
}


complex middlepryam2(int i1, int j1, int i2, int j2) {
    double		a1 = a + i1 * h1, b1 = a1 + h1, c1 = c + j1 * h2, d1 = c1 + h2,
        a2 = a + i2 * h1, b2 = a2 + h1, c2 = c + j2 * h2, d2 = c2 + h2;

    int nn = 8;
    double hh1, hh2, hh3, hh4, t1, t2, t3, t4;
    complex in(0.0, 0.0);
    hh1 = (b1 - a1) / nn;
    hh2 = (b2 - a2) / nn;
    hh3 = (d1 - c1) / nn;
    hh4 = (d2 - c2) / nn;

    for (int kk = 0; kk < nn; kk++)
    {
        for (int ll = 0; ll < nn; ll++)
        {
            for (int ii = 0; ii < nn; ii++)
            {
                for (int jj = 0; jj < nn; jj++) {

                    //!!!! phi2(double xi1, double xi2, int i, int j) * phi2(double xi1, double xi2, int i, int j)
                    t1 = a1 + (jj + 0.5) * hh1;
                    t2 = a2 + (ii + 0.5) * hh2;
                    t3 = c1 + (ll + 0.5) * hh3;
                    t4 = c2 + (kk + 0.5) * hh4;
                    in = in + Ker(t1, t2, t3, t4)/* * phi2(t1, t2, ii, jj) * phi2(t1, t2, ii, jj)*/;
                }
            }

        }

    }

    return in * hh1 * hh2 * hh3 * hh4;
}


//это первая половина матрицы, где двойной интеграл от фи
// могу ли я тут тоже на дельту это все заменить?
complex middlepryam2_phi(double a, double b, double a1, double b1, int i, int j) {
    double nn = 100, h, h1, x, x1; complex in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 < b1)
    {
        while (x < b) {
            in = in + (phi(x, i) * phi(x1, j));
            //!!!! phi2(double xi1, double xi2, int i, int j) * phi2(double xi1, double xi2, int i, int j)
            x = x + h;
        }
        x1 = x1 + h1;
    }
    return in * h * h1;
}

void Gauss(int k, complex Matrix[n * n][n * n + 1]) {
    if (Matrix[k][k] != complex(1.0, 0.0)) {
        complex T = Matrix[k][k];
        for (int j = k; j < n * n + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n * n; i++) {
        if ((Matrix[i][k] != complex(0.0, 0.0)) && (i != k)) {
            complex T = Matrix[i][k];
            Matrix[i][k] = complex(0.0, 0.0);
            for (int j = k + 1; j < n * n + 1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n * n - 1) {
        Gauss(k + 1, Matrix);
    }
}


complex un(double xi1, double xi2, complex c[n * n]) {

    complex s(0.0, 0.0);
    for (int i = 0; i < n * n; i++) {
        s = s + c[i] * phi2(xi1, xi2, i / n, i % n); // or vice versa ... i%n, i/n)

    }
    return(s);
}


int main() {
    //cout << " AAAAAAAAAAAAAA";

    double xi2[n + 1], xi1[n + 1]; complex c1[n * n];
    int i, j;

    for (int j = 0; j < n + 1; j++) {
        x1[j] = a + j * h1;
        x2[j] = c + j * h2;
        cout << x1[j] << "   " << x2[j] << endl;
    }

    std::cout << " ----------------------------------------------------- " << std::endl;

    // для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
    for (int i1 = 0; i1 < n; i1++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int i2 = 0; i2 < n; i2++) {
                for (int j2 = 0; j2 < n; j2++) {

                    i = i1 + n * j1;
                    j = i2 + n * j2;
                    A[i][j] = del2(i1, j1, i2, j2) * h1 * h2 - lymda * middlepryam2(i1, j1, i2, j2);
                    /*phi2(x1[i1], x2[j1], i1, j1)* phi2(x1[i2], x2[j2], i2, j2)*/    /*xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], xi2[j], xi2[j + 1], xi2[i], xi2[i + 1],*/
                }
                A[i][n * n] = middlepryam1(x1[i1], x1[i1 + 1], x2[j1], x2[j1 + 1]);
            }
        }
    }

    for (int i = 0; i < n * n; i++) {
        for (int j = 0; j < n * n + 1; j++) {
            printcomplex(A[i][j]);
        }
        cout << endl;
    }
    cout << endl;
    Gauss(0, A);
    cout << " ----------------------------------------------------- " << endl;
    for (int i = 0; i < n * n; i++) {
        for (int j = 0; j < n * n + 1; j++) {

            /* cout << real(A[i][j]) << " " << imag(A[i][j]) << "i";
             cout << " ";*/
             // printcomplex(A[i][j]);

        }
        c1[i] = A[i][n * n];
       // cout << endl;
    }



    //cout << " ----------------------------------------------------- " << endl;
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n + 1; j++) {
    //        printcomplex(un(xi[i][j], c1));
    //        cout << "  " /*<< "1" << endl*/;
    //
    //    }
    //    cout << endl;
    //}   
    cout << " ----------------------------------------------------- " << endl;
    for (int i = 0; i < n * n; i++)
    {
        printcomplex(c1[i]); cout << "  " << "1" << endl;

    }
    //system("pause");
    return 0;

}

