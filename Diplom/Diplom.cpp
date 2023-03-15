// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// REWENIE SKALYARNOY ZADA4i DIFRAKCII NA SYSTEME AKUSTI4ESKY MYAGKIH EKRANOV METODOM GALERKINA


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
const double lymda = 0.05;
const double a = 0;
const double b = 1;
const double c = 0;
const double d = 1;
complex1 A[n][n + 1];


//complex1 Ker(double x, double y) {
//    return(t.Exp(_i * k * (x - y)));
//}
//complex1 U0(double x) {
//
//    return 1 - lymda * t.Exp(_i * x) * _i * (t.Exp(-1 * _i) - 1);
//}
//complex1 Ker(double x, double y) {
//    return(_i *(x - y));
//}
//complex1 U0(double x) {
//
//    return 1 - lymda * _i *(x - 0.5);
//}
complex1 Ker(double x1, double y1, double x2, double y2) {
    return(_i * (x1 - y2));
}
complex1 U0(double x1, double x2) {

    return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
}
complex1 Ux(double x1, double x2) {

    complex1 ux(x1 * x2, 0);

    return ux;
}
complex1 Uy(double y1, double y2) {

    complex1 uy(y1 * y2, 0);

    return uy;
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

complex1 del(/*double xi,*/ int I, int J) { //vmesto etogo integral ot * itoi and jtoi basisnoy function
   /* complex1 S(0.0, 0.0);
    double h = (b - a) / n;
    for (double i = 0.0; i < 1.0; i += h)
    {
        S = S + phi(xi, I) * phi(xi, J);
    }
    return S;*/
    return complex1(I == J, 0.0);
}

//complex1 Mid(double xj, double xj1, double xi, double xi1 ) {   //gospodi, Ya je zabil daje chto ya doljen pisatb !!!
//                                    // nado bilo (b - a) / n
//    complex1 S(0.0, 0.0);
//    double h = (b - a) / n;
//    for (double i = a; i < b; i = i + h)
//    {
//        double a = i + (h / 2.0);
//        S = S + (Ker(a, 1.0));
//    }
//    return S * h;
//}
//


complex1 middlepryam2(double a, double b, double a1, double b1) { //noviy integral sdelat vmesto etogo
    double nn = 100, h, h1, x, x1; complex1 in(0.0, 0.0);
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
complex1 middlepryam1(double a, double b, double a1, double b1) {
    double nn = 10, h, h1, x, x1; complex1 in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x < b) {
        in = in + U0(x,x1);
        x = x + h;
    }
    return in * h;
}

int main() {
    //cout << " AAAAAAAAAAAAAA";

    double h1, h2, x[n + 1][n+1], xi[n][n]; complex1 c1[n];
    int i, j;

    h1 = (b - a) / n;
    h2 = (d - c) / n;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {
            x[i][j] = a + j * h1;
            x[j][i] = a + j * h1;

            //cout<< x[i]<<endl; 
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n+1; j++) {
            xi[i][j] = x[i][j] + (h1 / 2.0);
             cout<< xi[i][j] << " ";
        }
        cout << endl;
    }

   // printcomplex1(middlepryam2(xi[9][9], xi[9][10], xi[9][9], xi[9][10]));
    cout << "\n";
    for (i = 0; i < n; i++) {  //peredelat na galerkina// vrode teper on
        for (j = 0; j < n; j++) {
            //cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
            
           //A[i][j] = вот тут другое - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);

            A[i][j] = Ux(xi[i][j], xi[i][j+1]) * h1 - lymda * middlepryam2(xi[i][j], xi[i][j + 1], xi[i][j], xi[i][j + 1]);
        }
        //double Temp = (xi[i] * xi[i]) - lymda * ((xi[i] / 3) - 0.25);
       // complex Temp1(Temp, 0.0);
        //A[i][n] = middlepryam1(xi[i][j], xi[i][j + 1], xi[i][j], xi[i][j + 1]);
        A[i][n] = xi[i][j] * xi[i][j + 1] - _i * lymda * (3 * xi[i][j] - 2) / 12;
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
    cout << " ----------------------------------------------------- " << endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

           /* cout << real(A[i][j]) << " " << imag(A[i][j]) << "i";
            cout << " ";*/
            printcomplex1(A[i][j]);

        }
        c1[i] = A[i][n];
        cout << endl;
    }
    //cout << " ----------------------------------------------------- " << endl;
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n + 1; j++) {
    //        printcomplex1(un(xi[i][j], c1));
    //        cout << "  " /*<< "1" << endl*/;
    //
    //    }
    //    cout << endl;
    //}   
    cout << " ----------------------------------------------------- " << endl;
    for (int i = 0; i < n; i++)
    {
       printcomplex1(c1[i]); cout << "  " << "1" << endl;

    }
    return 0;

}

