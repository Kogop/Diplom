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
const int n = 14;
const double lymda = 0.125;
const double a = 0.0;
const double b = 1.0;
const double c = 0.0;
const double d = 1.0;
const double h1 = (b - a) / n, h2 = (d - c) / n;
double x1[n + 1], x2[n + 1];
complex1 A[n*n][n*n + 1];


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

//ядро
complex1 Ker(double x1, double y1, double x2, double y2) {
    return(_i * (x1 - y2));
}
// правая часть
complex1 U0(double x1, double x2) {

    //return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
    return 1 - lymda * _i * (x1 - 0.5);
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

    /*if ((xi >= x[i]) && (xi <= x[i + 1])) {
        s = 1;
    }
    else {
        s = 0;
    }*/
    return((xi >= x1[i]) && (xi < x1[i + 1]));
}

//сюда приходит i и j как дабл, но здесь определены как инт, не будет ли проблемы? выглядит очень неправильно
double phi2(double xi1, double xi2, int i, int j) { 
    return((xi1 >= x1[i]) && (xi1 < x1[i + 1]) && (xi2 >= x2[j]) && (xi2 < x2[j + 1]));
}

// это интеграл от правой части, но надо ли её интегрировать вообще...
complex1 middlepryam1(double a1, double b1, double a2, double b2) {
    double nn = 20.0, h1, h2, t1, t2;
    complex1 in(0.0, 0.0);
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

//правильно?
complex1 del2(int I1, int J1, int I2, int J2) {

    return complex1(I1 == I2 && J1 == J2, 0.0);
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

// запасная копия если испорчу рабочее
complex1 middlepryam2_save_copy(double a, double b, double a1, double b1) { //noviy integral sdelat vmesto etogo
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



//это интеграл с ядром, почему-то не хочет работать с 4 интегралами
// возможно потому что нету элементов? тк х это х1(из ядра) а х1 в данном случает это у2 (из ядра)?
//как правильно для каждого интеграла брать границы. Передавать их все очень муторно, передать начало конец считать внутри фунции?
//возможно ли передать фунцию в фунцию, чтобы не приходилось переделывать интеграл каждый раз как меняю ядро?

// ДОДЕЛАТЬ!!!
//теперь работает, чутка медленно, но работает, правильно? не знаю.
complex1 middlepryam2(/*double a1, double b1, double c1, double d1,
    double a2, double b2, double c2, double d2,*/
    int i1, int j1,
    int i2, int j2)

    //    complex1 middlepryam2(int i1, int j1, int i2, int j2) {
    //    a1 = a + i1 * h1, b1 = a1 + h1,  c1 = c + j1 * h2,  d1 = c1 + h2,
    //        a2 = a + i2 * h1, b2 = a2 + h1, c2 = c + j2 * h2, d2 = c2 + h2,
    //}

{
    double a1 = a + i1 * h1, b1 = a1 + h1, c1 = c + j1 * h2, d1 = c1 + h2,
        a2 = a + i2 * h1, b2 = a2 + h1, c2 = c + j2 * h2, d2 = c2 + h2;

    double nn = 10.0, hh1, hh2, hh3, hh4, t1, t2, t3, t4;
    complex1 in(0.0, 0.0);
    hh1 = (b1 - a1) / nn;
    hh2 = (b2 - a2) / nn;
    hh3 = (d1 - c1) / nn;
    hh4 = (d2 - c2) / nn;
    //t1 = a1 + (h1 / 2.0);
    //t2 = c1 + (h2 / 2.0);
    //t3 = a2 + (h3 / 2.0);
    //t4 = c2 + (h4 / 2.0);
    //cout<<" x= "<<x<<endl;
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
            /*t1 = a1 + (h1 / 2);
            t2 = c1 + (h2 / 2);*/
        }
        /*t1 = a1 + (h1 / 2);
        t2 = c1 + (h2 / 2);*/
    }

    return in * hh1 * hh2 * hh3 * hh4;
}


//это первая половина матрицы, где двойной интеграл от фи
// могу ли я тут тоже на дельту это все заменить?
complex1 middlepryam2_phi(double a, double b, double a1, double b1, int i, int j) {
    double nn = 100, h, h1, x, x1; complex1 in(0.0, 0.0);
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

void Gauss(int k, complex1 Matrix[n * n][n * n + 1]) {
    if (Matrix[k][k] != complex1(1.0, 0.0)) {
        complex1 T = Matrix[k][k];
        for (int j = k; j < n * n + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n * n; i++) {
        if ((Matrix[i][k] != complex1(0.0, 0.0)) && (i != k)) {
            complex1 T = Matrix[i][k];
            Matrix[i][k] = complex1(0.0, 0.0);
            for (int j = k + 1; j < n * n + 1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n * n - 1) {
        Gauss(k + 1, Matrix);
    }
}


complex1 un(double xi1, double xi2, complex1 c[n*n]) {

    complex1 s(0.0, 0.0);
    for (int i = 0; i < n*n; i++) {
        s = s + c[i] * phi2(xi1, xi2, i/n, i%n); // or vice versa ... i%n, i/n)

    }
    return(s);
}


int main() {
    //cout << " AAAAAAAAAAAAAA";

    double xi2[n + 1], xi1[n + 1]; complex1 c1[n * n];
    int i, j;

    for (int j = 0; j < n + 1; j++) {
        x1[j] = a + j * h1;
        x2[j] = c + j * h2;
        cout << x1[j] << "   " << x2[j] << endl;
    }

    //system("pause");
    //printcomplex1(middlepryam1(x1[0], x1[1], x2[0], x2[1])); cout << endl;
    //system("pause");


    cout << " ----------------------------------------------------- " << endl;
    //for (int j = 0; j < n + 1; j++) {
    //    xi1[j] = x1[j] + (h1 / 2.0);
    //    xi2[j] = x2[j] + (h2 / 2.0);
    //    std::cout << xi1[j] << " " << xi2[j] << " ";
    //    std::cout << endl;
    //    // 
    //}

    std::cout << endl;
    std::cout << " ----------------------------------------------------- " << std::endl;
    //printcomplex1(middlepryam2(xi[9][9], xi[9][10], xi[9][9], xi[9][10]));
    //cout << "\n";

    // для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
    //добавить 2 цикла 
    for (int i1 = 0; i1 < n; i1++) {  //peredelat na galerkina// vrode teper on
        for (int j1 = 0; j1 < n; j1++) {
            for (int i2 = 0; i2 < n; i2++) {  //peredelat na galerkina// vrode teper on
                for (int j2 = 0; j2 < n; j2++) {
                    //cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
                   //A[i][j] = вот тут другое - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
                    //A[i][j] = Ux(xi[j], xi[j+1]) * h1*h1 - lymda * middlepryam2(xi[j], xi[j + 1], xi[j], xi[j + 1]);

                    i = i1 + n * j1;
                    j = i2 + n * j2;
                    // закоментированное это то что было, снизу на дельту переделал. Все еще не понимаю почему i и j одни и теже используются.
                    //A[i][j] = middlepryam2_phi(xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], i, j) * h1 * h1 - lymda * middlepryam2(xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], xi2[j], xi2[j + 1], xi2[i], xi2[i + 1], i, j,i,j);
                    A[i][j] =  del2(i1, j1, i2, j2) * h1 * h2 - lymda * middlepryam2( i1, j1, i2, j2);
                    /*phi2(x1[i1], x2[j1], i1, j1)* phi2(x1[i2], x2[j2], i2, j2)*/    /*xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], xi2[j], xi2[j + 1], xi2[i], xi2[i + 1],*/
                }
                //double Temp = (xi[i] * xi[i]) - lymda * ((xi[i] / 3) - 0.25);
               // complex Temp1(Temp, 0.0);
                A[i][n * n] = middlepryam1(x1[i1], x1[i1 + 1], x2[j1], x2[j1 + 1]);

            }
        }
    }

    //for (int i = 0; i < n * n; i++) {
    //    for (int j = 0; j < n * n + 1; j++) {

    //        /*  cout << real(A[i][j]) << " " << imag(A[i][j])<< "i";
    //          cout << " ";*/
    //        printcomplex1(A[i][j]);
    //    }
    //    cout << endl;
    //}
    //cout << endl;
    Gauss(0, A);
    cout << " ----------------------------------------------------- " << endl;
    for (int i = 0; i < n * n; i++) {
        for (int j = 0; j < n * n + 1; j++) {

            /* cout << real(A[i][j]) << " " << imag(A[i][j]) << "i";
             cout << " ";*/
           // printcomplex1(A[i][j]);

        }
        c1[i] = A[i][n * n];
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
    for (int i = 0; i < n * n; i++)
    {
        printcomplex1(c1[i]); cout << "  " << "1" << endl;

    }
    //system("pause");
    return 0;

}

