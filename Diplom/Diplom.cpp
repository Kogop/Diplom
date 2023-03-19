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

//ядро
complex1 Ker(double x1, double y1, double x2, double y2) {
    return(_i * (x1 - y2));
}
// правая часть
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
double phi2(double xi1, double xi2, int i,int j) { //ETO 2 mernoe??? ILI NADO 2 x proveryat
    double x[n + 1], h, h1, s;
    
    h = (b - a) / n;
    h1 = (d - c) / n;
    for (int k = 0; k < n + 1; k++) {

        x[k] = a + k * h;

    }
    /*if ((xi >= x[i]) && (xi <= x[i + 1])) {
        s = 1;
    }
    else {
        s = 0;
    }*/
    return((xi1 >= x[i]) && (xi1 <= x[i + 1]) && (xi2 >= x[j]) && (xi2 <= x[j + 1]));
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
complex1 middlepryam2(double a, double b, double a1, double b1, int i, int j) { 
    double nn = 100, h, h1, x, x1; complex1 in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    //for (double kk = a1; kk < b1; kk = kk + h1)
    //{
    //    for (double ll = a; ll < b; ll = ll + h) 
    //    {
            for (double ii = a1; ii < b1; ii = ii + h1)
            {
                for (double jj = a; jj < b; jj = jj + h) {
                    in = in + (Ker(x, 0, 0, x1)) * phi(x, i) * phi(x1, j);
                    x = x + h;
                }
                x1 = x1 + h1;
            }
    //    }
    //}
    //while (x1 < b1)
    //{
    //    while (x < b) {
    //        in = in + (Ker(x1, x, x1, x))*phi(x,i)*phi(x1,j);
    //        x = x + h;
    //    }
    //    x1 = x1 + h1;
    //}
    return in * h * h1;
}


//это первая половина матрицы, где двойной интеграл от фи
// могу ли я тут тоже на дельту это все заменить?
complex1 middlepryam2_phi(double a, double b, double a1, double b1,int i,int j) { 
    double nn = 100, h, h1, x, x1; complex1 in(0.0, 0.0);
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 < b1)
    {
        while (x < b) {
            in = in + (phi(x,i)*phi(x1,j));
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
    
    complex1 s(0.0, 0.0);
    for (int i = 0; i < n; i++) {
        s = s + c[i] * phi(xi, i);

    }
    return(s);
}

// это интеграл от правой части, но надо ли её интегрировать вообще...
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

    double h1, h2, x[n + 1], xi[n+1]; complex1 c1[n];
    int i, j;

    h1 = (b - a) / n;
    h2 = (d - c) / n;

   
    for (j = 0; j < n + 1; j++) {
        x[j] = a + j * h1;

        //cout<< x[i]<<endl; 
    }

    for (j = 0; j < n + 1; j++) {
        xi[j] = x[j] + (h1 / 2.0);
        std::cout << xi[j] << " ";

       // 
    }
    std::cout << endl;
    std::cout << " ----------------------------------------------------- " << std::endl;
    //printcomplex1(middlepryam2(xi[9][9], xi[9][10], xi[9][9], xi[9][10]));
    //cout << "\n";
 
    // для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
    for (i = 0; i < n; i++) {  //peredelat na galerkina// vrode teper on
        for (j = 0; j < n; j++) {
            //cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
            
           //A[i][j] = вот тут другое - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);

            //A[i][j] = Ux(xi[j], xi[j+1]) * h1*h1 - lymda * middlepryam2(xi[j], xi[j + 1], xi[j], xi[j + 1]);
            A[i][j] = middlepryam2_phi(xi[j], xi[j + 1], xi[i], xi[i + 1], i, j) * h1 * h1 - lymda * middlepryam2(xi[j], xi[j + 1], xi[i], xi[i + 1], i, j);
        }
        //double Temp = (xi[i] * xi[i]) - lymda * ((xi[i] / 3) - 0.25);
       // complex Temp1(Temp, 0.0);
        A[i][n] = middlepryam1(xi[i], xi[i + 1], xi[i], xi[i + 1]);
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

