// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// REWENIE SKALYARNOY ZADA4i DIFRAKCII NA SYSTEME AKUSTI4ESKY MYAGKIH EKRANOV METODOM GALERKINA

#include <fstream>
//#include "Complex.h"
#include "Complex (2).h"
//#include "IntegraL.h"
//#include <ccomplex>
//#include "Complex_v2.h"
//#include "Complex_v1.h"
#include <iostream>
#include <cmath>

using namespace std;

complex t;
double k = 1, pi = 4.0 * atan(1.0);
const int n = 10;
const double lymda = 0.5;
const double a1 = 0.0, a2 = 0.0;
const double b1 = 1.0, b2 = 1.0;
const double c1 = 0.0, c2 = 0.0;
const double d1 = 1.0, d2 = 1.0;
const double H11 = (b1 - a1) / n, H12 = (d1 - c1) / n;
const double H21 = (b2 - a2) / n, H22 = (d2 - c2) / n;
double x11[n + 1], x12[n + 1];
double x21[n + 1], x22[n + 1];
complex A1[n * n][n * n + 1], C1[n * n];
complex A2[n * n][n * n + 1], C2[n * n];

void printcomplex(complex z) {
    printf("(%5.3f, %5.3f) ", real(z), imag(z));
    //cout <<"(" << z.real << "+" << z.imag<< ""<< ") ";
}

//ядро
complex Ker(double x1, double y1, double x2, double y2) {
    //return(_i * (x1 - y2));
    double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    return exp(_i * k * rho) / (4.0 * pi * rho);
}
// правая часть
complex U0(double x1, double x2) {

    //return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
    //return 1 - lymda * _i * (x1 - 0.5);
    return exp(_i * k * x1);
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
    return((xi >= x11[i]) && (xi < x11[i + 1]));
}

//сюда приходит i и j как дабл, но здесь определены как инт, не будет ли проблемы? выглядит очень неправильно
double phi2(double xi1, double xi2, int i, int j) {
    return((xi1 >= x11[i]) && (xi1 < x11[i + 1]) && (xi2 >= x12[j]) && (xi2 < x12[j + 1]));
}

// это интеграл от правой части, но надо ли её интегрировать вообще...
complex middlepryam1(int i1, int j1, int i2, int j2) {
    double		aa1 = a1 + i1 * H11, bb1 = a1 + H11, cc1 = c1 + j1 * H12, dd1 = c1 + H12,
                aa2 = a2 + i2 * H21, bb2 = a2 + H21, cc2 = c2 + j2 * H22, dd2 = c2 + H22;
    double nn = 20.0, h1, h2, t1, t2;
    complex in(0.0, 0.0);
    h1 = (bb1 - aa1) / nn;
    h2 = (bb2 - aa2) / nn;

    for (int i1 = 0; i1 < nn; i1++) {
        for (int i2 = 0; i2 < nn; i2++) {
            t1 = aa1 + (i1 + 0.5) * h1;
            t2 = aa2 + (i2 + 0.5) * h2;
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
    double		aa1 = a1 + i1 * H11, bb1 = a1 + H11, cc1 = c1 + j1 * H12, dd1 = c1 + H12,
                aa2 = a2 + i2 * H21, bb2 = a2 + H21, cc2 = c2 + j2 * H22, dd2 = c2 + H22;

    int nn = 8;
    double h11, h12, h21, h22, t11, t12, t21, t22, rho;
    complex in(0.0, 0.0);
    h11 = (bb1 - aa1) / nn;
    h12 = (dd1 - cc1) / nn;
    h21 = (bb2 - aa2) / nn;
    h22 = (dd2 - cc2) / nn;

    for (int kk = 0; kk < nn; kk++)
    {
        for (int ll = 0; ll < nn; ll++)
        {
            for (int ii = 0; ii < nn; ii++)
            {
                for (int jj = 0; jj < nn; jj++) {

                    //!!!! phi2(double xi1, double xi2, int i, int j) * phi2(double xi1, double xi2, int i, int j)
                    t11 = aa1 + (jj + 0.5) * h11;
                    t12 = cc1 + (ll + 0.5) * h12;
                    t21 = aa2 + (ii + 0.5) * h21;
                    t22 = cc2 + (kk + 0.5) * h22;
                    rho = sqrt((t11 - t21) * (t11 - t21) + (t12 - t22) * (t12 - t22));
                    if (rho > 1e-7) in = in + Ker(t11, t12, t21, t22)/* * phi2(t1, t2, ii, jj) * phi2(t1, t2, ii, jj)*/;
                }
            }

        }

    }

    return in * h11 * h12 * h21 * h22;
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


complex un(double x1, double x2) {

    complex s(0.0, 0.0);
    for (int i = 0; i < n * n; i++) {
        s = s + C1[i] * phi2(x1, x2, i / n, i % n); // or vice versa ... i%n, i/n)

    }
    return(s);
}


void print_un(int pn) {
   // std::ofstream File1("../Matrix_1.txt");
    double t1, t2;
    printf("\n");
    for (int i1 = 0; i1 < pn; i1++) {
        for (int i2 = 0; i2 < pn; i2++) {
            t1 = a1 + (b1 - a1) / pn * i1;
            t2 = c1 + (d1 - c1) / pn * i2;
            printf("%6.3f ", abs(un(t1, t2)));
           /// File1 << abs(un(t1, t2)) << " ";
        }
        printf("\n");
       // File1 << std::endl;
    }
   // File1.close();
}

void Zapis_v_File(int pn) {
    std::ofstream File1("./Matrix_1.txt");
    FILE* tab_file;
    fopen_s(&tab_file, "result1.xls", "w");
    double t1, t2;
   // printf("\n");
    for (int i1 = 0; i1 < pn; i1++) {
        for (int i2 = 0; i2 < pn; i2++) {
            t1 = a1 + (b1 - a1) / pn * i1;
            t2 = c1 + (d1 - c1) / pn * i2;
           // printf("%6.3f ", abs(un(t1, t2)));
            File1 << abs(un(t1, t2)) << "\t";
            fprintf(tab_file, "%5.5f\t", abs(un(t1, t2)));
        }
        //printf("\n");
        File1 << "\n";
        fprintf(tab_file, "\n");
    }
    File1.close();
    fclose(tab_file);

    //for (int i1 = 0; i1 < pn; i1++) {
    //    for (int i2 = 0; i2 < pn; i2++) {
    //        fprintf(tab_file, "%5.5f\t", abs(un(t1, t2))); 
    //    }
    //    fprintf(tab_file, "\n");
    //}
    //fclose(tab_file);
}

int main() {
    //cout << " AAAAAAAAAAAAAA";

    //double xi2[n + 1], xi1[n + 1];
   
   

    for (int j = 0; j < n + 1; j++) {
        x11[j] = a1 + j * H11;
        x12[j] = c1 + j * H12;
        cout << x11[j] << "   " << x12[j] << endl;
    }

    std::cout << " ----------------------------------------------------- " << std::endl;
    int i, j;
    // для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
    for (int i1 = 0; i1 < n; i1++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int i2 = 0; i2 < n; i2++) {
                for (int j2 = 0; j2 < n; j2++) {

                    i = i1 + n * j1;
                    j = i2 + n * j2;
                    //A[i][j] = del2(i1, j1, i2, j2) * h1 * h2 - lymda * middlepryam2(i1, j1, i2, j2);

                    A1[i][j] = middlepryam2(i1, j1, i2, j2);
                    //printcomplex(A1[i][j]);
                    A1[i][n * n] = middlepryam1(i1, j1, i2, j2); // перенес сюда потому что не было в j2 вне цикла
                    /*phi2(x1[i1], x2[j1], i1, j1)* phi2(x1[i2], x2[j2], i2, j2)*/    /*xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], xi2[j], xi2[j + 1], xi2[i], xi2[i + 1],*/
                }

            }
        }
    }

    //for (int i = 0; i < n * n; i++) {
    //    for (int j = 0; j < n * n + 1; j++) {
    //        printcomplex(A[i][j]);
    //    }
    //    cout << endl;
    //}
    cout << endl;
    Gauss(0, A1);
    cout << " ----------------------------------------------------- " << endl;

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
        C1[i] = A1[i][n * n];
        printcomplex(C1[i]); cout << "  " << "1" << endl;

    }

    print_un(15);
    Zapis_v_File(15);

    //system("pause");
    return 0;

}

