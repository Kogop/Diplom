// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// REWENIE SKALYARNOY ZADA4i DIFRAKCII NA SYSTEME AKUSTI4ESKY MYAGKIH EKRANOV METODOM GALERKINA


//начать переводить в трехмерную задачу, параметры, все такое, начинать дипломную записку.
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
const int n = 7, N = 2 * n * n;
const double lymda = 0.5;
const double GranA1 = 0.0, GranA2 = 2.0;
const double GranB1 = 1.0, GranB2 = 4.0;
const double GranC1 = 0.0, GranC2 = 2.0;
const double GranD1 = 1.0, GranD2 = 4.0;
const double H11 = (GranB1 - GranA1) / n, H12 = (GranD1 - GranC1) / n;
const double H21 = (GranB2 - GranA2) / n, H22 = (GranD2 - GranC2) / n;
double X11[n + 1], X12[n + 1];
double X21[n + 1], X22[n + 1];
//complex A1[n * n][n * n + 1], C1[n * n];
//complex A2[n * n][n * n + 1], C2[n * n];

complex A[N][N + 1], C[N]; // !!!
double u1[n + 1], v1[n + 1];
double u2[n + 1], v2[n + 1];


void printcomplex(complex z) {
	printf("(%5.3f, %5.3f) ", real(z), imag(z));
	//cout <<"(" << z.real << "+" << z.imag<< ""<< ") ";
}

//ядро
complex Ker(double x1, double y1, double x2, double y2) {   //добавить z1, z2, тк теперь будет объемное тело
	//return(_i * (x1 - y2));
	double rho = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	return exp(_i * k * rho) / (4.0 * pi * rho);
}
// правая часть
complex U0(double x1, double x2) {
	//return (x1 * x2) - (_i * lymda * (3.0 * x1 - 2.0)) / 12.0;
	//return 1 - lymda * _i * (x1 - 0.5);
	return exp(_i * k * x1);
	//return _i;
}
complex Ux(double x1, double x2) {
	complex ux(x1 * x2, 0);

	return ux;
}
complex Uy(double y1, double y2) {
	complex uy(y1 * y2, 0);

	return uy;
}

complex X_Param(double u, double v) {
	return u;
}
complex Y_Param(double u, double v) {
	return v;
}
complex Z_Param(double u, double v) {
	return 0;
}
complex DX_Param(double u, double v) {
	return (X_Param(u + 0.0000001, v) - X_Param(u - 0.0000001, v)) / (2 * 0.0000001);
}
complex DY_Param(double u, double v) {
	return (Y_Param(u + 0.0000001, v) - Y_Param(u - 0.0000001, v)) / (2 * 0.0000001);
}
complex DZ_Param(double u, double v) {
	return (Z_Param(u + 0.0000001, v) - Z_Param(u - 0.0000001, v)) / (2 * 0.0000001);
}

complex sqrtEGF2(double u1, double v1, double u2, double v2) {

	complex E, G, F, root;
	E = DX_Param(u1, v1) + DY_Param(u1, v1) + DZ_Param(u1, v1);
	G = DX_Param(u2, v2) + DY_Param(u2, v2) + DZ_Param(u2, v2);
	F = DX_Param(u1, v2) + DY_Param(u1, v2) + DZ_Param(u1, v2);
	return Root(E * G - F * F);
}



//pravilno?
double phi(double xi, int i) { //poka odnomernoe potom peredelat nado na 2 mernoe
	return((xi >= X11[i]) && (xi < X11[i + 1]));
}

//сюда приходит i и j как дабл, но здесь определены как инт, не будет ли проблемы? выглядит очень неправильно
double phi2(double xi1, double xi2, int i, int j, int f) {
	if (f)
	{
		return((xi1 >= X11[i]) && (xi1 < X11[i + 1]) && (xi2 >= X12[j]) && (xi2 < X12[j + 1]));
	}
	else {
		return((xi1 >= X21[i]) && (xi1 < X21[i + 1]) && (xi2 >= X22[j]) && (xi2 < X22[j + 1]));
	}
}

// f = 1 если первая матрица, если нет то вторая
complex middlepryam1(int i1, int j1, int f) {
	double aa1, bb1, cc1, dd1;
	if (f)
	{
		aa1 = GranA1 + i1 * H11; bb1 = aa1 + H11; cc1 = GranC1 + j1 * H12; dd1 = cc1 + H12;
	}
	else {
		aa1 = GranA2 + i1 * H21; bb1 = aa1 + H21; cc1 = GranC2 + j1 * H22; dd1 = cc1 + H22;
	};

	double nn = 20.0, h1, h2, t1, t2;
	complex in(0.0, 0.0);
	h1 = (bb1 - aa1) / nn;
	h2 = (dd1 - cc1) / nn;

	for (int i1 = 0; i1 < nn; i1++) {
		for (int i2 = 0; i2 < nn; i2++) {
			t1 = aa1 + (i1 + 0.5) * h1;
			t2 = cc1 + (i2 + 0.5) * h2;
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


// РАБОЧАЯ КОПИЯ ИНТЕГРАЛА

//// f = 1 если первая матрица, если нет то вторая
//complex middlepryam2(int i1, int j1, int i2, int j2, int f1, int f2) {
//	double aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2;
//	if (f1)
//	{
//		aa1 = GranA1 + i1 * H11; bb1 = aa1 + H11; cc1 = GranC1 + j1 * H12; dd1 = cc1 + H12;
//	}
//	else {
//		aa1 = GranA2 + i1 * H21; bb1 = aa1 + H21; cc1 = GranC2 + j1 * H22; dd1 = cc1 + H22;
//	};
//	if (f2)
//	{
//		aa2 = GranA1 + i2 * H11; bb2 = aa1 + H11; cc2 = GranC1 + j2 * H12; dd2 = cc2 + H12;
//	}
//	else {
//		aa2 = GranA2 + i2 * H21; bb2 = aa1 + H21; cc2 = GranC2 + j2 * H22; dd2 = cc2 + H22;
//	};
//
//	int nn = 8;
//	double h11, h12, h21, h22, t11, t12, t21, t22, rho;
//	complex in(0.0, 0.0);
//	h11 = (bb1 - aa1) / nn;
//	h12 = (dd1 - cc1) / nn;
//	h21 = (bb2 - aa2) / nn;
//	h22 = (dd2 - cc2) / nn;
//
//	for (int kk = 0; kk < nn; kk++)
//	{
//		for (int ll = 0; ll < nn; ll++)
//		{
//			for (int ii = 0; ii < nn; ii++)
//			{
//				for (int jj = 0; jj < nn; jj++) {
//					//!!!! phi2(double xi1, double xi2, int i, int j) * phi2(double xi1, double xi2, int i, int j)
//					t11 = aa1 + (jj + 0.5) * h11;
//					t12 = cc1 + (ll + 0.5) * h12;
//
//					t21 = aa2 + (ii + 0.5) * h21;
//					t22 = cc2 + (kk + 0.5) * h22;
//					rho = sqrt((t11 - t21) * (t11 - t21) + (t12 - t22) * (t12 - t22));
//					if (rho > 1e-7) in = in + Ker(t11, t12, t21, t22)/* * phi2(t1, t2, ii, jj) * phi2(t1, t2, ii, jj)*/;
//				}
//			}
//		}
//	}
//
//	return in * h11 * h12 * h21 * h22;
//}


// f = 1 если первая матрица, если нет то вторая
complex middlepryam2(int i1, int j1, int i2, int j2, int f1, int f2) {
	double aa1, aa2, bb1, bb2, cc1, cc2, dd1, dd2;
	if (f1)
	{
		aa1 = GranA1 + i1 * H11; bb1 = aa1 + H11; cc1 = GranC1 + j1 * H12; dd1 = cc1 + H12;
	}
	else {
		aa1 = GranA2 + i1 * H21; bb1 = aa1 + H21; cc1 = GranC2 + j1 * H22; dd1 = cc1 + H22;
	};
	if (f2)
	{
		aa2 = GranA1 + i2 * H11; bb2 = aa1 + H11; cc2 = GranC1 + j2 * H12; dd2 = cc2 + H12;
	}
	else {
		aa2 = GranA2 + i2 * H21; bb2 = aa1 + H21; cc2 = GranC2 + j2 * H22; dd2 = cc2 + H22;
	};

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
					if (rho > 1e-7) in = in + Ker(t11, t12, t21, t22) * sqrtEGF2();
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

void Gauss(int k, complex Matrix[N][N + 1]) {
	if (Matrix[k][k] != complex(1.0, 0.0)) {
		complex T = Matrix[k][k];
		for (int j = k; j < N + 1; j++) {
			Matrix[k][j] = Matrix[k][j] / T;
		}
	}
	for (int i = 0; i < N; i++) {
		if ((Matrix[i][k] != complex(0.0, 0.0)) && (i != k)) {
			complex T = Matrix[i][k];
			Matrix[i][k] = complex(0.0, 0.0);
			for (int j = k + 1; j < N + 1; j++) {
				Matrix[i][j] -= Matrix[k][j] * T;
			}
		}
	}
	if (k < N - 1) {
		Gauss(k + 1, Matrix);
	}
}

complex un(double x1, double x2, int f) {
	complex s(0.0, 0.0);
	for (int i = 0; i < N; i++) {
		if (f)
		{
			s = s + C[i] * phi2(x1, x2, i / n, i % n, f); // or vice versa ... i%n, i/n)
		}
		else {
			s = s + C[i] * phi2(x1, x2, i / n, i % n, f); // or vice versa ... i%n, i/n)
		}
	}
	return(s);
}

void print_un(int pn, int f) {
	// double localA1, localB1, localC1, localD1;
	double t1, t2;
	if (f) {
		printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA1 + (GranB1 - GranA1) / pn * i1;
				t2 = GranC1 + (GranD1 - GranC1) / pn * i2;
				printf("%5.5f ", abs(un(t1, t2, f)));
				/// File1 << abs(un(t1, t2)) << " ";
			}
			printf("\n");
			// File1 << std::endl;
		}
	}
	else {
		printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA2 + (GranB2 - GranA2) / pn * i1;
				t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
				printf("%5.5f ", abs(un(t1, t2, f)));
				/// File1 << abs(un(t1, t2)) << " ";
			}
			printf("\n");
			// File1 << std::endl;
		}
	}

}

void Zapis_v_File(int pn, int f) {
	//double localA1, localB1, localC1, localD1;
	if (f) {
		std::ofstream File1("./Matrix_1.txt");
		FILE* tab_file;
		fopen_s(&tab_file, "result1.xls", "w");
		double t1, t2;
		// printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA1 + (GranB1 - GranA1) / pn * i1;
				t2 = GranC1 + (GranD1 - GranC1) / pn * i2;
				// printf("%6.3f ", abs(un(t1, t2)));
				File1 << abs(un(t1, t2, f)) << "\t";
				fprintf(tab_file, "%5.5f\t", abs(un(t1, t2, f)));
			}
			//printf("\n");
			File1 << "\n";
			fprintf(tab_file, "\n");
		}
		File1.close();
		fclose(tab_file);
	}
	else {
		std::ofstream File1("./Matrix_2.txt");
		FILE* tab_file;
		fopen_s(&tab_file, "result2.xls", "w");
		double t1, t2;
		// printf("\n");
		for (int i1 = 0; i1 < pn; i1++) {
			for (int i2 = 0; i2 < pn; i2++) {
				t1 = GranA2 + (GranB2 - GranA2) / pn * i1;
				t2 = GranC2 + (GranD2 - GranC2) / pn * i2;
				// printf("%6.3f ", abs(un(t1, t2)));
				File1 << abs(un(t1, t2, f)) << "\t";
				fprintf(tab_file, "%5.5f\t", abs(un(t1, t2, f)));
			}
			//printf("\n");
			File1 << "\n";
			fprintf(tab_file, "\n");
		}
		File1.close();
		fclose(tab_file);
	}
}


int main() {
	//cout << " AAAAAAAAAAAAAA";

	//double xi2[n + 1], xi1[n + 1];

	for (int j = 0; j < n + 1; j++) {
		X11[j] = GranA1 + j * H11;
		X12[j] = GranC1 + j * H12;
		cout << X11[j] << "   " << X12[j] << endl;
	}
	for (int j = 0; j < n + 1; j++) {
		X21[j] = GranA2 + j * H21;
		X22[j] = GranC2 + j * H22;
		cout << X21[j] << "   " << X22[j] << endl;
	}

	std::cout << " ----------------------------------------------------- " << std::endl;
	int i, j;
	// для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
	for (int i1 = 0; i1 < n; i1++) 
	{
		for (int j1 = 0; j1 < n; j1++)
		{
			for (int i2 = 0; i2 < n; i2++) 
			{
				for (int j2 = 0; j2 < n; j2++)
				{
					i = i1 + n * j1;
					j = i2 + n * j2;
					A[i][j] = middlepryam2(i1, j1, i2, j2, 1, 1);
					A[i][j + n * n] = middlepryam2(i1, j1, i2, j2, 1, 0);
					A[i + n * n][j] = middlepryam2(i1, j1, i2, j2, 0, 1);
					A[i + n * n][j + n * n] = middlepryam2(i1, j1, i2, j2, 0, 0);
				}
				A[i][N] = middlepryam1(i1, j1, 1);
				A[i + n * n][N] = middlepryam1(i1, j1, 0);
			}
		}
	}
	// i = 0, j = 0;
	//// для галеркина двухмерного шаг должен быть уже в квадрате, и дополнительные 2 интеграла должны быть на каждом элементе
	//for (int i1 = 0; i1 < n; i1++) {
	//    for (int j1 = 0; j1 < n; j1++) {
	//        for (int i2 = 0; i2 < n; i2++) {
	//            for (int j2 = 0; j2 < n; j2++) {
	//                i = i1 + n * j1;
	//                j = i2 + n * j2;
	//                //A[i][j] = del2(i1, j1, i2, j2) * h1 * h2 - lymda * middlepryam2(i1, j1, i2, j2);

	//                A2[i][j] = middlepryam2(i1, j1, i2, j2, 0);
	//                //printcomplex(A1[i][j]);
	//                A2[i][n * n] = middlepryam1(i1, j1, i2, j2, 0); // перенес сюда потому что не было в j2 вне цикла
	//                /*phi2(x1[i1], x2[j1], i1, j1)* phi2(x1[i2], x2[j2], i2, j2)*/    /*xi1[j], xi1[j + 1], xi1[i], xi1[i + 1], xi2[j], xi2[j + 1], xi2[i], xi2[i + 1],*/
	//            }

	//        }
	//    }
	//}

	//for (int i = 0; i < N; i++) {
	//    for (int j = 0; j < N + 1; j++) {
	//        printcomplex(A[i][j]);
	//    }
	//    cout << endl;
	//}
	//cout << endl;

	cout << " ----------------------------------------------------- " << endl;
	//for (int i = 0; i < n * n; i++) {
	//    for (int j = 0; j < n * n + 1; j++) {
	//        printcomplex(A2[i][j]);
	//    }
	//    cout << endl;
	//}
	cout << endl;
	Gauss(0, A);
	//Gauss(0, A2);
	cout << " ----------------------------------------------------- " << endl;
	//for (int i = 0; i < n * n; i++) {
	//    for (int j = 0; j < n * n + 1; j++) {
	//        printcomplex(A1[i][j]);
	//    }
	//    cout << endl;
	//}
	//cout << endl;
	cout << " ----------------------------------------------------- " << endl;

	for (int i = 0; i < N; i++)
	{
		C[i] = A[i][N];
		printcomplex(C[i]); cout << "  " << "1" << endl;
	}
	/*cout << " ----------------------------------------------------- " << endl;
	for (int i = 0; i < n * n; i++)
	{
		C2[i] = A2[i][n * n];
		printcomplex(C2[i]); cout << "  " << "1" << endl;
	}*/

	print_un(15, 1);
	Zapis_v_File(15, 1);
	//print_un(15, 0);
	//Zapis_v_File(15, 0);

	//system("pause");
	return 0;
}