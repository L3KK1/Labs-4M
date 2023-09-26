//#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <clocale>
using namespace std;

double** initial(int n, int m)
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[m];
	return A;
}
void out(double** A, int str, int stolb)
{
	cout << endl;
	for (int i = 0; i < str; i++)
	{
		for (int j = 0; j < stolb; j++)
			cout << setw(10) << A[i][j];
		cout << endl;
	}
	cout << endl;
}

double* gauss(double** matrix, int n, int m)
{
	//prjamoj
	double elem;
	for (int j = 0; j < n; j++)
	{
		int max = 0;
		int coord_str = 0;
		for (int t = j; t < n; t++)
		{
			if (abs(matrix[t][j]) > max)
			{
				max = abs(matrix[t][j]); coord_str = t;
			}
		}
		if (max > abs(matrix[j][j]))
		{
			double* ptr = matrix[j];
			matrix[j] = matrix[coord_str];
			matrix[coord_str] = ptr;
		}
		elem = matrix[j][j];
		for (int c = j; c < m; c++)
		{
			matrix[j][c] /= elem;   //delenije stroki na elem
		}
        out(matrix, n,m);
		for (int i2 = j + 1; i2 < n; i2++)
		{
			elem = matrix[i2][j];
			for (int k = j; k < m; k++)
				matrix[i2][k] -= elem * matrix[j][k];
		}
		//		out(matrix, n, m);

	}
	//obratnyj
	double* xx = new double[m];
	xx[n - 1] = matrix[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		xx[i] = matrix[i][n];
		for (int j = i + 1; j < n; j++) xx[i] -= matrix[i][j] * xx[j];
	}
    cout << "xx: " << endl;
	for (int i = 0; i < n; i++)
		cout << xx[i] << " ";
	cout << endl;
	delete[] matrix;

	return xx;
}

void mult(double** x, int n, int m, double* y, double* mt)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
				mt[i] += x[i][j] * y[j];
}

int main()
{
    setlocale(LC_ALL, "rus");
	int n = 3, m = 3;
	
	m++;
	double** matrix;
	matrix = initial(n, m);

	matrix[0][0] = 21.547;
	matrix[0][1] = -95.510;
	matrix[0][2] = -96.121;
	matrix[1][0] = 10.223;
	matrix[1][1] = -91.065;
	matrix[1][2] = -7.343;
	matrix[2][0] = 51.218;
	matrix[2][1] = 12.264;
	matrix[2][2] = 86.457;

	matrix[0][3] = -49.930;
	matrix[1][3] = -12.465;
	matrix[2][3] = 60.812;

	double** copy = initial(n, n);
	double* svoboda = new double[n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			copy[i][j] = matrix[i][j];
			svoboda[i] = matrix[i][n];
		}

	cout << "Матрица:" << endl;
	out(matrix, n, m);
	double* resh = new double[n];
	resh = gauss(matrix, n, m);

/*	out(copy, n, n);

	for (int i = 0; i < n; i++)
		cout << svoboda[i] << " ";
	cout << endl;
*/

	double* nevjaz = new double[n];
	for (int i = 0; i < n; i++)
		nevjaz[i] = 0;

	mult(copy, n, n, resh, nevjaz);

	for (int i = 0; i < n; i++)
		nevjaz[i] -= svoboda[i];
    cout << "Невязка: " << endl;
	for (int i = 0; i < n; i++)
		cout << nevjaz[i] << " ";
	cout << endl;cout << endl;

	return 0;
}
