#include <iostream>
#include <cmath>
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
        out(matrix, n, m);
        for (int i2 = j + 1; i2 < n; i2++)
        {
            elem = matrix[i2][j];
            for (int k = j; k < m; k++)
                matrix[i2][k] -= elem * matrix[j][k];
        }
        //out(matrix, n, m);
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

double* auxiliary_system(double** matrix, int n, int m) {
    // Solve the auxiliary system Ax = Ax.
    double* x = gauss(matrix, n, m);

    // Compute the relative error of the solution.
    double* relative_error = new double[n];
    for (int i = 0; i < n; i++) {
        relative_error[i] = abs(x[i] - 1) / abs(x[i]);
    }

    return relative_error;
}

void mult(double** x, int n, int m, double* y, double* mt)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            mt[i] += x[i][j] * y[j];
}

void solve_auxiliary_system(double **A, double *b, int n) {
    double *x = new double[n];
    for (int i = 0; i < n; i++) {
        x[i] = 1.0;
    }
    gauss(A, n, n);
    for (int i = 0; i < n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;
}

double calculate_residual(double *x, double *x_aux, int n) {
    double delta = 0.0;
    for (int i = 0; i < n; i++) {
        delta = max(delta, abs(x[i] - x_aux[i]));
    }
    return delta / max(abs(x_aux[0]), 1e-10);
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

    double* discrepancy = new double[n];
    for (int i = 0; i < n; i++)
        discrepancy[i] = 0;

    mult(copy, n, n, resh, discrepancy);

    for (int i = 0; i < n; i++)
        discrepancy[i] -= svoboda[i];
    cout << "Невязка: " << endl;
    for (int i = 0; i < n; i++)
        cout << discrepancy[i] << " ";
    cout << endl;
    cout << endl;

    double **A = initial(n, n);
    double *b = new double[n];
    for (int i = 0; i < n; i++) {
        b[i] = discrepancy[i];
        for (int j = 0; j < n; j++) {
            A[i][j] = copy[i][j];
        }
    }

    double *x_aux = auxiliary_system(A, n, n);
    cout << "Решение вспомогательной системы: " << endl;
    for (int i = 0; i < n; i++) {
        cout << x_aux[i] << " ";
    }
    cout << endl;

    double delta = calculate_residual(resh, x_aux, n);
    cout << "Относительная погрешность δ: " << delta << endl;

    // Clean up dynamically allocated memory
    for (int i = 0; i < n; i++) {
        delete[] A[i];
        delete[] copy[i];
    }
    delete[] A;
    delete[] copy;
    delete[] b;
    delete[] resh;
    delete[] discrepancy;
    delete[] x_aux;

    return 0;
}
