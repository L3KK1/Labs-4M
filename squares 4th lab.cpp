#include <iostream>
#include <math.h>
#include <clocale>
#include <time.h>
using namespace std;
const int N = 7;
double *Gauss(double **Array_A, double *Array_B, int n);

double* initialX(){
    double* X = new double[N];
X[0] = 19.1; X[1] = 25; X[2] = 30.1; X[3] = 36; X[4] = 40; X[5] = 45.1; X[6] = 50;
return X;
}
double* initialY(){
    double* Y = new double[N];
Y[0] = 76.3; Y[1] = 77.8; Y[2] = 79.75; Y[3] = 80.8; Y[4] = 82.35; Y[5] = 83.9; Y[6] = 85;
return Y;
}

int main()
{
    setlocale(LC_ALL, "rus");
double* X; X = initialX();
double* Y; Y = initialY();
int sum = 0;
int m = 2;
double* POWERX = new double[2*m];
cout << "POWERX: "; 

for (int k = 0; k < 2 * m; k++)
{   
    POWERX[k] = 0;
    for (int i = 0; i < N; i++){
        POWERX[k] += pow(X[i], k+1);
    }

}

double** SUMX = new double*[m+1];
for (int i = 0; i < m+1; i++){
    SUMX[i] = new double[m+1];
}
for (int l = 0; l < m+1; l++){
    for (int j = 0; j < m+1; j++){
        if (j+l){
            SUMX[l][j] = POWERX[l+j-1];
        }
        else {
            SUMX[l][j] = N;
        }
    }
}
double* PRAW = new double [m+1];
for (int l = 0; l < m+1; l++){
    PRAW[l] = 0;
    for (int i = 0; i < N; i++ ){
        PRAW[l] += Y[i] * pow(X[i], l);
    }
}
double * a = Gauss(SUMX, PRAW, m+1);
double S2 = 0;
for (int i = 0; i < N; i++){
    double sum = Y[i];
    for (int j = 0; j < m+1; j++){
        sum -= a[j] * pow(X[i], j);
    }
    S2 += pow(sum, 2);
}
S2 /= N - m - 1;
double sigma = sqrt(S2);
cout << 'Коэффіценты a: ' << endl;
for (int i = 0; i < m + 1; i++){
    cout << a[i] << " ";
}
cout << endl; cout << "Cреднеквадратічное отклоненіе: " << sigma << endl;

return 0;

}

double *Gauss(double **Array_A, double *Array_B, int n)
{
	double *X = new double[n];  //ìàññèâ îòâåòîâ
	for (int k = 0; k < n; k++) // ïðÿìîé õîä
	{
		for (int i = k + 1; i < n; i++)
		{
			if (abs(Array_A[i][k]) > abs(Array_A[k][k]))
			{
				swap(Array_A[i], Array_A[k]);  //ìåíÿþòñÿ àäðåñà ò.ê. äâóìåðíûé ìàññèâ
				swap(Array_B[i], Array_B[k]);   //ìåíÿþòñÿ çíà÷åíèÿ
			}
		}
		double A_Main = Array_A[k][k];
		if (A_Main == 0)
		{
			cout << "error\n";
			system("pause");
			exit(0);
		}
		for (int i = k; i < n; i++)
		{
			Array_A[k][i] /= A_Main;
		}
		Array_B[k] /= A_Main;
		for (int i = k + 1; i < n; i++)
		{
			double s = Array_A[i][k];
			for (int j = k; j < n; j++)
			{
				Array_A[i][j] -= s * Array_A[k][j];
			}
			Array_B[i] -= s * Array_B[k];
		}
	}
	for (int k = n - 1; k >= 0; k--)  //îáðàòíûé õîä
	{
		X[k] = Array_B[k];
		for (int i = n - 1; i > k; i--)
		{
			X[k] -= Array_A[k][i] * X[i];
		}
	}
	return X;
}
