#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdlib>  
#include <cmath>    
using namespace std;


//Initializing the matrix
vector<vector<double>> initial(int n, int m) {
    vector<vector<double>> A(n, vector<double>(m));
    for (int i = 0; i < n; i++) {
        A[i][i] = 1.0; 
    }
    return A;
}
//Output of matrix
void out(const vector<vector<double>>& A, int str, int stolb) {
    cout << endl;
    for (int i = 0; i < str; i++) {
        for (int j = 0; j < stolb; j++) {
            cout << "x[" << i << "][" << j << "] = " << setw(15) << A[i][j];
        }
        cout << endl;
    }
    cout << endl;
}
//Matrix solution using Gauss Method
vector<double> gauss(vector<vector<double>>& matrix, int n, int m) {
    vector<double> output(n);

    for (int j = 0; j < n; j++) {
        double max = abs(matrix[j][j]);
        int coord_str = j;
        for (int t = j; t < n; t++) {
            if (abs(matrix[t][j]) > max) {
                max = abs(matrix[t][j]);
                coord_str = t;
            }
        }
        if (max == 0.0) {
            cout << "Error 1" << endl;
            exit(1);
        }

        if (coord_str != j) {
            swap(matrix[j], matrix[coord_str]);
        }
        double elem = matrix[j][j];
        for (int c = j; c < m; c++) {
            matrix[j][c] /= elem;
        }

        for (int i2 = j + 1; i2 < n; i2++) {
            elem = matrix[i2][j];
            for (int k = j; k < m; k++) {
                matrix[i2][k] -= elem * matrix[j][k];
            }
        }
    }

    output[n - 1] = matrix[n - 1][n];
    for (int i = n - 2; i >= 0; i--) {
        output[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++) {
            output[i] -= matrix[i][j] * output[j];
        }
    }

    return output;
}
//Solution of the symetric matrix using LDL^T factorization
vector<double> LDLTFactorization(vector<vector<double>>& matrix, int n, int m) {
    //Initializing matrix L and vector D
    vector<vector<double>> L(n, vector<double>(n));
    vector<double> D(n, 0.0);

    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        //Finding(vichislenie to be correct =)) diagonal element D
        for (int k = 0; k < j; k++) {
            sum += L[j][k] * L[j][k] * D[k];
        }
        D[j] = matrix[j][j] - sum;

        for (int i = j + 1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k] * D[k];
            }
            if (D[j] == 0.0) {
                cout << "Error 2" << endl;
                exit(1);
            }
            L[i][j] = (matrix[i][j] - sum) / D[j];
        }
    }
    //SLAU solution with bottom triangle matrix L
    vector<double> y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (matrix[i][m - 1] - sum) / L[i][i];
    }
    //SLAU solution with diagonal matrix D
    vector<double> z(n);
    for (int i = 0; i < n; i++) {
        z[i] = y[i] / D[i];
    }
    //SLAU solution with upper matrix L^T
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += L[j][i] * x[j];
        }
        x[i] = z[i] - sum;
    }

    return x;
}

void mult(const vector<vector<double>>& x, int n, int m, const vector<double>& y, vector<double>& mt) {
    for (int i = 0; i < n; i++) {
        mt[i] = 0.0;
        for (int j = 0; j < m; j++) {
            mt[i] += x[i][j] * y[j];
        }
    }
}

double calculateNorm(const vector<double>& vec) {
    double maxNorm = 0.0;
    for (double value : vec) {
        double absValue = abs(value);
        if (absValue > maxNorm) {
            maxNorm = absValue;
        }
    }
    return maxNorm;
}

int main() {
    int n=3, m=4;
    
    vector<vector<double>> matrix = initial(n, m + 1);
    vector<vector<double>> symetricMatrix = initial(n, m + 1);

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
    
    symetricMatrix[0][0] = 6;
    symetricMatrix[0][1] = 13;
    symetricMatrix[0][2] = -17;
    symetricMatrix[1][0] = 13;
    symetricMatrix[1][1] = 29;
    symetricMatrix[1][2] = -38;
    symetricMatrix[2][0] = -17;
    symetricMatrix[2][1] = -38;
    symetricMatrix[2][2] = 50;
    
    symetricMatrix[0][3] = 2;
    symetricMatrix[1][3] = 4;
    symetricMatrix[2][3] = -5;

    int method;
    cout << "Enter method:" << endl;
    cout << "1. Method of Gauss" << endl;
    cout << "2. LDL^T factorization" << endl;
    cin >> method;

     vector<vector<double>> copy(n, vector<double>(n));
    vector<double> rightHandSide(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            copy[i][j] = matrix[i][j];
        }
        rightHandSide[i] = matrix[i][n];
    }

    vector<double> solution;

    if (method == 1) {
        solution = gauss(matrix, n, m);

        cout << "Gauss solution:" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << solution[i] << endl;
        }
    } else if (method == 2) {
        solution = LDLTFactorization(symetricMatrix, n, m);

        cout << "LDL^T factorization solution" << endl;
        for (int i = 0; i < n; i++) {
            cout << "x[" << i << "] = " << solution[i] << endl;
        }
    } else {
        cout << "Error." << endl;
    }
    
    vector<double> residualVector(n);

    mult(copy, n, n, solution, residualVector);

    for (int i = 0; i < n; i++) {
        residualVector[i] -= rightHandSide[i];
    }

    cout << "Residual vector:" << endl;
    for (int i = 0; i < n; i++) {
        cout << residualVector[i] << setw(5) << " ";
    }
    cout << endl;
    
double residualNorm = calculateNorm(residualVector);
cout << "Norm of the Residual Vector: " << residualNorm << endl;

 for(int i = 0; i < n; i++){
    copy[i][3] = solution[i];
 }
  vector<double> auxiliarySolution = gauss(copy,n,m);
double relativeErrorAuxiliary = 0.0;

for (int i = 0; i < n; i++) {
    double originalValue = solution[i];
    double auxiliaryValue = auxiliarySolution[i];
    
    if (originalValue != 0.0) {
        double relative = abs(auxiliaryValue - originalValue) / abs(originalValue);
        if (relative > relativeErrorAuxiliary) {
            relativeErrorAuxiliary = relative;
       }
    }
}

cout << "Relative error for auxiliary system solution: " << relativeErrorAuxiliary << endl;
    return 0;
}
//Something wrong, but something also correct. ;-(
