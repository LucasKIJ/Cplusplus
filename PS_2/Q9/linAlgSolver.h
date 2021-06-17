#ifndef LINALGSOLVER_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define LINALGSOLVER_H
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

void print_mat(double** A, int* sizeA){
    std::cout << "\n"; 
    for (int i = 0; i < sizeA[0]; i++){
        for(int j=0; j < sizeA[1]; j++)
        {
            std::cout << A[i][j]<<" ";
        }
        std::cout << "\n";
    }
    std:: cout <<"\n";
}

double* GaussianElimination(double** A, int* sizeA, double* b){
    // Check strict assumption of gaussian elimination without pivots
    if (A[0][0] == 0){
        std::cout << "A(0,0) = 0. Cannot compute Gaussian Elimination without pivots";
        return 0;
    }

    double M;

    // Perform elimination
    for (int k = 0; k < sizeA[1]; k++){
        for (int i = k+1; i < sizeA[0]; i++){

            M = A[i][k] / A[k][k];
            for (int j = k; j < sizeA[1]; j++){
                
                A[i][j] = A[i][j] - M * A[k][j];
            }
            b[i] = b[i] - M * b[k];
        }
     }

    // Declare x
    double* x;
    x = new double [sizeA[1]];

    for (int k = sizeA[1] - 1; k >= 0; k--){
        double sum = 0;
        for (int j = k; j < sizeA[1]; j++){
            sum += A[k][j] * x[j];
        }
        x[k] = (1. / A[k][k]) * (b[k] - sum);
    }

    return x;
}
#endif
