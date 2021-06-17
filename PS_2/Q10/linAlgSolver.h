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

double GetColMax(double** A, int* sizeA, int col){
    int i_max = 0;
    double max_val = 0;
    for (int i = col; i < sizeA[0]; i++){
        if (max_val < fabs(A[i][col])){
            max_val = fabs(A[i][col]);
            i_max = i;
        }
    }
    return i_max;
}

void swap_row(double** A, int* sizeA, int row1, int row2){
    double temp;
    for (int i = 0; i < sizeA[1]; i++){
        temp = A[row1][i];
        A[row1][i] = A[row2][i];
        A[row2][i] = temp;
    }
}

double* GaussianElimination(double** A, int* sizeA, double* b){


    double M;
    // Perform elimination
    for (int k = 0; k < sizeA[0]; k++){
        
        // Get pivot position
        int i_max = GetColMax(A ,sizeA,k);

        // No pivot in this column.
        if (A[i_max,k] == 0){
            continue;
        }

        // If pivot exists and not equal to current column
        if (i_max != k){
            swap_row(A, sizeA, k, i_max);
         
            // Swap vector entries
            double temp;
            temp = b[k];
            b[k] = b[i_max];
            b[i_max] = temp;
        }



        for (int i = k+1; i < sizeA[0]; i++){

            M = A[i][k] / A[k][k];
            for (int j = k; j < sizeA[0]; j++){
                
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
