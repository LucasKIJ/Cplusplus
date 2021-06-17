#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include "linAlgSolver.h"

int main(){
    int* sizeA;
    sizeA = new int [2];

    std::cout << "Enter dimension of first matrix (number of rows then columns):" << std::endl;
    std::cin >> sizeA[0];   
    std::cin >> sizeA[1];

    // Declare matrices
    double** A;
    
    // Allocate memory
    A = new double* [sizeA[0]];

    for (int i = 0; i < sizeA[0]; i++){
        A[i] = new double [sizeA[1]];
    }

    // Get matrix values
    // Input first matrix
    std::cout << "Enter the first matrix (size " << sizeA[0] << "x" << sizeA[1] << ") \n";
    for (int i = 0; i < sizeA[0]; i++){
        for (int j = 0; j < sizeA[1]; j++){
             std::cout << "Enter the enter position [" << i << "," << j <<  "] \n";
             std::cin >> A[i][j];
        }
    }

    // Declare vector b
    double* b;
    b = new double [sizeA[0]];

    std::cout << "Enter the vector b (size " << sizeA[0] << ")" << std::endl;
    // Input b
    for (int i = 0; i < sizeA[0]; i++){
             std::cout << "Enter the enter position [" << i << "] \n";
             std::cin >> b[i];
    }

    double* x = GaussianElimination(A, sizeA, b);

    std::cout << "The solution to the linear system is: " <<std::endl;
    for (int i = 0; i < sizeA[0]; i++){
        std::cout << "x(" << i << ") = " << x[i] <<std::endl;
    }

    // Garbage collection
    for (int i = 0; i < sizeA[1]; i++){
        delete[] A[i];
    } 
    delete[] A;
    delete[] b;
    delete[] x;
    delete sizeA;

    return 0;
}