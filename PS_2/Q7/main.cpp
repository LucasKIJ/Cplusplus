#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>


double* get_row(double** A, int i, int j, int n){
    // Get 
    double* row;
    row = new double [n];
    for (int j = 0; j < n; j++){
        row[j] = A[i][j];
    }
    return row;
}

double* get_col(double** A, int i, int j, int n){
    double* col;
    col = new double [n];
    for (int i = 0; i < n; i++){
        col[i] = A[i][j];
    }
    return col;


}

// Multiply

// Vector dot procduct
double multiply(double* x, double* y, int n){
    // Vector dot product function
    double product = 0;
    for (int i = 0; i < n; i++){
        product = product + x[i] * y[i];
    }
    return product;
}

// Matrix product
double** multiply(double** A, double** B, int* sizeA, int* sizeB){
    double* tempRow;
    double* tempCol;
    double** prod;
    prod = new double* [sizeA[0]];
    for (int i = 0; i < sizeA[0]; i++){
        prod[i] = new double [sizeB[1]];
    }

    for (int i = 0; i < sizeA[0]; i++){
        for (int j = 0; j < sizeB[1]; j++){
            
            tempRow = get_row(A,i,j,sizeA[1]);
            tempCol = get_col(B,i,j,sizeB[0]);
            prod[i][j] = multiply(tempRow,tempCol, sizeA[1]);
        }
    }

    // Garbage collection
    delete tempRow, tempCol;

    return prod;
}
// Matrix scalar product
double** multiply(double** A, double a, int* sizeA){
    double** prod;
    prod = new double* [sizeA[0]];
    for (int i = 0; i < sizeA[0]; i++){
        prod[i] = new double [sizeA[1]];
    }

    for (int i = 0; i < sizeA[0]; i++){    
        for (int j = 0; j < sizeA[1]; j++){
            prod[i][j] = A[i][j] * a;
        }
    }
    return prod;
}

double** multiply(double a, double** A, int* sizeA){
    return multiply(A, a, sizeA);
}

// Vector scalar product
double* multiply(double* x, double a, int sizeX){
    double* prod;
    prod = new double [sizeX];
    for (int i = 0; i < sizeX; i++){
        prod[i] = x[i] * a;
    }
    return prod;
}

double* multiply(double a, double* x, int sizeX){
    return multiply(x, a, sizeX);
}

// Matrix vector product
double* multiply(double** A, double* x, int* sizeA){
    double* prod;
    double* tempRow;
    prod = new double [sizeA[0]];

    for (int i = 0; i < sizeA[0]; i++){    
        tempRow = get_row(A,i,1,sizeA[1]);
        prod[i] = multiply(tempRow, x, sizeA[0]);
    }

    // Garbage collection
    delete tempRow;
    return prod;
}

double* multiply(double* x, double** A, int* sizeA){
    double* prod;
    double* tempCol;
    prod = new double [sizeA[1]];

    for (int j = 0; j < sizeA[1]; j++){    
        tempCol = get_col(A,1,j,sizeA[0]);
        prod[j] = multiply(tempCol, x, sizeA[0]);
    }

    // Garbage collection
    delete tempCol;
    return prod;
}


void print_mat(double** A, int* sizeA){

    for(int i=0; i < sizeA[0]; i++)
    {
        for(int j=0; j < sizeA[1]; j++)
        {
            std::cout << A[i][j]<<" ";
        }
        std::cout << "\n";
    }
    std:: cout <<"\n";
}

// ------------ Main --------------- //

int main(){
    return 0;
}