#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

double multiply(double* x, double* y, int n){
    // Vector dot product function
    double product = 0;
    for (int i = 0; i < n; i++){
        product = product + x[i] * y[i];
    }
    return product;
}

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
    int* sizeA;
    int* sizeB;
    int* sizeProd;
    sizeA = new int [2];
    sizeB = new int [2];
    sizeProd = new int [2];

    while (true){
        std::cout << "Enter dimension of first matrix (number of rows then columns):" << std::endl;
        std::cin >> sizeA[0];   
        std::cin >> sizeA[1];


        std::cout << "Enter dimension of second matrix (number of rows then columns):" << std::endl;
        std::cin >> sizeB[0];
        std::cin >> sizeB[1];

        if (sizeA[1] == sizeB[0]){
            break;
        }
        

        std::cout << "Number if rows in the first matrix must match number of columns in the second \nEnter then again." << std::endl;
    }

    // Declare size of product
    sizeProd[0] = sizeA[0];
    sizeProd[1] = sizeB[1];

    // Declare matrices
    double** A;
    double** B;
    double** prod;
    
    // Allocate memory
    A = new double* [sizeA[0]];
    B = new double* [sizeB[0]];
    prod = new double* [sizeA[0]];


    for (int i = 0; i < sizeA[0]; i++){
        A[i] = new double [sizeA[1]];
    }

    for (int i = 0; i < sizeB[0]; i++){
        B[i] = new double [sizeB[1]];
    }

    for (int i = 0; i < sizeA[0]; i++){
        prod[i] = new double [sizeB[1]];
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

    std::cout << "Enter the second matrix (size " << sizeB[0] << "x" << sizeB[1] << ") \n";
    // Input second matrix
    for (int i = 0; i < sizeB[0]; i++){
        for (int j = 0; j < sizeB[1]; j++){
             std::cout << "Enter the enter position [" << i << "," << j <<  "] \n";
             std::cin >> B[i][j];
        }
    }

    std::cout << "The matrix multiplaction A and B is " << std::endl;
    prod = multiply(A,B,sizeA,sizeB);
    print_mat(prod,sizeProd);
    
    
    // Garbage collection
    for (int i = 0; i < sizeA[1]; i++){
        delete[] A[i];
    } 

    for (int i = 0; i < sizeB[1]; i++){
        delete[] B[i];
    } 

    for (int i = 0; i < sizeB[1]; i++){
        delete[] prod[i];
    } 


    delete[] A;
    delete[] B;
    delete[] prod;
    delete sizeA, sizeB, sizeProd;

    return 0;
}