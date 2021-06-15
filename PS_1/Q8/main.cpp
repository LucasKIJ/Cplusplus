  
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#define size 3

double vec_mult(double* x1, double* x2){
    // Vector dot product function
    double product = 0;
    for (int i = 0; i < size; i++){
        product = product + x1[i] * x2[i];
    }
    return product;
}

double* get_row(double** A, int i, int j){
    // Get 
    double* row;
    row = new double [size];
    for (int j = 0; j < size; j++){
        row[j] = A[i][j];
    }
    return row;
}

double* get_col(double** A, int i, int j){
    double* col;
    col = new double [size];
    for (int i = 0; i < size; i++){
        col[i] = A[i][j];
    }
    return col;
}


double** mat_mult(double** A, double** B){
    double* temp1;
    double* temp2;
    double** prod;
    prod = new double* [size];
    for (int i = 0; i < size; i++){
        prod[i] = new double [size];
    }

    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            
            temp1 = get_row(A,i,j);
            temp2 = get_col(B,i,j);
            prod[i][j] = vec_mult(temp1,temp2);
        }
    }

    // Garbage collection
    delete temp1, temp2;

    return prod;
}


void print_mat(double** A){

    for(int i=0; i < size; i++)
    {
        for(int j=0; j < size; j++)
        {
            std::cout << A[i][j]<<" ";
        }
        std::cout << "\n";
    }
    std:: cout <<"\n";
}

// *********** Main *********** //

int main(){
    // Define variables
    double** A;
    double** B;
    double** prod;

    // Initialise memory
    A = new double* [size];
    B = new double* [size];
    for (int i = 0; i < size; i++){
        A[i] = new double [size];
        B[i] = new double [size];
    }

    std::cout << "Enter the first matrix (size " << size << "x" << size << ") \n";
    // Input first matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
             std::cout << "Enter the enter position [" << i << "," << j <<  "] \n";
             std::cin >> A[i][j];
        }
    }

    std::cout << "Enter the second matrix (size " << size << "x" << size << ") \n";
    // Input second matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
             std::cout << "Enter the enter position [" << i << "," << j <<  "] \n";
             std::cin >> B[i][j];
        }
    }
    
    std::cout << "The first array you entered is: \n";
    print_mat(A);

    std::cout << "The second array you entered is: \n";
    print_mat(B);

    prod = mat_mult(A,B);
    std::cout << "\nThe matrix product is: \n";
    print_mat(prod);

    // Garbage collection
    for (int i = 0; i < size; i++){
        delete[] A[i];
        delete[] B[i];
        delete[] prod[i];
    }   

    delete[] A;
    delete[] B;
    delete[] prod;

    return 0;
}