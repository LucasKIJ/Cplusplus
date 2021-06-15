#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

double vec_mult(double* x1, double* x2, int n){
    int product = 0;
    for (int i = 0; i < n; i++){
        product = product + x1[i] * x2[i];
    }
    return product;
}

int main(){
    double* pX;
    int n;

    std::cout << "How long is your array?" << std::endl;
    std::cin >> n;

    pX = new double [n];

    std::cout << "Enter values of the array" << std::endl;
    for (int i = 0; i < n; i++){
        std::cin >> pX[i];
    }


    std::cout << "Its mean is: " << vec_mult(pX,pY,n) << std::endl;

    //Garbage collection
    delete[] pX;

    return 0;
}