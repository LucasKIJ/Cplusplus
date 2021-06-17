#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

double vec_mult(double* x1, double* x2, int n){
    double product = 0;
    for (int i = 0; i < n; i++){
        product = product + x1[i] * x2[i];
    }
    return product;
}

int main(){
    double* pX;
    double* pY;
    int n;

    std::cout << "How long is your vector?" << std::endl;
    std::cin >> n;

    pX = new double [n];
    pY = new double [n];

    std::cout << "Enter values the first vector" << std::endl;
    for (int i = 0; i < n; i++){
        std::cin >> pX[i];
    }

    std::cout << "Enter values the second vector" << std::endl;
    for (int i = 0; i < n; i++){
        std::cin >> pY[i];
    }

    std::cout << "Their product is: " << vec_mult(pX,pY,n) << std::endl;

    //Garbage collection
    delete[] pX;
    delete[] pY;

    return 0;
}