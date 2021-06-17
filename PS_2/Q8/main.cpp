#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

double norm(double* x, int n, int p = 2){
    double sum = 0;

    if (p < 1){
        p = 2;
    }

    for (int i = 0; i < n; i++){
        sum = sum + pow(fabs(x[i]),p);
    }
    sum = pow(sum, 1./p);
    return sum;
}

int main(){
    double* pX;
    int p;
    int n;

    std::cout << "How long is your array?" << std::endl;
    std::cin >> n;

    pX = new double [n];

    std::cout << "Enter values of the array" << std::endl;
    for (int i = 0; i < n; i++){
        std::cin >> pX[i];
    }

    std::cout << "What norm would you like to take?" <<std::endl;
    std::cin >> p;

    std::cout << "The norm is: " << norm(pX,n,p) << std::endl; 

    // Garbage collection
    delete[] pX;

    return 0;
}
