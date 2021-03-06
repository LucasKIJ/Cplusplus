#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

double mean(double* x, int n){
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum = sum + x[i];
    }
    
    return sum / n;
}

double sd(double* x, int n){
    double xbar = mean(x,n);
    double sum = 0;
    for (int i = 0; i < n; i++){
        sum = sum + pow((x[i] - xbar),2);
    }
    
    return pow(sum / (n-1),0.5);
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

    
    std::cout << "Its mean is: " << mean(pX,n) << std::endl;
    std::cout << "Its standard deviation is: " << sd(pX,n) << std::endl;
    
    //Garbage collection
    delete[] pX;

    return 0;
}