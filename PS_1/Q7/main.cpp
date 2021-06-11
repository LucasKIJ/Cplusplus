#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#define size 3

double vec_mult(double x1[], double x2[]){
    int product = 0;
    for (int i = 0; i < size; i++){
        product = product + x1[i] * x2[i];
    }
    return product;
}

int main(){
    double x1[size];
    double x2[size];
    std::cout << "Enter the first vector (size " << size << ") \n";
    std::cin >> x1[0];
    std::cin >> x1[1];
    std::cin >> x1[2];

    std::cout << "Enter the second vector (size " << size << ") \n";
    std::cin >> x2[0];
    std::cin >> x2[1];
    std::cin >> x2[2];

    double prod;
    prod = vec_mult(x1,x2);
    std::cout << "\nThe dot product is: " << prod;
}