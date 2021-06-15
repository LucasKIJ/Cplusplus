#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

void swap(double* pX, double* pY){
    double temp = *pX;
    *pX = *pY;
    *pY = temp;
}

int main(){
    // Initialise pointer
    double* pX;
    double* pY;
    pX = new double;
    pY = new double;

    std::cout << "Enter the first number: " << std::endl;
    std::cin >> *pX;
    
    std::cout << "Enter the second number: " << std::endl;
    std::cin >> *pY;

    std::cout << "The memory address is: " << pX << std::endl;
    std::cout << "The value in memory is: " << *pX << std::endl;

    std::cout << "The memory address is: " << pY << std::endl;
    std::cout << "The value in memory is: " << *pY << std::endl;

    std::cout << "Now the values are swapped" << std::endl;

    swap(pX,pY);
    
    std::cout << "The memory address is: " << pX << std::endl;
    std::cout << "The value in memory is: " << *pX << std::endl;

    std::cout << "The memory address is: " << pY << std::endl;
    std::cout << "The value in memory is: " << *pY << std::endl;

    // Garbage collection
    delete pX, pY;
    
    return 0;
}