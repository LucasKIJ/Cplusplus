#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>


int main(){
    // Initialise pointer
    int* pNum;
    pNum = new int;

    std::cout << "Enter an integer: " << std::endl;
    std::cin >> *pNum;

    std::cout << "The memory address is: " << pNum << std::endl;
    std::cout << "The value in memory is: " << *pNum << std::endl;

    return 0;
}