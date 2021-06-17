#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

void change_val(int* pNum, int new_num){
    *pNum = new_num;
}

int main(){
    // Initialise pointer
    int* pNum;
    pNum = new int;

    int new_num;

    std::cout << "Enter an integer: " << std::endl;
    std::cin >> *pNum;

    std::cout << "The memory address is: " << pNum << std::endl;
    std::cout << "The value in memory is: " << *pNum << std::endl;


    std::cout << "Enter a new integer: " << std::endl;
    std::cin >> new_num;

    change_val(pNum, new_num);
    
    std::cout << "The memory address is: " << pNum << std::endl;
    std::cout << "The value in memory is: " << *pNum << std::endl;

    // Garbage collection
    delete pNum, new_num;

    return 0;
}