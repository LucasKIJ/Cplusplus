#include "graduate_student.hpp"
#include <iostream>

graduate_student::graduate_student() : student()
{
    mStatus = "graduate";
}

void graduate_student::TotalBattels(){
    std::cout << "Total amount of money owed by student in Battels: " << library_fines << std::endl;
}
