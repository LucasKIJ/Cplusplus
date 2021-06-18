#include "dtc_student.hpp"
#include <iostream>

dtc_student::dtc_student() : graduate_student()
{
    mStatus = "graduate";
}

void dtc_student::TotalBattels(){
    std::cout << "DTC students don't pay fees or library fines." << std::endl;
}
