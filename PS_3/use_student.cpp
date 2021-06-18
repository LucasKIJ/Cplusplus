#include "student.hpp"
#include "graduate_student.hpp"
#include "dtc_student.hpp"
#include <iostream>

int main(int argc, char* argv[])
{   
    std::cout << "Undergraduate Student test:" << std::endl;
    student joe;
    joe.tuition_fees = 9000;
    joe.library_fines = 10;
    joe.TotalBattels();
    joe.SetStatus("ggraduate");
    joe.SetStatus("undergraduate");
    std::cout << joe.GetStatus() << std::endl;


    std::cout << "Graduate Student test:" << std::endl;
    graduate_student alice;
    std::cout << alice.GetStatus() << std::endl;
    alice.tuition_fees = 9000;
    alice.library_fines = 10;
    alice.TotalBattels();

    std::cout << "DTC Student test:" << std::endl;
    dtc_student dan;
    dan.TotalBattels();

    std::cout << "Question 4 test:" << std::endl;
    student* alex = new dtc_student;
    alex -> TotalBattels();
    return 0;
}
