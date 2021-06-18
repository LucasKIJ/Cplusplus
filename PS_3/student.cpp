#include "student.hpp"
#include <iostream>

student::student()
{
	name = "unspecified";
	college = "unspecified";
	status = "unspecified";
	degree = "unspecified";
}

void student::total_battels(){
		std::cout << "Total amount of money owed by student in Battels: " << library_fines + tuition_fees << std::endl;
}
