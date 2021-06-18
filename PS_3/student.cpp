#include "student.hpp"
#include <iostream>

student::student()
{
	name = "unspecified";
	college = "unspecified";
	mStatus = "unspecified";
	degree = "unspecified";
}

void student::TotalBattels(){
		std::cout << "Total amount of money owed by student in Battels: " << library_fines + tuition_fees << std::endl;
}

void student::SetStatus(std::string status){
	if (status == "graduate" || status == "undergraduate"){
		mStatus = status;
	}
	else{
		std::cout << "Invalid status" << std::endl;
	}
}

std::string student::GetStatus(){
	return mStatus;
}
