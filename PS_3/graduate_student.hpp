#ifndef GRADUATE_STUDENTHEADERDEF
#define GRADUATE_STUDENTHEADERDEF

#include <iostream>
#include <string>
#include "student.hpp"

class graduate_student : public student{
	public:
		graduate_student();
		void TotalBattels();
};
#endif
