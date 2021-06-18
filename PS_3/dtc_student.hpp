#ifndef DTC_STUDENTHEADERDEF
#define DTC_STUDENTHEADERDEF

#include <iostream>
#include <string>
#include "graduate_student.hpp"

class dtc_student : public graduate_student{
	public:
		dtc_student();
		void TotalBattels();
};
#endif