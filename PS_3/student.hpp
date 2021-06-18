#ifndef STUDENTHEADERDEF
#define STUDENTHEADERDEF

#include <iostream>
#include <string>
class student{
	public:
		std::string name, college, status, degree;
		double library_fines, tuition_fees;

		student();
		void total_battels();
};
#endif

