#ifndef STUDENTHEADERDEF
#define STUDENTHEADERDEF

#include <iostream>
#include <string>
class student{
	public:
		std::string name, college, degree;
		double library_fines, tuition_fees;

		student();
		virtual void TotalBattels();
		void SetStatus(std::string status);
		std::string GetStatus();

	protected:
		std::string mStatus;
};
#endif

