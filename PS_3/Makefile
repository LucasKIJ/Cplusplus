default: use_student
all: default

student.o : student.cpp student.hpp
	g++ -c -O student.cpp

graduate_student.o : graduate_student.cpp graduate_student.hpp student.hpp
	g++ -c -O graduate_student.cpp

dtc_student.o : dtc_student.cpp dtc_student.hpp graduate_student.hpp
	g++ -c -O dtc_student.cpp

use_student : student.o use_student.cpp
	g++ -O -o use_student student.o graduate_student.o dtc_student.o use_student.cpp

clean:
	rm -f *.o *~ use_student