default: use_student
all: default

student.o : student.cpp student.hpp
	g++ -c -O student.cpp

graduate_student.o : graduate_student.cpp graduate_student.hpp student.hpp
	g++ -c -O graduate_student.cpp

use_student : student.o use_student.cpp
	g++ -O -o use_student student.o use_student.cpp

clean:
	rm -f *.o *~ use_student