#include <stdlib.h>
#include <iostream>
#include <cassert>
#include "Exception.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "MatrixTests.hpp"
#include "Examples.hpp"

int main(int argc, char *argv[]) {
  // This would  produce a compiler warning (there is no default constructor)
  try {
    // MatrixTests::RunTests();
    Examples::RunExamples();
  } catch (Exception &ex) {
    ex.DebugPrint();
  }

  /*
    Matrix A_matrix = Matrix(3,3);

    //Show that friends and methods can be used to do the same thing
    assert ( a_vector.norm() == norm(a_vector));
    assert ( a_vector.norm(3) == norm(a_vector, 3));


    Matrix B_matrix(3,3);
    B_matrix = A_matrix * A_matrix;
    print(B_matrix);

    std::cout << "a_vector = " << a_vector << "\n";
    std::cout << "a_vector = " << a_vector(1) << "\n";

    Vector bigger_vector(3);
    Vector smaller_vector(1);
    std::cout << "The following produces a warning\n";
    // This produces a warning
    bigger_vector = a_vector;
    std::cout << "The following throws an exception\n";
    // This throws an exception
    try
    {
        smaller_vector = a_vector;
    }
    catch (Exception &ex)
    {
        ex.DebugPrint();
    }
    std::cout << "The following throws an exception\n";
    // This throws an exception
    try
    {
        a_vector/2.0;//Okay
        a_vector/0.0;
    }
    catch (Exception &ex)
    {
        ex.DebugPrint();
    }
    exit(0);
    */
  return 0;
}
