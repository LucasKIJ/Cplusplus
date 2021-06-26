#include "MatrixTests.hpp"
namespace MatrixTests
{
void RunTests()
{
  Equality();
  //Add();
}

void Equality()
{
  Matrix A(3,3);
  Matrix B(3,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
      B.setValue(temp, i,j);
    }
  }
  print(A);
  print(B);

  if (&A== &B && &B== &A)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void Add()
{
  Matrix A(3,3);
  Matrix B(3,3);
  Matrix C(3,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
      B.setValue(temp*temp, i ,j);
      C.setValue(temp + temp*temp, i,j);
    }
  }
  if (C == A + B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}
}