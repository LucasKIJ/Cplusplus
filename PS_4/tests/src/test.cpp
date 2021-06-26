#include "test.hpp"
#include "Matrix.hpp"

// Function to add two integers 
int addIntegers(const int& a, const int& b)
{
  return a + b;
}


// Test case to check commutative property of addIntegers 
BEGIN_TEST(AddProps, Commutative)
{
  EXPECT_EQ(addIntegers(3, 5), addIntegers(5, 3))
}END_TEST

BEGIN_TEST(MatrixAdd, Commutative)
{
  Matrix A = rand(10);
  Matrix B = rand(10);

}END_TEST

int main()
{
  std::cout << "Test" <<std::endl;
  return 0;
}