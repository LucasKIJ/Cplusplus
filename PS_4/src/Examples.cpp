#include "Examples.hpp"
#include <tuple>

namespace Examples {
void RunExamples() {
  std::cout << " ---- Running All Examples ---- \n" << std::endl;
  GaussianElimination();
  LUDecomposition();
  Det();
  QRFactorisation();
  LeastSquares();
  HessReduction();
  EigenVal();
}

void GaussianElimination() {

  /* Gaussian elimination is used to solve full rank systems
  of the form:  Ax = b.
  For examples we may wish to solve the following system:
       2x + y -  z =   8
      -3x - y + 2z = -11
      -2x + y + 2z =  -3
  In matrix notation we write this as:
  [ 2   1  -1] [x]   [  8]
  [-3  -1   2] [y] = [-11]
  [-2   1   2] [z]   [ -3]
  Now we define Matrix A and vector b using our class
  */

  // Declare matrix A
  Matrix A(3, 3);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;

  // Declare vector b
  Vector b(3);
  b(1) = 8;
  b(2) = -11;
  b(3) = -3;

  // Solve using Gaussion Elimination
  Vector x(3);
  x = gaussianElimination(A, b);

  std::cout << "\n ---- Gaussian Elimination Example ---- \n" << std::endl;

  std::cout << "Gaussian elimination is used to solve full rank systems\n"
            << "of the form:  Ax = b.\n"
            << "For examples we may wish to solve the following system:\n"
            << "     2x + y -  z =   8\n"
            << "    -3x - y + 2z = -11\n"
            << "    -2x + y - 2z =  -3\n"
            << "In matrix notation we write this as:\n"
            << "[ 2   1  -1] [x]   [  8]\n"
            << "[-3  -1   2] [y] = [-11]\n"
            << "[-2   1   2] [z]   [ -3]\n"
            << "Now we define Matrix A and vector b using our class\n \n"
            << "Code: \n"
            << "// Declare matrix A\n"
            << "Matrix A(3,3);\n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2;\n \n"
            << "// Declare vector b\n"
            << "Vector b(3);\n"
            << "b(1) =   8;\n"
            << "b(2) = -11;\n"
            << "b(3) =  -3;\n \n"
            << "// Solve using Gaussion Elimination\n"
            << "Vector x(3);\n"
            << "x = gaussianElimination(A,b);\n"
            << "print(x);\n"
            << "x = \n" << x << std::endl;
}

void LUDecomposition() {

  /* Find the LU Decomposition of a square matrix A.
  Let
          [ 2   1  -1]
   A =    [-3  -1   2]
          [-2   1   2]
  */

  // Declare matrix A
  Matrix A(3, 3);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;

  // Find PLU decompositon
  int swaps;
  Matrix P(3, 3);
  Matrix L(3, 3);
  Matrix U(3, 3);
  std::tie(P, L, U, swaps) = lu(A);

  std::cout << "\n ---- LU Decomposition Example ---- \n" << std::endl;

  std::cout << "Find the LU Decomposition of a square matrix A.\n"
            << "Let\n"
            << "       [ 2   1  -1] \n"
            << "A =    [-3  -1   2] \n"
            << "       [-2   1   2] \n \n"
            << "Code : \n"
            << "// Declare matrix A \n"
            << "Matrix A(3,3); \n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2;\n \n"
            << "// Find PLU decompositon \n"
            << "int swaps; \n"
            << "Matrix P(3,3); Matrix L(3,3); Matrix U(3,3);\n"
            << "std::tie(P,L,U,swaps) = lu(A);\n \n"
            << "We have then that\n"
            << "P = \n" << P << "\n L = \n" << L << "\n U = \n" << U
            << std::endl;
}

void Det() {

  /* Find the  Determinat of a square matrix A.
  Let
          [ 2   1  -1]
   A =    [-3  -1   2]
          [-2   1   2]
  */

  // Declare matrix A
  Matrix A(3, 3);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;

  // Find Determinant
  double D = det(A);

  std::cout << "\n ---- Determinant Example ---- \n" << std::endl;

  std::cout << "Find the determinant of a square matrix A.\n"
            << "Let\n"
            << "       [ 2   1  -1] \n"
            << "A =    [-3  -1   2] \n"
            << "       [-2   1   2] \n \n"
            << "Code : \n"
            << "// Declare matrix A \n"
            << "Matrix A(3,3); \n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2;\n \n"
            << "// Find determinant \n"
            << "double D = det(A); \n \n"
            << "Determinant of A = " << D << std::endl;
}

void QRFactorisation() {

  /* Find the QR Factorisation of  matrix A.
  Let
          [ 2   1  -1]
   A =    [-3  -1   2]
          [-2   1   2]
          [ 5   7  13]
  */

  // Declare matrix A
  Matrix A(4, 3);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;
  A(4, 1) = -5;
  A(4, 2) = 7;
  A(4, 3) = 13;

  // Find QR factorisation
  Matrix Q(4, 3);
  Matrix R(4, 3);
  Matrix Q_reduced(4, 3);
  Matrix R_reduced(3, 3);
  std::tie(Q, R) = qr(A);
  std::tie(Q_reduced, R_reduced) = qr(A, true);

  std::cout << "\n ---- QR Factorisation Example ---- \n" << std::endl;

  std::cout << "Find the QR Factorisation of matrix A.\n"
            << "Let\n"
            << "       [ 2   1  -1] \n"
            << "A =    [-3  -1   2] \n"
            << "       [-2   1   2] \n"
            << "       [ 5   6  13] \n \n"
            << "Code : \n"
            << "// Declare matrix A \n"
            << "Matrix A(3,3); \n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2;\n"
            << "A(4,1) = -5; A(4,2) =  7; A(4,3) = 13;\n \n"
            << "// Find QR Factorisation \n"
            << "Matrix Q(4,3); Matrix R(4,3);\n"
            << "Matrix Q_reduced(4,3); Matrix R_reduced(3,3);\n"
            << "std::tie(Q,R) = qr(A);\n"
            << "std::tie(Q_reduced,R_reduced) = qr(A, true);\n \n"
            << "The full QR Factorisation is \n"
            << "Q = \n" << Q << "\n R = \n" << R << std::endl
            << "The reduce QR Factorisation is \n"
            << "Q = \n" << Q_reduced << "\n R = \n" << R_reduced << "\n"
            << std::endl;
}

void LeastSquares() {

  /* Find the solution to the least squares problem
       min ||A x = b||_2
  where
          [ 2   1  -1] [x1]   [1]
          [-3  -1   2] [x2]   [2]
          [-2   1   2] [x3] = [3]
          [ 5   7  13] [x4]   [4]
          [10  20  30] [x5]   [5]
  */

  // Declare matrix A
  Matrix A(5, 3);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;
  A(4, 1) = -5;
  A(4, 2) = 7;
  A(4, 3) = 13;
  A(5, 1) = 10;
  A(5, 2) = 20;
  A(5, 3) = 30;

  // Declare vector b
  Vector b(5);
  b(1) = 1;
  b(2) = 2;
  b(3) = 3;
  b(4) = 4;
  b(5) = 5;

  // Find Least squares solution
  Vector x(3);
  x = lsq(A, b);

  std::cout << "\n ---- Least squares Example ---- \n" << std::endl;

  std::cout << "Find the solution to the least squares\n"
            << "problem  min ||Ax-b||_2"
            << "Let\n"
            << "       [ 2   1  -1] [x1]   [1]\n"
            << "A =    [-3  -1   2] [x2]   [2]\n"
            << "       [-2   1   2] [x3] = [3]\n"
            << "       [ 5   7  13] [x4]   [4]\n"
            << "       [10  20  30] [x5]   [5]\n \n"
            << "Code :\n"
            << "// Declare matrix A \n"
            << "Matrix A(3,3); \n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2;\n"
            << "A(4,1) = -5; A(4,2) =  7; A(4,3) = 13;\n"
            << "A(5,1) = 10; A(5,2) = 20; A(5,3) = 30;\n \n"
            << "// Declare vector b"
            << "Vector b(5);\n"
            << "b(1) = 1; b(2) = 2; b(3) = 3; b(4) = 4; b(5) = 5;\n \n"
            << "// Find Least squares solution\n"
            << "Vector x(3); \n"
            << "x = lsq(A);\n \n"
            << "The solution to the least squares problem is \n"
            << "\n x = \n" << x << std::endl;
}

void HessReduction() {

  /* Find the upper Hessenberg reduction of a square matrix A.
  Let
          [ 2   1  -1   3]
      A = [-3  -1   2   4]
          [-2   1   2   5]
          [-6   2   8   9]
  */

  // Declare matrix A
  Matrix A(4, 4);
  A(1, 1) = 2;
  A(1, 2) = 1;
  A(1, 3) = -1;
  A(1, 4) = 3;
  A(2, 1) = -3;
  A(2, 2) = -1;
  A(2, 3) = 2;
  A(2, 4) = 4;
  A(3, 1) = -2;
  A(3, 2) = 1;
  A(3, 3) = 2;
  A(3, 4) = 5;
  A(4, 1) = -6;
  A(4, 2) = 2;
  A(4, 3) = 8;
  A(4, 4) = 9;

  // Find upper hessenberg Reduction
  Matrix H(4, 4);
  H = hessenbergReduction(A);

  std::cout << "\n ---- Hessenberg Reduction Example ---- \n" << std::endl;

  std::cout << "Find the upper Hessenberg Reduction of a square matrix A.\n"
            << "Let\n"
            << "       [ 2   1  -1  3] \n"
            << "A =    [-3  -1   2  4] \n"
            << "       [-2   1   2  5] \n"
            << "       [-6   2   8  9] \n \n"
            << "Code :\n"
            << "// Declare matrix A \n"
            << "Matrix A(4,4); \n"
            << "A(1,1) =  2; A(1,2) =  1; A(1,3) = -1; A(1,4) = 3;\n"
            << "A(2,1) = -3; A(2,2) = -1; A(2,3) =  2; A(2,4) = 4;\n"
            << "A(3,1) = -2; A(3,2) =  1; A(3,3) =  2; A(3,4) = 5;\n"
            << "A(4,1) = -6; A(4,2) =  2; A(4,3) =  8; A(4,4) = 9;\n \n"
            << "// Find upper Hessenberg Reduction \n"
            << "Matrix H(3,3); \n"
            << "H = hessenbergReduction(A); \n \n"
            << "The upper hessenberg Reduction of A \n \n H = \n" << H
            << std::endl;
}

void EigenVal() {
  /* Find the eigenvalues of a square positive definite matrix A.
  Let's find the eigenvalues of the Magic matrix size 5
          [ 17    24     1     8    15 ]
          [ 23     5     7    14    16 ]
      A = [  4     6    13    20    22 ]
          [ 10    12    19    21     3 ]
          [ 11    18    25     2     9 ]
  */

  // Declare matrix A
  Matrix A(5, 5);
  A(1, 1) = 17;
  A(1, 2) = 24;
  A(1, 3) = 1;
  A(1, 4) = 8;
  A(1, 5) = 15;
  A(2, 1) = 23;
  A(2, 2) = 5;
  A(2, 3) = 7;
  A(2, 4) = 14;
  A(2, 5) = 16;
  A(3, 1) = 4;
  A(3, 2) = 6;
  A(3, 3) = 13;
  A(3, 4) = 20;
  A(3, 5) = 22;
  A(4, 1) = 10;
  A(4, 2) = 12;
  A(4, 3) = 19;
  A(4, 4) = 21;
  A(4, 5) = 3;
  A(5, 1) = 11;
  A(5, 2) = 18;
  A(5, 3) = 25;
  A(5, 4) = 2;
  A(5, 5) = 9;

  // Find eigenvalues
  Vector eigs(5);
  eigs = eigenVal(A);

  std::cout << "\n ---- Eigenvalue Example ---- \n" << std::endl;

  std::cout
      << "Find the eigenvalues of a square positive definite matrix A.\n"
      << "Let's find the eigenvalues of the Magic matrix size 5\n"
      << "       [ 17  24   1   8  15 ] \n"
      << "A =    [ 23   5   7  14  16 ] \n"
      << "       [  4   6  13  20  22 ] \n"
      << "       [ 10  12  19  21   3 ]"
      << "       [ 11  18  25   2   9 ]\n \n"
      << "Code :\n"
      << "// Declare matrix A \n"
      << "Matrix A(5,5); \n"
      << "A(1,1) = 17; A(1,2) = 24; A(1,3) =  1; A(1,4) =  8; A(1,5) = 15;\n"
      << "A(2,1) = 23; A(2,2) =  5; A(2,3) =  7; A(2,4) = 14; A(2,5) = 16;\n"
      << "A(3,1) =  4; A(3,2) =  6; A(3,3) = 13; A(3,4) = 20; A(3,5) = 22;\n"
      << "A(4,1) = 10; A(4,2) = 12; A(4,3) = 19; A(4,4) = 21; A(4,5) =  3;\n"
      << "A(5,1) = 11; A(5,2) = 18; A(5,3) = 25; A(5,4) =  2; A(5,5) =  9;\n "
         "\n "
      << "// Find eigenvalues \n"
      << "Vector eigs(5); \n"
      << "eigs = eigenVal(A); \n \n"
      << "The eigenvalues of A \n \n eigs = \n" << eigs << std::endl;
}
}