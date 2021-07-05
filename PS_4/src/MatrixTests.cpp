#include "MatrixTests.hpp"
#include <tuple>
namespace MatrixTests
{
void RunTests()
{
  Equality();
  Add();
  Minus();
  MultiplySquare();
  MultiplyRect();
  MultiplyConst();
  UnaryMinus();
  Transpose();
  RowMax();
  ColMax();
  SwapRow();
  SwapCol();
  Swap();
  Eye();
  LU();
  Det();
  //QR();
  //LSQ();
  //HessReduction();
  //EigenVal();
  //NormAll();
  //GaussianElimination();
}

void Equality()
{
  // Test matrix equality

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

  if (A==B && B==A)
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
  // Test matrix add
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
  if (C == A + B && C == B + A)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void Minus()
{
  // Test matrix minus
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
      C.setValue(temp - temp*temp, i,j);
    }
  }
  if (C == A - B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void MultiplySquare()
{
  // Test matrix multiply with square
  Matrix A(3,3);
  Matrix B(3,3);
  Matrix C(3,3);
  Matrix D(3,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
      B.setValue(temp*temp, i ,j);
    }
  }
  C(1,1) = 180; C(1,2) = 246; C(1,3) = 324;
  C(2,1) = 378; C(2,2) = 525; C(2,3) = 702;
  C(3,1) = 576; C(3,2) = 804; C(3,3) = 1080;

  D(1,1) =  80; D(1,2) =  94; D(1,3) = 108;
  D(2,1) = 368; D(2,2) = 445; D(2,3) = 522;
  D(3,1) = 872; D(3,2) =1066; D(3,3) = 1260;
  if (C == A * B && D == B * A)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void MultiplyRect()
{
  // Test matrix multiply with rectangular matrices
  Matrix A(3,2);
  Matrix B(2,3);
  Matrix C(3,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }

  temp = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      temp += 1;
      B.setValue(temp*temp, i ,j);
    }
  }
  C(1,1) =  33; C(1,2) =  54; C(1,3) =  81;
  C(2,1) =  67; C(2,2) = 112; C(2,3) = 171;
  C(3,1) = 101; C(3,2) = 170; C(3,3) = 261;
  if (C == A * B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void MultiplyConst()
{
  // Test double times matrix
  Matrix A(3,2);
  Matrix B(3,2);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp,  i,j);
      B.setValue(temp*5, i,j);
    }
  }
  
  if (A*5 == B && 5*A == B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}


void UnaryMinus()
{
  // Test unary operator
  Matrix A(3,2);
  Matrix B(3,2);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp,  i,j);
      B.setValue(-temp, i,j);
    }
  }
  
  if (A == -B && -A == B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void Transpose()
{
  // Test transpose
  Matrix A(3,2);
  Matrix B(2,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }

  temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      B.setValue(temp, j ,i);
    }
  }
  if (A == B.T() && A == transpose(B))
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void RowMax()
{
  // Test row max
  Matrix A(3,2);
  Matrix B(2,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }

  temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      B.setValue(-temp, j ,i);
    }
  }
  if (A.rowMax(1)==4 && B.argRowMax(0,1) == 2 && B.rowMax(0,1,false) == -3)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void ColMax()
{
  // Test col max
  Matrix A(3,2);
  Matrix B(2,3);
  double temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }

  temp = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      B.setValue(-temp, j ,i);
    }
  }
  if (A.colMax(1)==6 && B.argColMax(0,1) == 1 && B.colMax(0,1,false) == -2)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void SwapRow()
{
  // Test swapping rows
  Matrix A(2,2);
  Matrix B(2,2);
  double temp = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }
  A.swapRow(0,1);
  B(1,1) = 3; B(1,2) = 4;
  B(2,1) = 1; B(2,2) = 2;

  if (A == B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void SwapCol()
{
  // test swapping columns
  Matrix A(2,2);
  Matrix B(2,2);
  double temp = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }
  A.swapCol(0,1);
  B(1,1) = 2; B(1,2) = 1;
  B(2,1) = 4; B(2,2) = 3;

  if (A == B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}


void Eye()
{
  // Test identity construction
  // Check tall
  Matrix A = eye(3,2);
  Matrix B(3,2);
  // Check fat
  Matrix C = eye(2,3);
  Matrix D(2,3);
  // Check square
  Matrix E = eye(2);
  Matrix F(2);
  // Check add
  Matrix G(3,2);

  B(1,1) = 1; B(2,2) = 1;
  D(1,1) = 1; D(2,2) = 1;
  F(1,1) = 1; F(2,2) = 1;

  if (A == B && C == D && E == F && G + eye(3,2) == eye(3,2))
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void Swap()
{
  // test swap function
  Matrix A(2,2);
  Matrix B(2,2);
  double temp = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      temp += 1;
      A.setValue(temp, i,j);
    }
  }
  A.swap(0,0,0,1);
  B(1,1) = 2; B(1,2) = 1;
  B(2,1) = 3; B(2,2) = 4;

  if (A == B)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}

void LU()
{
  // Test LU decomposition
  Matrix A(3,3);
  Matrix P_true(3,3);
  Matrix L_true(3,3);
  Matrix U_true(3,3);
  double elementsA [9] = {1,-2.,-2.,-1.,3,0,-1,1,-5};
  double elementsU [9] = {1,-2,-2,0,1,-2,0,0,-9};
  double elementsP [9]  = {1,0,0,0,1,0,0,0,1};
  double elementsL [9]  = {1,0,0,-1,1,0,-1,-1,1};
  int count = 0;
  for (int i = 1; i <= 3; i++)
  {
    for (int j = 1; j <= 3; j++)
    {
      A(i,j) = elementsA[count];
      U_true(i,j) = elementsU[count];
      P_true(i,j) = elementsP[count];
      L_true(i,j) = elementsL[count];
      count += 1;
    }
  }
  int sign;
  Matrix P(3,3);
  Matrix U(3,3);
  Matrix L(3,3);
  std::tie(P, L, U,sign) = lu(A);
  
  if (P==P_true && L==L_true && U==U_true)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }

}

void Det()
{
  // Test det()
  Matrix A(3,3);
  int count = 0;
  double elements [9] = {1,4,3,2,1,5,3,2,1};
  for (int i = 1; i <= 3; i++)
  {
    for (int j = 1; j <= 3; j++)
    {
      A(i,j) = elements[count];
      count += 1;
    }
  }

  double det_true = 46;
  double det = A.det();

  if (det_true == det)
  {
    std::cout << "Test - " << __func__ << ", Result: Passed." << std::endl;
  }
  else
  {
    std::cout << "Test - " << __func__ << ", Result: Failed." << std::endl;
  }
}


// void QR()
// {
//   std::cout << "QR" <<std::endl;
//   Matrix A(4,3);
//   int count = 0;
//   double elements [12] = {1,4,3,
//                          4,1,5,
//                          3,5,1,
//                          5,6,7};
//   for (int i = 1; i <= 4; i++)
//   {
//     for (int j = 1; j <= 3; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }
  

//   Matrix Q(4,3);
//   Matrix R(3,3);
//   std::tie(Q, R) = qr(A, true);
//   std::cout << Q << std::endl;
//   std::cout << R << std::endl;
//  // double detR = 1;
//   //std::cout << Q.T() * Q << std::endl;
//   //std::cout << det(R) << " = " << detR << std::endl;
//   std::cout << Q*R << std::endl;

// }

// void LSQ()
// {
//   std::cout << "LSQ" <<std::endl;
//   Matrix A(4,3);
//   int count = 0;
//   double elements [12] = {1,4,3,
//                          4,1,5,
//                          3,5,1,
//                          5,6,7};
//   for (int i = 1; i <= 4; i++)
//   {
//     for (int j = 1; j <= 3; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }

//   Matrix b(4,1);
//   b(1,1) = 1; b(2,1) = 2; b(3,1) = 3; b(4,1) = 4;
  
//   Matrix x = lsq(A, b);

// }


// void HessReduction()
// {
//   Matrix A(4,4);
//   int count = 0;
//   double elements [16] = {1,4,3,5,
//                          4,1,5,6,
//                          3,5,1,7,
//                          5,6,7,1};
//   for (int i = 1; i <= 4; i++)
//   {
//     for (int j = 1; j <= 4; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }

//   std::cout <<  hessenbergReduction(A) << std::endl;

// }

// void EigenVal()
// { 
//   {
//     int n = 4;
//   Matrix A(n,n);
//   int count = 0;
//   double elements [n*n] = {16,2,3,13,
//                            5,11,10,8,
//                            9,7,6,12,
//                            4,14,15,1};;
//   for (int i = 1; i <= n; i++)
//   {
//     for (int j = 1; j <= n; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }

//   Matrix vals(n,1);
//   vals = eigenVal(A);
//   print(vals);
//   }
// }

// void NormAll()
// {
//   int n = 4;
//   Matrix A(n,n);
//   int count = 0;
//   double elements [n*n] = {16,2,3,13,
//                            5,11,10,8,
//                            9,7,6,12,
//                            4,14,15,1};
//   for (int i = 1; i <= n; i++)
//   {
//     for (int j = 1; j <= n; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }
//   try
//   {
//     {
//       norm(A);
//     }
//   }
//   catch(Exception &e)
//   {
//     e.DebugPrint();
//   }
  
//   std::cout << norm(A) << std::endl;
// }

// void GaussianElimination()
// {
//   int n = 3;
//   Matrix A(n,n);
//   int count = 0;
//   double elements [n*n] = {2,1,-1,
//                            -3,-1,2,
//                            -2,1,2};
//   for (int i = 1; i <= n; i++)
//   {
//     for (int j = 1; j <= n; j++)
//     {
//       A(i,j) = elements[count];
//       count += 1;
//     }
//   }
// }


}




