#ifndef MATRIXDEF
#define MATRIXDEF
//#include "Vector.hpp"



//  **********************
//  *  Class of matrix  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Exception.hpp"//  This class throws errors using the class "error"

class Matrix
{
protected:
   // member variables
   double* mData;   // data stored in Matrix
   const int mRow;      // number of rows of Matrix
   const int mCol;      // number of columns
   const int mSize;     // number of elements

public:
   // constructors
   // No default constructor
   // overridden copy constructor
   Matrix(const Matrix& m);
   // construct Matrix of given size
   Matrix(int row, int col);
   // construst square matrix
   Matrix(int n);

   // destructor
   ~Matrix();
   
   // Set value
   void setValue(double value, int row, int col); 
   // Get value
   double& getValue(int row, int col) const;
   


   // All "friend" external operators and functions are declared as friend inside the class (here)
   // but their actual prototype definitions occur outside the class.
   // Binary operators
   friend Matrix operator+(const Matrix& m1, const Matrix& m2);
   friend Matrix operator-(const Matrix& m1, const Matrix& m2);
   friend Matrix operator*(const Matrix& m1, const Matrix& m2);
   friend Matrix operator*(const Matrix& m, const double& a);
   friend Matrix operator*(const double& a, const Matrix& m);
   friend Matrix operator/(const Matrix& m, const double& a);
   // Unary operator
   friend Matrix operator-(const Matrix& m);

   //other operators

   //assignment
   Matrix& operator=(const Matrix& m);
   //indexing
   double& operator() (int i, int j);
   //transposing
   Matrix T() const;
   friend Matrix transpose(const Matrix& m);


   //output
   friend void print(const Matrix& m);
   friend std::ostream& operator<<(std::ostream& output, const Matrix& m);

   //norm (as a member method)
   double norm(int p) const;
   double norm(std::string p) const;
   // functions that are friends
   friend double norm(Matrix& m, int p);
   friend double norm(Matrix& m, std::string p);
   friend int* size(const Matrix& m);



};


// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Matrix operator+(const Matrix& m1, const Matrix& m2);
Matrix operator-(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator*(const Matrix& m, const double& a);
Matrix operator*(const double& a, const Matrix& m);
Matrix operator/(const Matrix& m, const double& a);
Matrix operator/(const Matrix& m1, const Matrix& m2);
// Unary operator
Matrix operator-(const Matrix& m);

//print
void print(const Matrix& m);

//transpose
Matrix transpose(const Matrix& m);

// function prototypes
double norm(Matrix& m, int p=2);
// Prototype signature of length() friend function
int length(const Matrix& m);

const Matrix eye(int row, int col);
const Matrix eye(int n);

#endif
