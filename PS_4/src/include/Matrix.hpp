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
#include <tuple>
#include "Exception.hpp" //  This class throws errors using the class "error"

class Matrix {
protected:
  // member variables
  double *mData;   // data stored in Matrix
  const int mRow;  // number of rows of Matrix
  const int mCol;  // number of columns
  const int mSize; // number of elements

public:
  /////// ---------- Constructors and Destructors -----------//////
  // constructors
  // No default constructor
  // overridden copy constructor
  Matrix(const Matrix &m);
  // construct Matrix of given size
  Matrix(int row, int col);
  // construst square matrix
  Matrix(int n);

  // destructor
  ~Matrix();

  //////////---------- Setters and getters ------------- /////////////
  // Set value
  void setValue(double value, int row, int col) const;
  // Get value
  double &getValue(int row, int col) const;

  // Get number of rows
  int getNumRow() const;
  int getNumCol() const;

  ///////////// ------------- Operators ------------- ////////////////
  // All "friend" external operators and functions are declared as friend inside
  // the class (here)
  // but their actual prototype definitions occur outside the class.
  // Binary operators
  friend Matrix operator+(const Matrix &m1, const Matrix &m2);
  friend Matrix operator+(const Matrix &m, const double &a);
  friend Matrix operator+(const double &a, const Matrix &m);
  friend Matrix operator-(const Matrix &m1, const Matrix &m2);
  friend Matrix operator-(const Matrix &m, const double &a);
  friend Matrix operator-(const double &a, const Matrix &m);
  friend Matrix operator*(const Matrix &m1, const Matrix &m2);
  friend Matrix operator*(const Matrix &m, const double &a);
  friend Matrix operator*(const double &a, const Matrix &m);
  friend Matrix operator/(const Matrix &m, const double &a);
  // assignment
  Matrix &operator=(const Matrix &m);
  // equality
  friend bool operator==(const Matrix &m1, const Matrix &m2);
  // inequality
  friend bool operator!=(const Matrix &m1, const Matrix &m2);
  // indexing
  double &operator()(int i, int j);
  // Unary operator
  friend Matrix operator-(const Matrix &m);

  ////////// ----------- Linear Algebra helps ------------ ///////
  friend double dot(const Matrix &m1, const Matrix &m2);
  Matrix T() const;
  friend Matrix transpose(const Matrix &m);
  // Max
  double rowMax(int row, int col = 0, bool abs_true = true) const;
  double colMax(int col, int row = 0, bool abs_true = true) const;
  int argRowMax(int row, int col = 0, bool abs_true = true) const;
  int argColMax(int col, int row = 0, bool abs_true = true) const;
  // Row and Column swapping
  void swapRow(int row1, int row2) const;
  void swapCol(int col1, int col2) const;
  // swapping
  void swap(int row1, int col1, int row2, int col2) const;
  // eye
  friend const Matrix eye(int row, int col);
  friend const Matrix eye(int n);
  // random matrix
  friend const Matrix rand(int row, int col);
  friend const Matrix rand(int n);
  // output
  friend void print(const Matrix &m);
  friend std::ostream &operator<<(std::ostream &output, const Matrix &m);
  // norm (as a member method)
  double norm(int p) const;
  double norm(std::string p) const;
  // functions that are friends
  friend double norm(Matrix &m, int p);
  friend double norm(Matrix &m, std::string p);
  friend int *size(const Matrix &m);

  //////// --------- Linear Algebra Algorithms ---------- //////////

  friend Matrix forwardSub(const Matrix &L, const Matrix &c);
  friend Matrix backwardSub(const Matrix &U, const Matrix &b);
  // Gaussian elimination solver
  friend Matrix gaussianElimination(const Matrix &A, const Matrix &b);
  friend Matrix solverLU(const Matrix &A, const Matrix &b);
  // LU Decomposition
  friend std::tuple<Matrix, Matrix, Matrix, int> lu(const Matrix &A);
  // Determinant
  double det() const;
  friend double det(const Matrix &A);
  // QR factorisation
  friend std::tuple<Matrix, Matrix> qr(const Matrix &A, bool reduced = false);
  // Least-squares
  friend Matrix lsq(const Matrix &A, const Matrix &b);
  // QR Algorithm
  friend Matrix hessenbergReduction(const Matrix &A);
  friend Matrix eigenVal(const Matrix &A, double tol = 1.0e-10);
};

// All "friend" external operators and functions are declared as friend inside
// the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Matrix operator+(const Matrix &m1, const Matrix &m2);
Matrix operator+(const Matrix &m, const double &a);
Matrix operator+(const double &a, const Matrix &m);
Matrix operator-(const Matrix &m1, const Matrix &m2);
Matrix operator-(const Matrix &m, const double &a);
Matrix operator-(const double &a, const Matrix &m);
Matrix operator*(const Matrix &m1, const Matrix &m2);
Matrix operator*(const Matrix &m, const double &a);
Matrix operator*(const double &a, const Matrix &m);
Matrix operator/(const Matrix &m, const double &a);

double dot(const Matrix &m1, const Matrix &m2);

// Solver
Matrix gaussianElimination(const Matrix &A, const Matrix &b);
Matrix solverLU(const Matrix &A, const Matrix &b);

// LU Decomposition
std::tuple<Matrix, Matrix, Matrix, int> lu(const Matrix &A);
Matrix forwardSub(const Matrix &L, const Matrix &c);
Matrix backwardSub(const Matrix &U, const Matrix &b);
// QR
std::tuple<Matrix, Matrix> qr(const Matrix &A, bool reduced);
// Least-squares
Matrix lsq(const Matrix &A, const Matrix &b);
// QR Algorithm
Matrix hessenbergReduction(const Matrix &A);
Matrix eigenVal(const Matrix &A, double tol);

// Det
double det(const Matrix &A);

// Unary operator
Matrix operator-(const Matrix &m);

// equality
bool operator==(const Matrix &m1, const Matrix &m2);

// print
void print(const Matrix &m);

// transpose
Matrix transpose(const Matrix &m);

// function prototypes
double norm(Matrix &m, int p = 2);
// Prototype signature of length() friend function
int length(const Matrix &m);

// eye
const Matrix eye(int row, int col);
const Matrix eye(int n);
// random matrix
const Matrix rand(int row, int col);
const Matrix rand(int n);

#endif
