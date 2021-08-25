#ifndef VECTORDEF
#define VECTORDEF

//  **********************
//  *  Class of vectors  *
//  **********************

//  Class written in such a way that code similar to Matlab
//  code may be written

#include <cmath>
#include "Exception.hpp" //  This class throws errors using the class "error"
#include "Matrix.hpp"

class Vector : public Matrix {
public:
  // constructors
  // No default constructor
  // overridden copy constructor
  Vector(const Vector &v1);
  // construct vector of given size
  Vector(int sizeVal);

  // destructor
  ~Vector();

  // Set value
  void setValue(double value, int row) const;
  // Get value
  double &getValue(int row) const;

  // assignment
  Vector &operator=(const Vector &m);
  Vector &operator=(const Matrix &m);

  // indexing
  double &operator()(int i);

  // norm (as a member method)
  double norm(int p = 2) const;
  // functions that are friends
  friend double norm(Vector &v, int p);
  friend int length(const Vector &v);
};

// All "friend" external operators and functions are declared as friend inside
// the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators

// function prototypes
double norm(Vector &v, int p = 2);
// Prototype signature of length() friend function
int length(const Vector &v);

#endif
