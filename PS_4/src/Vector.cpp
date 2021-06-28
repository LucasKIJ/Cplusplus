#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"


// constructor that creates vector of given size with
// double precision entries all initially set to zero
Vector::Vector(int sizeVal) : Matrix(sizeVal,1)
{
}





// copy constructor - creates vector with the same entries as v1
Vector::Vector(const Vector& v1) : Matrix(v1)
{
}






// destructor - deletes pointer
Vector::~Vector()
{
  delete[] mData;
}


////// ************ GETTERS AND SETTERS ************ ////// 
// Set value using 0-indexing
void Vector::setValue(double value, int row) const
{
  int col = 0;
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Vector::setValue - Out of range", "The element coordinates given when\
                      setting an element value are out of range (use 0-indexing).");
  }
  mData[mCol * row + col] = value;
}

// Get value using 0-indexing
double& Vector::getValue(int row) const
{
  int col = 0;
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Vector::getValue - Out of range", "The element coordinates given when\
                      getting an element value are out of range (use 0-indexing).");
  }
  return mData[mCol * row + col];
}

Vector& Vector::operator=(const Vector& v)
{
//  check both Matrixs have same dimension
  
  if (mRow != v.mRow && mCol != v.mCol)
  {
    throw Exception("Dimension mismatch", "Matrix assignment operator - Matrices\
                        have different dimensions");
  }
  else
  {
    for (int i=0; i<v.mSize; i++)
	  {
	    mData[i] = v.mData[i];
	  }
  }
  return *this;
}

Vector& Vector::operator=(const Matrix& m)
{
//  check both Matrixs have same dimension
  
  if (mRow != m.getNumRow() && mCol != m.getNumCol())
  {
    throw Exception("Dimension mismatch", "Matrix assignment operator - Matrices\
                        have different dimensions");
  }
  if (m.getNumRow() != 1 && m.getNumCol() != 1)
  {
    throw Exception("Dimension mismatch", "Input is not a vector");
  }

  else
  {
    if (m.getNumRow() == 1)
    {
      for (int i=0; i<  m.getNumCol(); i++)
      {
        mData[i] = m.getValue(0,i);
      }
    }
    else
    {
      for (int i=0; i< m.getNumRow(); i++)
      {
        mData[i] = m.getValue(i,0);
      }
    }
  }
  return *this;
}



double& Vector::operator()(int i)
{

  if (i < 1)
    {
      throw Exception("out of range",
		  "accessing vector through () - index too small");
    }
  else if (i > mSize)
    {
      std::cout << mSize;
      throw Exception("length mismatch",
		  "accessing vector through () - index too high");
    }


  return mData[i-1];
}

// Friend function
// calculate p-norm of a vector v
// default value for p is 2
double norm(Vector& v, int p)
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=1; i<=length(v); i++)
    {
      temp = fabs(v(i));
      norm_val += pow(temp, p);
    }

  return pow(norm_val, 1.0/((double) (p)));
}
// Member method
// calculate p-norm of a vector v
// default value for p is 2
double Vector::norm(int p) const
{
  double temp, norm_val;

  norm_val = 0.0;
  for (int i=0; i<mSize; i++)
    {
      temp = fabs(mData[i]);
      norm_val += pow(temp, p);
    }

  return pow(norm_val, 1.0/((double) (p)));
}




// return length of a vector
int length(const Vector& v)
{
  return v.mSize;
}