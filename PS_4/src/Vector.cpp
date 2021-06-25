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