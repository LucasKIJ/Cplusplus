#include "Matrix.hpp"
//#include "Vector.hpp"


// constructor that creates Matrix of given size with
// double precision entries all initially set to zero
Matrix::Matrix(int row, int col) : mRow(row), mCol(col), mSize(row * col)
{
  mData=new double[mSize];
  for (int i=0; i< mSize; i++)
  {
     mData[i] = 0.0;
  }
}


Matrix::Matrix(int n) : mRow(n), mCol(n), mSize(n * n)
{
  mData=new double[mSize];
  for (int i=0; i< mSize; i++)
  {
     mData[i] = 0.0;
  }
}





// copy constructor - creates Matrix with the same entries as v
Matrix::Matrix(const Matrix& m) : mRow(m.mRow), mCol(m.mCol), mSize(m.mSize)
{
  mData=new double[mSize];
  for (int i=0; i< mSize; i++)
  {
    mData[i] = m.mData[i];
  }
}




// destructor - deletes pointer
Matrix::~Matrix()
{
  delete[] mData;
}


////// ************ GETTERS AND SETTERS ************ ////// 
// Set value using 0-indexing
void Matrix::setValue(double value, int row, int col)
{
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Matrix::setValue - Out of range", "The element coordinates given when\
                      setting an element value are out of range (use 0-indexing).");
  }
  mData[mCol * row + col] = value;
}

// Get value using 0-indexing
double& Matrix::getValue(int row, int col) const
{
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Matrix::getValue - Out of range", "The element coordinates given when\
                      getting an element value are out of range (use 0-indexing).");
  }
  return mData[mCol * row + col];
}
///////////////////////////////////////////////////////////

//////// ************ Identity matrix ************ //////// 
// Create identity matrix
const Matrix eye(int row, int col)
{
  Matrix identity(row,col);
  int n;
  if (row < col)
  {
    n = row;
  }
  else{
    n = col;
  }

  for (int i = 0; i < n; i++)
  {
    identity.setValue(1,i,i);
  }
  return identity;
}

// Square identity matrix
const Matrix eye(int n)
{
  return eye(n,n);
}
//////////////////////////////////////////////////////////




//////// ******** Arithmatic operators ******** //////////
// definition of + between two Matrixs
Matrix operator+(const Matrix& m1, 
                        const Matrix& m2)
{
  //  add the Matrixs
  //  if one Matrix is shorter than the other assume missing entries are 0
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    Matrix ma(m1.mRow, m1.mCol);
    for (int i=0; i<m1.mSize; i++)
	  {
	    ma.mData[i] = m1.mData[i] + m2.mData[i];
	  }
    return ma;
  }
  else
  {
    throw Exception("Different dimensions", "Matrix add - Matrices are different dimensions\n");
  }
}



// definition of - between two Matrices
Matrix operator-(const Matrix& m1, 
                        const Matrix& m2)
{
  
  // Check that number of rows and columns are equal
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    Matrix ms(m1.mRow, m1.mCol);
    for (int i=0; i<m1.mSize; i++)
	  {
	    ms.mData[i] = m1.mData[i] - m2.mData[i];
	  }
    return ms;
  }
  else
  {
    throw Exception("Different dimensions", "Matrix subtract - Matrices are different dimensions\n");
  }
}



// definition of matrix procture
Matrix operator*(const Matrix& m1, const Matrix& m2)
{
//  check Matrice are of the same dimension,
//  i.e m1 number of col = m2 number of row
  if (m1.mCol == m2.mRow)
    {
      
      int n = m1.mCol;
      double prod = 0;
      Matrix mp(m1.mRow, m2.mCol);
      for (int i = 0; i < m1.mRow; i++)
      {
        for (int j = 0; j < m2.mCol; j++)
        {
          prod = 0;
          for (int k = 0; k < n; k++)
          {
            prod += m1.getValue(i,k) * m2.getValue(k,j);
          }
          mp.setValue(prod,i,j);
        }
      }
      return mp;
    }
  else
    {
      throw Exception("Different dimensions", "Matrix product - Number of columns in the first matrix\
                        does not equal number of rows in the second matrix\n");
    }
}


// definition of multiplication between a Matrix and a scalar
Matrix operator*(const Matrix& m, const double& a)
{

//  create a Matrix of the same dimension as m with entries equal to a*m
  Matrix mp(m.mRow, m.mCol);

  for (int i=0; i<m.mSize; i++)
    {
      mp.mData[i] = a * m.mData[i];
    }

  return mp;
}


// definition of multiplication between a scalar and a Matrix
Matrix operator*(const double& a, const Matrix& m)
{

//  create a Matrix of the same dimension as m with entries equal to a*m
  Matrix mp(m.mRow, m.mCol);

  for (int i=0; i<m.mSize; i++)
    {
      mp.mData[i] = a * m.mData[i];
    }
  return mp;
}



// definition of division of a Matrix by a scalar
Matrix operator/(const Matrix& m, const double& a)
{
  if (a == 0.0)
  {
     throw Exception("Division by 0", "Attempt to divide by zero in the\
                       Matrix-scalar division operator");
  }
//  create a Matrix of the same length as v with entries equal to v/a

  Matrix md(m.mRow, m.mCol);

  for (int i=0; i< m.mSize; i++)
    {
      md.mData[i] = m.mData[i] / a;
    }

  return md;
}
///////////////////////////////////////////////////////////////////////////////

/////////////// *********** Unary operator *************** ////////////////////
// definition of the unary operator -
Matrix operator-(const Matrix& m)
{
//  create a Matrix w with entries equal to -v
  Matrix ms(m.mRow, m.mCol);

  for (int i=0; i< m.mSize; i++)
    {
      ms.mData[i] = -m.mData[i];
    }

  return ms;
}



// definition of Matrix operator ()
// Accessing this way is done with 1-indexing

double& Matrix::operator() (int i, int j)
{
  if (i <= 0 || i > mRow || j <= 0 || j > mCol)
  {
    throw Exception("(i,j) - Out of range", "The element coordinates given when\
                      getting an element value are out of range (use 1-indexing).");
  }
  return mData[mCol * (i-1) + (j-1)];
}


//////////// ******** Print functions ************ /////////////////
void print(const Matrix& m)
{
  std::cout << std::fixed << std::setprecision(1);
  for (int i = 0; i < m.mRow; i++) 
  {
    std::cout << "[";
    for (int j = 0; j < m.mCol; j++) 
    {
      std::cout << (j > 0 ? " " : "") << std::setw(4);
      std::cout << m.getValue(i,j);
    }
  std::cout << " ]" << std::endl;
    }
}

std::ostream& operator<<(std::ostream& output, const Matrix& m) 
{
  output << std::fixed << std::setprecision(1);
  for (int row = 0; row < m.mRow; row++) 
  {
    output << "[";
    for (int col = 0; col < m.mCol; col++) 
    {
      output << (col > 0 ? " " : "") << std::setw(4);
      output << m.mData[m.mCol * row + col];
    }
    output << " ]" << std::endl;
  }
  return output;
}
///////////////////////////////////////////////////////////////////



///////////// ****** Assignment operator ******** /////////////////
// definition of Matrix operator =
Matrix& Matrix::operator=(const Matrix& m)
{
//  check both Matrixs have same dimension
  
  if (mRow != m.mRow && mCol != m.mCol)
    {
      throw Exception("Dimension mismatch", "Matrix assignment operator - Matrices\
                          have different dimensions");
    }
  else
  {
    for (int i=0; i<m.mSize; i++)
	  {
	    mData[i] = m.mData[i];
	  }
  }
  return *this;
}
///////////////////////////////////////////////////////////////////

//////////// ******* Transpose operator ******* ///////////////////
Matrix Matrix::T() const
{
    Matrix mt(mCol, mRow);

    for(int i=0; i < mRow; ++i)
    {
        for(int j=0; j < mCol; ++j)
        {
            mt.setValue(getValue(i,j), j,i);
        }
    }

    return mt;
}

Matrix transpose(const Matrix& m)
{
  return m.T();
}

///////////////////////////////////////////////////////////////////



///////// ********* Matrix norms *********** ////////////////
// Friend function
// calculate p-norm of a Matrix v
// default value for p is 2
double norm(Matrix& m, int p)
{
  return m.norm(p);
}

double norm(Matrix& m, std::string p)
{
  return m.norm(p);
}
// Member method
// calculate p-norm of a Matrix v
double Matrix::norm(int p) const
{
  if (p==2)
  {
    std::cout << "Unfinished" << std::endl;
  }
  else if (p==1)
  {
    double max_sum = 0;
    double sum;

    for (int j=0; j < mCol; j++)
    {
      sum = 0;
      for (int i=0; i < mRow; i++)
      {
        sum += fabs(getValue(i,j));
      }

      if (sum > max_sum)
      {
        max_sum = sum;
      }
    }

    return max_sum;
  }
  else
  {
  throw Exception("Matrix::norm - norm selected not available","Choice of norm not available,\
                       please choose from the 1-norm (p=1), 2-norm (p=2), infinity-norm (p=\"inf\")\
                       or the Frobenius-norm (p=\"frob\")");
  }
}

double Matrix::norm(std::string p) const
{
  if (p=="frob")
  {
    double sum = 0;
    for (int i=0; i < mRow; i++)
    {
      for (int j=0; j < mCol; j++)
      {
        sum += pow(getValue(i,j),2.);
      }
    }

    return pow(sum,1./2.);
  }
  else if (p=="inf")
  {
    double max_sum = 0;
    double sum;

    for (int i=0; i < mRow; i++)
    {
      sum = 0;
      for (int j=0; j < mCol; j++)
      {
        sum += fabs(getValue(i,j));
      }

      if (sum > max_sum)
      {
        max_sum = sum;
      }
    }
    return max_sum;
  }
  else
  {
  throw Exception("Matrix::norm - norm selected not available","Choice of norm not available,\
                       please choose from the 1-norm (p=1), 2-norm (p=2), infinity-norm (p=\"inf\")\
                       or the Frobenius-norm (p=\"frobenius\")");
  }
}

/////////////////////////////////////////////////////////////////////////////////////



// return size of a Matrix
int* size(const Matrix& m)
{
  int* dims;
  dims = new int [2];
  dims[0] = m.mRow;
  dims[1] = m.mCol;
  return dims;
}

