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
void Matrix::setValue(double value, int row, int col) const
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


////////////// ******** Random Matrix ********* //////////
const Matrix rand(int row, int col)
{
  Matrix R(row,col);
  for (int i=0; i < R.mSize; i++)
  {
     R.mData[i] = rand();
  }
  return R;
}

const Matrix rand(int n)
{
  return rand(n,n);
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
  std::cout << std::endl;
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


/////////////////// ********** Matrix Size ************ //////////////////////////
// return size of a Matrix
int* size(const Matrix& m)
{
  int* dims;
  dims = new int [2];
  dims[0] = m.mRow;
  dims[1] = m.mCol;
  return dims;
}
//////////////////////////////////////////////////////////////////////////////////

////////////////// ********** Max functions *********** //////////////////////////
double Matrix::rowMax(int row, int col, bool abs_true) const
{
  double max_val ;
  double temp;
  max_val = getValue(row,col);
  for (int j = col; j < mCol; j++)
  {
    if (abs_true)
    {
      temp = fabs(getValue(row,j));
    }
    else
    {
      temp = getValue(row,j);
    }
    if (max_val < temp)
    {
      max_val = temp;
    }
  }
  return max_val;
}


double Matrix::colMax(int col, int row, bool abs_true) const
{
  double max_val;
  double temp;

  max_val = getValue(row,col);
  for (int i = row; i < mRow; i++)
  {
    if (abs_true)
    {
      temp = fabs(getValue(i,col));
    }
    else
    {
      temp = getValue(i,col);
    }
    if (max_val < temp)
    {
      max_val = temp;
    }
  }
  return max_val;
}

int Matrix::argRowMax(int row, int col, bool abs_true) const
{
  int arg_max = col;
  double max_val ;
  double temp;

  max_val = getValue(row,col);
  for (int j = col; j < mCol; j++)
  {
    if (abs_true)
    {
      temp = fabs(getValue(row,j));
    }
    else
    {
      temp = getValue(row,j);
    }
    if (max_val < temp)
    {
      max_val = temp;
      arg_max = j;
    }
  }
  return arg_max;
}

int Matrix::argColMax(int col, int row, bool abs_true) const
{
  int arg_max = row;
  double max_val;
  double temp;
  max_val = getValue(row,col);
  for (int i = row; i < mRow; i++)
  {
    if (abs_true)
    {
      temp = fabs(getValue(i,col));
    }
    else
    {
      temp = getValue(i,col);
    }
    if (max_val < temp)
    {
      max_val = temp;
      arg_max = i;
    }
  }
  return arg_max;
}
//////////////////////////////////////////////////////////////////////////////////


/////////////////// ********** Boolean Operators *********** /////////////////////
bool operator==(const Matrix& m1, const Matrix& m2)
{
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    for (int i = 0; i < m1.mSize; i++)
    {
      if (m1.mData[i] != m2.mData[i])
      {
        return false;
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

bool operator!=(const Matrix& m1, const Matrix& m2)
{
  return !(m1 == m2);
}
//////////////////////////////////////////////////////////////////////////////////

////////////////// *************** Swap elements ************* ////////////////////
void Matrix::swap(int row1, int col1, int row2, int col2) const
{
  double temp;
  temp = getValue(row1,col1);
  setValue(getValue(row2,col2),row1,col1);
  setValue(temp,row2,col2);
}


void Matrix::swapRow(int row1, int row2) const
{
  double temp;
  for (int j = 0; j < mCol; j++){
    temp = getValue(row1,j);
    setValue(getValue(row2,j),row1,j);
    setValue(temp,row2,j);
  }
}

void Matrix::swapCol(int col1, int col2) const
{
  double temp;
  for (int i = 0; i < mRow; i++){
    temp = getValue(i,col1);
    setValue(getValue(i,col2),i,col1);
    setValue(temp,i,col2);
  }
}

//////////////////////////////////////////////////////////////////////////////////

///////////// =#=#=#=#=#=#=#=#=#  Linear Solvers #=#=#=##=#=#=#=# ////////////////
//////////////////////////////////////////////////////////////////////////////////

////////////// *********** LU Decomposition ************* ///////////////////////

std::tuple <Matrix, Matrix, Matrix> LU(const Matrix& A)
{
  int m = A.mCol;
  int n = A.mRow;
  Matrix U(A);
  Matrix L = eye(n);
  Matrix P = eye(n);

  int i;
  double temp;

  for (int k = 0; k < m - 1; k++)
  {
    i = U.argColMax(k,k);

    for (int j = k; j < m; j++)
    {
      U.swap(k,j,i,j);
    }

    for (int j = 0; j < k - 1; j++)
    {
      L.swap(k,j,  i,j);
    }

    P.swapRow(k,i);

    for (int i = k + 1; i < m; i++)
    {
      temp = U.getValue(i,k) / U.getValue(k,k);
      L.setValue(temp, i,k);
      for (int j = k; j < m; j++)
      {
        temp = U.getValue(i,j) - L.getValue(i,k) * U.getValue(k,j);
        U.setValue(temp, i,j);
      }
    }
  }

  return {P, L, U};
}



/////////////////////////////////////////////////////////////////////////////////



///////////// ********* Gaussian Elimination ************ ////////////////////////
Matrix gaussianElimination(const Matrix A, const Matrix b)
{
  double M;
  // Perform elimination
  for (int k = 0; k < A.mRow; k++)
  {      
      // Get pivot position
      int i_max = A.argColMax(k,k);
      // No pivot in this column.
      if (A.getValue(i_max,k) == 0)
      {
        continue;
      }

      // If pivot exists and not equal to current column
      if (i_max != k)
      {
        A.swapRow(k, i_max);  
        // Swap vector entries
        b.swapRow(k, i_max);
      }



      for (int i = k+1; i < A.mRow; i++)
      {
          M = A.getValue(i,k) / A.getValue(k,k);
          for (int j = k; j < A.mRow; j++)
          {    
            A.setValue(A.getValue(i,j) - M * A.getValue(k,j), i, j);
          }
          for (int j = 0; b.mRow; j++)
          {
            b.setValue(b.getValue(i,j)- M * b.getValue(k,j), i, j);
          }
      }
    }

  // Declare x
  Matrix x(A.mRow,b.mCol);
  double sum;
  for (int j = 1; j < x.mCol; j++)
  {
    for (int k = A.mCol - 1; k >= 0; k--)
    {
      sum = 0;
      for (int i = k; i < A.mCol; i++)
      {
        sum += A.getValue(k,j) * x.getValue(i,j);
      }
      x.setValue((1. / A.getValue(k,k)) * (b.getValue(k,j) - sum), k, j);
    }
  }
  return x;
}


