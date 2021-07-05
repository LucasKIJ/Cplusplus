#include "Matrix.hpp"
#include "Vector.hpp"

// constructor that creates Matrix of given size with
// double precision entries all initially set to zero
Matrix::Matrix(int row, int col) : mRow(row), mCol(col), mSize(row * col)
{
/** @brief Constructor for the matrix class
 * @details 
 *    Creates a matrix object with the desired number
 *    of rows and columns with element values initialised
 *    as zero. Inputs row and col become attributes mRow and mCol.
 *    mSize := mRow * mSize.
 * @param row The desired number of row in the matrix
 * @param col The desired number of row in the matrix
 * @return    Matrix object with desired number of row and columns
 * 
 */

  // If invalid dimensions entered throw warning
  if (mCol < 1 || mRow < 1)
  {
    throw Exception("Matrix::Matrix - Invalid dimensions", "Please ender dimension greater or equal"
                      "to one.");
  }
  // If valid, create new pointer to an array of size mRow * mCol
  mData=new double[mSize];
  // Initialise values as zero
  for (int i=0; i< mSize; i++)
  {
     mData[i] = 0.0;
  }
}


Matrix::Matrix(int n) : mRow(n), mCol(n), mSize(n * n)
{
/** @brief Constructor for the matrix class
 * @details 
 *    Creates a square matrix object with the desired number
 *    of rows and columns with element values initialised
 *    as zero. Input n become attributes mRow and mCol.
 *    mSize := mRow * mSize.
 * @param n The desired number of rows and columns in the square matrix
 * @return    Matrix object with desired number of row and columns
 */

  // If invalid dimensions entered throw warning
  if (n < 1)
  {
    throw Exception("Matrix::Matrix - Invalid dimensions", "Please ender dimension"
                      " greater or equal to one.");
  }

  // If valid, create new pointer to an array of size mRow * mCol
  mData=new double[mSize];
  // Initialise values as zeros
  for (int i=0; i< mSize; i++)
  {
     mData[i] = 0.0;
  }
}


// // copy constructor - creates Matrix with the same entries as v
Matrix::Matrix(const Matrix& m) : mRow(m.mRow), mCol(m.mCol), mSize(m.mSize)
{
  /** @brief Copy constructor for the matrix class
 * @details 
 *    Creates a copy of matrix object with the same dimensions and 
 *    element values.
 * @param m The matrix object from which a copy is made
 * @return    Matrix object with desired number of row and columns
 * 
 */

  // Create new pointer to an array of mSize = mRow * mCol
  mData=new double[mSize];
  // Initialise data values with data values from input Matrix object
  for (int i=0; i< mSize; i++)
  {
    mData[i] = m.mData[i];
  }
}




// destructor - deletes pointer
Matrix::~Matrix()
{
/** @brief Destructor for the matrix class
 * @details 
 *    Once a matrix object goes out of scope the destructor
 *    is called. It deletes the memory array used to store
 *    the matrix element data.
 */
  delete[] mData;
}


////// ************ GETTERS AND SETTERS ************ ////// 
// Set value using 0-indexing
void Matrix::setValue(double value, int row, int col) const
{
  /** @brief Setter: Used to set values of matrix elements.
 * @details 
 *    Given a value and valid coordinates this method sets the
 *    value at the coordinate location to the given value.
 * @param value The value assign to the specified coordinates
 * @param row The row position of the element to be updated
 * @param col The row position of the element to be updated
 */

  // Throw error if element position entered is out of range
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Matrix::setValue - Out of range", "The element coordinates given when"
                      "setting an element value are out of range (use 0-indexing).");
  }
  // mData saved in 1-D array, element (i,j) = mData[mCol * i + j]
  // Change mData value at position entered.
  mData[mCol * row + col] = value;
}

// Get value using 0-indexing
double& Matrix::getValue(int row, int col) const
{
  /** @brief Getter: Used to get values of matrix elements.
   * @details 
   *    Given valid coordinates this method gets the
   *    value at the coordinate location (using 0-indexing)
   * @param row The row position of the element to be returned
   * @param col The row position of the element to be returned
   * @returns Value of element at specified coordinate location
   */

  // Throw error if element position entered is out of range 
  if (row < 0 || row >= mRow || col < 0 || col >= mCol)
  {
    throw Exception("Matrix::getValue - Out of range", "The element coordinates given when"
                      "getting an element value are out of range (use 0-indexing).");
  }
  // mData saved in 1-D array, element (i,j) = mData[mCol * i + j]
  // return mData value at position entered.
  return mData[mCol * row + col];
}

int Matrix::getNumRow() const
{
  /** @brief Getter: Gets the number of rows in Matrix object.
   * @return The number of rows in matrix object
 */
  return mRow;
}

int Matrix::getNumCol() const
{
  /** @brief Getter: Gets the number of columns in Matrix object.
   * @return The number of columns in matrix object
 */
  return mCol;
}


///////////////////////////////////////////////////////////

//////// ************ Identity matrix ************ //////// 
// Create identity matrix
const Matrix eye(int row, int col)
{
  /** @brief Identity matrix constructor
 * @details 
 *    Creates an identity matrix with the desired number
 *    of rows and columns with ones on the diagonal and zeros
 *    everywhere else
 * @param row The desired number of row in the matrix
 * @param col The desired number of row in the matrix
 * @return    Matrix object with desired number of row and columns
 *            with ones on the diagonal
 * 
 */
  // Initialise Matrix object with desired dimensions
  Matrix identity(row,col);
  
  // Determine lenth of diagonal if non square
  // matrix desired.
  int n;
  if (row < col)
  {
    n = row;
  }
  else{
    n = col;
  }

  // Change diagonal elements to one.
  for (int i = 0; i < n; i++)
  {
    identity.setValue(1,i,i);
  }
  return identity;
}

// Square identity matrix
const Matrix eye(int n)
{
  /** @brief Square identity matrix constructor
   * @details 
   *    Creates a square identity matrix with n rows and n columns. 
   * @param n The desired number of rows and columns in the square 
   *          identity matrix
   * @return   Square identity matrix with desired number of row and columns
   */

  // Call eye with dimensions n x n.
  return eye(n,n);
}
//////////////////////////////////////////////////////////


////////////// ******** Random Matrix ********* //////////
const Matrix rand(int row, int col)
{
  /** @brief Random matrix with desired number of rows
   *        and columns
   * @details Creates matrix with dimensions row x col,
   *      elements are given values sampled from uniform
   *      distribution between 0 and 1.
   * @param row Number of rows
   * @param col Number of columns
   * @return Matrix object with random element values
   *         between 0 and 1 on a uniform distribution.
   */

  // Initialise matrix
  Matrix R(row,col);
  // Updated element values with random values
  for (int i=0; i < R.mSize; i++)
  {
     R.mData[i] = rand();
  }
  return R;
}

const Matrix rand(int n)
{
  /** @brief Random square matrix
   * @details Creates square matrix with dimensions n x n,
   *      elements are given values sampled from uniform
   *      distribution between 0 and 1.
   * @param n Number of rows of square matrix
   * @return Matrix object with random element values
   *         between 0 and 1 on a uniform distribution.
   */

  // Call rand with dimension n x n
  return rand(n,n);
}
//////////////////////////////////////////////////////////


//////// ******** Arithmatic operators ******** //////////
// definition of + between two Matrixs
Matrix operator+(const Matrix& m1, 
                        const Matrix& m2)
{
  /** @brief Addition operator for matrix-matrix addition
   * @details Overloads addition operator to support addition
   * of two matrix objects with the same dimensions
   * @param m1 First matrix object
   * @param m1 Second matrix object
   * @return Matrix object that results from element-wise
   *        addition of m1 and m2.
   */

  //  add the Matrixs
  //  If demnsion are equal proceed with addition
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    // Initialise addition matrix
    Matrix ma(m1.mRow, m1.mCol);
    for (int i=0; i<m1.mSize; i++)
	  { 
      // Set addition matrix element values to
      // addition of two input matrices
	    ma.mData[i] = m1.mData[i] + m2.mData[i];
	  }
    return ma;
  }
  // If input matrices not of same dimensions throw error
  else
  {
    throw Exception("Different dimensions", "Matrix add - Matrices are different dimensions\n");
  }
}

Matrix operator+(const Matrix& m, 
                        const double& a)
{
  /** @brief Addition between doubles and Matrix object
   * @details Overloads addition operator to support 
   *      element-wise addition of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by addition
   *         of double value
   */

  //  add double to matrix
  Matrix ma(m);
  for (int i=0; i<m.mSize; i++)
  {
    ma.mData[i] = m.mData[i] + a;
  }
  return ma;
}

Matrix operator+(const double& a, const Matrix& m)
{
   /** @brief Addition between doubles and Matrix object
   * @details Overloads addition operator to support 
   *      element-wise addition of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by addition
   *         of double value
   */
  return m+a;
}

// definition of - between two Matrices
Matrix operator-(const Matrix& m1, 
                        const Matrix& m2)
{
    /** @brief Subraction operator for matrix-matrix subtraction
   * @details Overloads subtraction operator to support subtraction
   * of two matrix objects with the same dimensions
   * @param m1 First matrix object
   * @param m1 Second matrix object
   * @return Matrix object that results from element-wise
   *        subtraction of m1 and m2.
   */

  // Check that number of rows and columns are equal
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    // Initialise matrix subtract object
    Matrix ms(m1.mRow, m1.mCol);
    for (int i=0; i<m1.mSize; i++)
	  {
      // Element-wise subtraction
	    ms.mData[i] = m1.mData[i] - m2.mData[i];
	  }
    return ms;
  }
  else
  {
    throw Exception("Different dimensions", "Matrix subtract - Matrices are different dimensions\n");
  }
}

Matrix operator-(const Matrix& m, const double& a)
{
    /** @brief Subtraction between doubles and Matrix object
   * @details Overloads subtraction operator to support 
   *      element-wise subtraction of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by subtraction of double
   *        from element values.
   */
  return m + (-a);
}

Matrix operator-(const double& a, const Matrix& m)
{
  /** @brief Subtraction between doubles and Matrix object
   * @details Overloads subtraction operator to support 
   *      element-wise subtraction of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by addition of
   *         double value to negative element values,
   *         i.e. a - m[0,0]
   */
  return (-m)+a;
}







// definition of matrix procture
Matrix operator*(const Matrix& m1, const Matrix& m2)
{
  /** @brief Matrix-matrix multiplication
   * @details Overloads multiplication operator for
   * matrix-matrix multiplication. Supports non-square matrices.
   * @param m1 Left matrix
   * @param m2 Right matrix
   * @return Matrix multiplication of m1 * m2
   */

//  check matrices are of the same dimension,
//  i.e m1 number of col = m2 number of row
  if (m1.mCol == m2.mRow)
    {
      
      int n = m1.mCol;
      double prod = 0;
      Matrix mp(m1.mRow, m2.mCol);
      // Perform matrix multiplication
      // Computes row column dot product
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
      throw Exception("Different dimensions", "Matrix product - Number of columns in the first matrix"
                        " does not equal number of rows in the second matrix\n");
    }
}


// definition of multiplication between a Matrix and a scalar
Matrix operator*(const Matrix& m, const double& a)
{
 /** @brief Multiplication between doubles and Matrix object
   * @details Overloads multiplication operator to support 
   *      element-wise multiplication of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by multiplication of
   *         double value to element values,
   */
//  create a Matrix of the same dimension as m with entries equal to a*m
  Matrix mp(m.mRow, m.mCol);
  // Element-wise multiplication
  for (int i=0; i<m.mSize; i++)
    {
      mp.mData[i] = a * m.mData[i];
    }

  return mp;
}


// definition of multiplication between a scalar and a Matrix
Matrix operator*(const double& a, const Matrix& m)
{
   /** @brief Multiplication between doubles and Matrix object
   * @details Overloads multiplication operator to support 
   *      element-wise multiplication of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by multiplication of
   *         double value to element values,
   */

//  create a Matrix of the same dimension as m with entries equal to a*m
  Matrix mp(m.mRow, m.mCol);
// Element-wise multiplication
  for (int i=0; i<m.mSize; i++)
    {
      mp.mData[i] = a * m.mData[i];
    }
  return mp;
}



// definition of division of a Matrix by a scalar
Matrix operator/(const Matrix& m, const double& a)
{
   /** @brief Division between doubles and Matrix object
   * @details Overloads division operator to support 
   *      element-wise division of a double number.
   * @param m Matrix object
   * @param a double number
   * @return Matrix object with element updated by multiplication of
   *         double value to element values,
   */
  
  // Throw error if division by zero occurs
  if (a == 0.0)
  {
     throw Exception("Division by 0", "Attempt to divide by zero in the"
                       " Matrix-scalar division operator");
  }
//  create a Matrix of the same length as v with entries equal to v/a
  Matrix md(m.mRow, m.mCol);
// Element-wise division
  for (int i=0; i< m.mSize; i++)
    {
      md.mData[i] = m.mData[i] / a;
    }

  return md;
}


// dot product
double dot(const Matrix& m1, const Matrix& m2)
{
  /** @brief Dot product between two matrices (vectors) with
   *        one column or row
   * @details Performs the vector dot product on matrix or vector
   * objects (i.e those with one column or row)
   * @param m1 First vector (matrix)
   * @param m2 Second vector (matrix)
   * @return Dot product of vector m1 and vector m2.
   */

  // Check that at least row or column is equal to one
  // And that the two vectors have the same dimensions.
  if ((m1.mCol == m1.mSize || m1.mRow == m1.mSize) &&
        (m2.mCol == m2.mSize || m2.mRow == m2.mSize)
        && (m1.mSize == m2.mSize))
  {
    double prod = 1;
    
    // Sum of element-wise products
    for (int i = 0; i < m1.mSize; i++)
    {
      prod *= m1.mData[i] * m2.mData[i];
    }
    return prod;
  }
  // Throw exception if one row or column dimension not equal to one
  // and if vectors are not the same length.
  else
  {
    throw Exception("dot - incorrect dimensions","");
  }
}




///////////////////////////////////////////////////////////////////////////////

/////////////// *********** Unary operator *************** ////////////////////
// definition of the unary operator -
Matrix operator-(const Matrix& m)
{
  /** @brief Return element-wise negative of input matrix
   * @details Overloads unary operator to perform element-wise
   *      negation.
   * @param m Input matrix
   * @return Negative of input matrix
   */
//  create a Matrix w with entries equal to -v
  Matrix ms(m.mRow, m.mCol);
  // Element-wise negation
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
  /**
   * @brief Getter: Allows MATLAB like access of
   *        Matrix elements (1-indexing)
   * @details Allows matrix elements to be accessed 
   *      and indexed in a MATLAB fashion
   *      i.e. A_ij = A(i,j)
   * @param i the row of desired element (in 1-indexing)
   * @param j the column of the desired element (in 1-indexing)
   * @return The element at positon (i,j) (in 1-indexing)
   */

  // Throw exception if coordinates inputed are outside
  // of matrix dimensions
  if (i <= 0 || i > mRow || j <= 0 || j > mCol)
  {
    throw Exception("(i,j) - Out of range", "The element coordinates given when"
                      " getting an element value are out of range (use 1-indexing).");
  }
  // Return element
  return mData[mCol * (i-1) + (j-1)];
}


//////////// ******** Print functions ************ /////////////////
void print(const Matrix& m)
{
  /** @brief Prints matrix object in formatted way
   * @param m Matrix to be printed
   */

  // Set print precision to 3 dp
  std::cout << std::fixed << std::setprecision(3);
  // Loop through matrix elements
  for (int i = 0; i < m.mRow; i++) 
  {
    // If start value negative add extra space
    std::cout << "[";
    for (int j = 0; j < m.mCol; j++) 
    {
      // If non start add spacing between elements
      std::cout << " ";
      std::cout << (m.getValue(i,j) >= 0 ? " " : "") << m.getValue(i,j);
    }
    // Add end bracket at end of row
      std::cout << " ]" << std::endl;
    }
    // Return precision to 6 dp.
  std::cout << std::endl << std::setprecision(6);
}

std::ostream& operator<<(std::ostream& output, const Matrix& m) 
{
  /** @brief Print with in std::cout; 
   *  @details Overloads << operator so that matrices and correctly
   *          formatted when using std::cout
   * @param output current stream to cout
   * @param m matrix to be printed
   * @return Updated stream with matrix formating
   */

  // Set print precision to 3 dp
  output << std::fixed << std::setprecision(3);
  // Loop through matrix elements
  for (int row = 0; row < m.mRow; row++) 
  {
    // If start value negative add extra space
    output << "[";
    for (int col = 0; col < m.mCol; col++) 
    {
      // Add space in between elements
      output << " ";
      // If non negative add extra space
      output << (m.getValue(row,col) >= 0 ? " " : "") << m.getValue(row,col);
    }
    // Add ] at the end of the row
    output << " ]" << std::endl;
    
  }
  // Return precision to 6 dp
  output << std::setprecision(6);
  return output;
}
///////////////////////////////////////////////////////////////////



///////////// ****** Assignment operator ******** /////////////////
// definition of Matrix operator =
Matrix& Matrix::operator=(const Matrix& m)
{
  /** @brief Assignment operator
   * @details Assigns one matrix object to another, however,
   *      they are not linked. Essential creates a copy.
   * @param m Matrix object we wish to "copy" to self/this.
   * @return Updated self/this object with values of m.
   */

// check both Matrices have same dimension
  if (mRow != m.mRow && mCol != m.mCol)
    {
      throw Exception("Dimension mismatch", "Matrix assignment operator - Matrices"
                          " have different dimensions");
    }
  else
  {
    // Element-wise update mData entries
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
  /** @brief Transposes matrix
   * @details Create transpose matrix by swapping rows
   *          and columns.
   * @return transposed matrix
   */
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
  /** @brief Transposes matrix
   * @details Create transpose matrix by swapping rows
   *          and columns.
   * @param m Matrix to be transposed
   * @return transposed matrix
   */
  
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
  /** @brief Calculates matrix norm
   * @details Calculates induces matrix 1-norm and 
   *      2-norm (for symmetric matrices). If matrices are 
   *      vectors (1 row or col) calculates vector norm.
   * @param p choice of norm. For matrices this is 1 or 2-norm,
   *    for vectors this can be any positive integer
   * @return Norm calculate for self/this
   * 
   */
  // If ve
  if ((mRow==1 || mCol ==1) && (p > 0))
  {
    // Calculate vector norm
    double temp, norm_val;
    norm_val = 0.0;
    // Sum of abs value of element to the power p
    for (int i=0; i < mSize; i++)
      {
        temp = fabs(mData[i]);
        norm_val += pow(temp, p);
      }
    
    return pow(norm_val, 1.0/((double) (p)));
  }
  // if two norm selected
  if (p==2)
  {
    // Can only find singular values of symmetric matrices
    if (transpose(*this) == *this)
    {
      Matrix eigVals = eigenVal(*this);
      double max_eig = eigVals.colMax(0);
      return max_eig;
    }
    throw Exception("2-Norm : Non-symmetric", "The 2-norm is not implemented"
                    " for non-symmetric matrices.");
  }
  else if (p==1)
  {
    // Compute induced 1 norm, i.e max abs column sum
    double max_sum = 0;
    double sum;
    // Row vector get max abs row sum
    if (mRow == 1)
    {
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
    }
    else{
    // Calculate each abs column sum
    for (int j=0; j < mCol; j++)
    {
      sum = 0;
      for (int i=0; i < mRow; i++)
      {
        sum += fabs(getValue(i,j));
      }
      // Check if current sum greater than previous
      if (sum > max_sum)
      {
        max_sum = sum;
      }
    }
    }
  // Return max
    return max_sum;
  }
  else
  {
  throw Exception("Matrix::norm - norm selected not available","Choice of norm not available,"
                       " please choose from the 1-norm (p=1), 2-norm (p=2), infinity-norm (p=\"inf\")"
                       " or the Frobenius-norm (p=\"frob\")");
  }
}

double Matrix::norm(std::string p) const
{
  /** @brief Calculates matrix norm
   * @details Calculates induces matrix frobenius norm and 
   *      infinity -norm .
   * @param p choice of norm. 
   * @return Norm calculate for self/this
   * 
   */
  //  sum of squared elements 
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
  // Find max abs row sum
  else if (p=="inf")
  {
    double max_sum = 0;
    double sum;
    // Calculate abs row sum
    for (int i=0; i < mRow; i++)
    {
      sum = 0;
      for (int j=0; j < mCol; j++)
      {
        sum += fabs(getValue(i,j));
      }
      // Check if current row sum greater than previous
      if (sum > max_sum)
      {
        max_sum = sum;
      }
    }
    return max_sum;
  }
  else
  {
  throw Exception("Matrix::norm - norm selected not available","Choice of norm not available,"
                       " please choose from the 1-norm (p=1), 2-norm (p=2), infinity-norm (p=\"inf\")"
                       " or the Frobenius-norm (p=\"frobenius\")");
  }
}

/////////////////////////////////////////////////////////////////////////////////////


/////////////////// ********** Matrix Size ************ //////////////////////////
// return size of a Matrix
int* size(const Matrix& m)
{
  /** @brief Gets number of rows and columns
   *  @param m Matrix from which number of rows and columns
   *          want to be known
   *  @return Array of size 2 [mRow, mCol];
   */
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
  /** @brief Find maximum value in row, by default this is the absolute max
   * @details Find maximum value in row, and optionally from a starting column
   *       (default is the first column).
   * @param row Row to check for max
   * @param col Optional: Column to start checking for max from
   * @param abs_true Optional: return absolution max or regular max
   * @return Row max value.
   */
  double max_val ;
  double temp;
  max_val = getValue(row,col);
  // Check each element in row
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
    // if current greater than previous
    if (max_val < temp)
    {
      max_val = temp;
    }
  }
  return max_val;
}


double Matrix::colMax(int col, int row, bool abs_true) const
{
  /** @brief Find maximum value in column, by default this is the absolute max
   * @details Find maximum value in column, and optionally from a starting row
   *       (default is the first row).
   * @param col Column to check for max
   * @param row Optional: Row to start checking for max from
   * @param abs_true Optional: return absolution max or regular max
   * @return Row max value.
   */

  double max_val;
  double temp;
  // Check each element of column
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
    // if current greater than previous
    if (max_val < temp)
    {
      max_val = temp;
    }
  }
  return max_val;
}

int Matrix::argRowMax(int row, int col, bool abs_true) const
{
  /** @brief Find argument of maximum value in row, by default this is the absolute max
   * @details Find row coordinate with maximum value in column, 
   *      and optionally from a starting column (default is the first column).
   * @param row Row to check for max
   * @param col Optional: Columns to start checking for max from
   * @param abs_true Optional: return absolution max or regular max
   * @return Column positon of row max value.
   */


  int arg_max = col;
  double max_val ;
  double temp;
  // Check each element in row
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
  /** @brief Find argument of maximum value in column, by default this is the absolute max
   * @details Find column coordinate with maximum value in column, 
   *      and optionally from a starting row (default is the first row).
   * @param col Column to check for max
   * @param row Optional: Rows to start checking for max from
   * @param abs_true Optional: return absolution max or regular max
   * @return Row position of column max value.
   */


  int arg_max = row;
  double max_val;
  double temp;
  // check each element in column
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
  /** @brief Equality check between two matrices
   *  @details Overloads equaility boolean operator
   *           to check equality element-wise between two
   *          matrices
   * @param m1 First matrix
   * @param m2 Second matrix
   * @return True or false
   */

  // Check if dimension are equal
  if (m1.mRow == m2.mRow && m1.mCol == m2.mCol)
  {
    // Check equality element wise
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
    /** @brief non-equality check between two matrices
   *  @details Overloads non-equality boolean operator
   *           to check non-equality element-wise between two
   *          matrices
   * @param m1 First matrix
   * @param m2 Second matrix
   * @return True or false
   */
  return !(m1 == m2);
}
//////////////////////////////////////////////////////////////////////////////////

////////////////// *************** Swap elements ************* ////////////////////
void Matrix::swap(int row1, int col1, int row2, int col2) const
{
  /** @brief Swap two elements in matrix
   * @param row1 Element 1 row coordinate
   * @param col1 Element 1 column coordinate
   * @param row2 Element 2 row coordinate
   * @param col2 Element 2 column coordinate
   */
  // Create temp and swap values
  double temp;
  temp = getValue(row1,col1);
  setValue(getValue(row2,col2),row1,col1);
  setValue(temp,row2,col2);
}


void Matrix::swapRow(int row1, int row2) const
{
   /** @brief Swap two columns in matrix
   * @param row1 Row 1 position
   * @param row2 Row 2 position
   */
  // Create and swap every element in rows
  double temp;
  for (int j = 0; j < mCol; j++){
    temp = getValue(row1,j);
    setValue(getValue(row2,j),row1,j);
    setValue(temp,row2,j);
  }
}

void Matrix::swapCol(int col1, int col2) const
{
  /** @brief Swap two columns in matrix
   * @param col1 Row 1 position
   * @param col2 Row 2 position
   */
  // Create and swap every element in column
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

std::tuple <Matrix, Matrix, Matrix, int> lu(const Matrix& A)
{
  /** @brief Compute LU decompositon with partial pivoting of square matrix
   * @param A Target matrix for LU decompositon
   * @return Matrix P - permutation matrix, Matrix L - lower triangular matrix
   *       Matrix U - upper triangular matrix, int sign - number of permutation made
   */

  // Get number of rows and columns
  int m = A.mRow;
  int n = A.mCol;

  // If not equal throw exception
  if (m != n)
  {
    throw Exception("LU - Non-square matrix input", "LU decomposition is"
                    " designed to receive only square matrices as an input.");
  }

  

  // Initialise P, L, U matrices
  Matrix U(A);
  Matrix L = eye(n);
  Matrix P = eye(n);

  int i;
  double temp;
  int sign = 1;

  // Compute LU decomposition
  for (int k = 0; k < m - 1; k++)
  {
    // Find argument of col max from diagonal down
    i = U.argColMax(k,k);
    if (i != k)
    {
      // Permute rows of L
      for (int j = k; j < m; j++)
      {
        U.swap(k,j,i,j);
      }
     // Permute rows of U
      for (int j = 0; j < k - 1; j++)
      {
        L.swap(k,j,  i,j);
      }
      // Permute rows of P
      P.swapRow(k,i);
      // Record sign change
      sign *= -1;
    }

    // Perform gaussian elimination to achieve 
    // Upper and lower triangular matrices
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

  return {P, L, U,sign};
}
/////////////////////////////////////////////////////////////////////////////////

/////////////// ************ Determinant ************* /////////////////////////
double Matrix::det() const
{
  /** @brief Calculate determinant of square matrix
   * @details Find the determinant of a matrix using LU decomposition
   * @return the determinant of self/this
   * 
   */

  // Check if dimension are equal
  if (mRow != mCol)
  {
    throw Exception("Matrix::det() - Non-square matrix","Matrix entered is non-square."
                      " Square matrix required to find determinant.");
  }
  // Get LU decomposition
  int sign;
  Matrix P(mRow,mCol); Matrix L(mRow,mCol); Matrix U(mRow,mCol);
  std::tie(P,L,U,sign) = lu(*this);
  // Calculate product of diagonals of L and U
  double prod = 1;
  for (int i = 0; i < mRow; i++)
  {
    prod *= L.getValue(i,i) * U.getValue(i,i);
  }
  // Adjust determinant for number of permutation made in LU Decomposition
  return sign*prod;

}

double det(const Matrix& A)
{
  /** @brief Calculate determinant of square matrix
   * @details Find the determinant of a matrix using LU decomposition
   * @param A target matrix to find determinant of
   * @return the determinant of self/this
   * 
   */
  return A.det();
}

///////////////////////////////////////////////////////////////////////////////////

//////////////////////////// ********** QR Factorisation ************ /////////////
std::tuple <Matrix, Matrix> qr(const Matrix& A, bool reduced)
{
  /** @brief Calculates QR factorisation of matrix
   *  @details Performs a qr factorisation on m-by-n matrix A such that
      A = Q*R. The factor R is an m-by-n upper triangular matrix and Q is an
        m-by-m orthogonal matrix.
      @param A Target matrix for QR factorisation
      @param reduce Optional: by default false, if true returns reduced QR
      @return Matrix Q - orthogonal matrix, Matrix R - upper triangular
   */

  
  int m = A.mRow;
  int n = A.mCol;
  // Check columns greater than rows
  if (m < n)
  {
    throw Exception("qr - Columns greater than rows", "This implementation of QR factorisation"
                      " requires that the number of rows\n is greater than or equal to"
                      " the number of columns.");
  }
  
  // Initialise matrices
  Matrix R(A);
  Matrix Q = eye(m);

  // Householder vectors
  Vector x(m);
  Vector w(m);
  Vector u(n);
  
  // Sign and norm
  double g;
  double s;;
  for (int k=0; k < n; k++)
  {

    // Construct Householder x(1:k-1) = 0; x(k:m) = R
    for (int i = 0; i < k; i++)
    {
      x.setValue(0,   i);
    }

    for (int i = k; i < m; i++)
    {
      x.setValue(R.getValue(i,k),   i);
    }
    
    // Get sign(x(k)) * norm(x)
    g = (double) ((x.getValue(k) > 0) - (x.getValue(k) < 0)) * norm(x);
    
    // x(k) = x(k) + sign(x(k)) * norm(x)
    x.setValue(x.getValue(k) + g, k);
    s = norm(x);
    if (s != 0)
    {
      // Perform Householder reflection
      w = (x / s);
      u = 2 * (R.T() * w); // H 
      R = R - (w * u.T());  // R = H*R
      Q = Q - 2 * Q * (w * w.T());  // Q = Q*H
    }
  } 
  // Return reduced QR if reduced == true
  if (reduced)
  {
    // Create reduce matrices
    Matrix Q_red(m,n);
    Matrix R_red(n,n);
    // Copy reduced values from original QR matrices
    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        Q_red.setValue(Q.getValue(i,j),i,j);
        if (i < n)
        {
          R_red.setValue(R.getValue(i,j),i,j);
        }
      }
    }
    return {Q_red, R_red};
  }
  return {Q,R};
}

///////////////////////////////////////////////////////////////////////////////////

/////////////////// ************* Least Squares *********** ///////////////////////

Matrix lsq(const Matrix& A, const Matrix& b)
{
  /** @brief Finds least squares solutio to overdetermined 
   *        system |Ax-b|_2
   *  @details It attempts to solve the least squares solution x 
   * that minimizes norm(b-A*x,2). Uses reduced QR factorisation 
   * and then gaussian elimination.
   * @param A Overdetermined matrix in Ax = b
   * @param b Vector b
   * @return Solution x to least-squares problem
   */
  // Initialise matrices
  int m = A.mRow; int n = A.mCol;
  Matrix Q(m,n); Matrix R(n);
  std::tie(Q,R) = qr(A,true);
  // Solve R x = Q^T *b
  return  gaussianElimination(R, Q.T() * b );
}

//////////////////// ********** Upper Hessenberg ********** ///////////////////////
Matrix hessenbergReduction(const Matrix& A)
{
  /** @brief Finds upper Hessenberg reduction of matrix A
   *  @details Compute the upper Hessenberg reduction of A
   *      using Householder reflections
   * @param A Target matrix for upper Hessenberg reduction
   * @return Matrix H, upper Hessenberg reduction of A.
   * 
   */

  // Initialise Hessenberg and Householder vectors
  int m = A.mRow;
  Matrix H(A);
  Vector x(m);
  Vector w(m);
  Vector u(m);
  
  double g;
  double s;
  if (m > 2)
  {
    // Loop up to second last column of A
    for (int k=0; k < m-2 ; k++)
    {
      
      // Initialise x with zero up to current row
      for (int i = 0; i < k+1; i++)
      {
        x.setValue(0,   i);
      }
      // Assign values below current row with values from H(k+1:m,k+1)
      for (int i = k+1; i < m; i++)
      {
        x.setValue(H.getValue(i,k),   i);
      }
      
      // Calculate   sign(x(k+1)) * norm(x,2)
      g = (double) ((x.getValue(k+1) >= 0) - (x.getValue(k+1) < 0)) * norm(x);

      // x(k+1) = x(k+1) + sign(x(k+1)) * norm(x,2)
      x.setValue(x.getValue(k+1) + g, k+1);

      s = norm(x);
      if (s != 0)
      {
        w = (x / (double) s);
        H = H - 2. * w * (w.T() * H);
        H = H - 2. * (H * w) * w.T(); 
      }

    } 
    return H;
  }
  else
  {
    return H;
  }
}
//////////////////////////////////////////////////////////////////////////////////

////////////////////// ********* Eigenvalues by QR Algorithm ************ ///////////////////////
Matrix eigenVal(const Matrix& A, double tol)
{
  /** @brief Finds eigenvalues of symmetric matrix using QR Algorithm
   *  @details Uses QR Algorithm with shift, hessenberg reduction, and 
   *          deflation to finds eigenvalues.
   *  @param A Target matrix for eigenvalues
   *  @param tol Tolerance for when eigenvalue has been found.
   *            Determined by H(k,k-1) approx 0.
   * @return eigVals: Vecotr of eigenvalues of A.
   * 
   */

  // Check dimension equal
  if (A.mRow != A.mCol)
  {
    throw Exception("eigenSystem - Non-square matrix","Please enter a square matrix");
  }

  // Initialise eigVal vector and Hessenberg
  int m = A.mRow;
  Vector eigVals(m);
  Matrix H_old(A);   // Store hessenberg from previous iteration

  // Variables for Wilkinson shift
  double delta;
  double mu;


    for (int k = m; k > 0; k--)
    {
      if (k == 1)
      {
        eigVals.setValue(H_old.getValue(0,0),0);
        return eigVals;
      }

      // Initialise matrices for QR Alg
      Matrix I = eye(k);
      Matrix H(k);
      Matrix Q(k);
      Matrix R(k);

      // Deflate current H with values from H_old
      for (int i = 0; i < k; i++)
      {
        for (int j = 0; j < k; j++)
        {
          H.setValue(H_old.getValue(i,j), i,j);
        }
      }

      // Perform reduction to upper hessenberg
      H = hessenbergReduction(H);

      
      while (fabs(H.getValue(k-1,k-2)) > tol)
      {
        // Wilkinson shift
        delta = (H.getValue(k-2,k-2) - H.getValue(k-2,k-1)) / 2.;
        mu = H.getValue(k-2,k-1) - (double) ((delta > 0) - (delta < 0)) * pow(H.getValue(k-1,k-1),2) 
                / (fabs(delta) + pow(pow(delta,2)+pow(H.getValue(k-1,k-1),2),0.5));
        
        // QR factorisation with shift
        std::tie(Q,R) = qr(H - mu*I);

        // QR Alg with shift back
        H = R*Q + mu*I;
      }

      // Once converges add eigenvalue to vector
      eigVals.setValue(H.getValue(k-1,k-1),k-1);

      // Deflation
      for (int i = 0; i < k-1; i++)
      {
        for (int j = 0; j < k-1; j++)
        {
          H_old.setValue(H.getValue(i,j), i,j);
        }
      }
    }
  
  return eigVals;
}


///////////// ********* Gaussian Elimination ************ ////////////////////////
Matrix gaussianElimination(const Matrix& A, const Matrix& b)
{
  /** @brief Solves full-rank square system Ax = b
   *   with Gaussian Elimination.
   *  @param A Matrix of the linear system
   *  @param b vector of the linear system
   *  @return vector x: solution to the linear system
   */


  // Check if A is singular
  if (fabs(det(A)) < 1e-10)
  {
    throw Exception("gaussianElimination - No solutions","Matrix is singular");
  }
  // Check if A is square
  if (A.mCol > A.mRow)
  {
    throw Exception("gaussianElimination - Underdetermined", "Cannot find a solution");
  }
  // Check if b is a vector
  if (b.mCol != 1)
  {
    throw Exception("gaussianElimination - b has more than one column",
                           "Enter b as vector or one column matrix");
  }

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
        b.setValue(b.getValue(i,0)- M * b.getValue(k,0), i, 0); 

    }
  }
  // Declare x
  Matrix x(A.mRow,1);
  double sum;

  for (int k = A.mCol - 1; k >= 0; k--)
  {
    sum = 0;
    for (int i = k; i < A.mCol; i++)
    {
      sum += A.getValue(k,i) * x.getValue(i,0);
    }
    x.setValue((1. / A.getValue(k,k)) * (b.getValue(k,0) - sum), k, 0);
  }
  return x;
}