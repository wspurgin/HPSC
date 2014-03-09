/* Daniel R. Reynolds
   SMU Mathematics
   Math 3316
   9 August 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mat.h"
#include <iostream>

// This file implements the operations defined in the Mat class


///// General class routines /////

// constructor (initializes values to 0.0)
Mat::Mat(long int m, long int n) {
  nrows = m;
  ncols = n;
  data = new double[n*m];
  own_data = true;
  for (long int i=0; i<n*m; i++)  data[i] = 0.0;
}


// empty constructor (sets data to point to existing array)
Mat::Mat(long int m, long int n, double *data_array) {
  nrows = m;
  ncols = n;
  data = data_array;
  own_data = false;
}

// copy constructor
Mat::Mat(const Mat& A) {
  nrows = A.Rows();
  ncols = A.Cols();
  data = new double[ncols*nrows];
  own_data = true;
  for (long int i=0; i<nrows*ncols; i++)  data[i] = A.data[i];

}

// write myself to stdout
int Mat::Write() const {
  // return with failure if data array isn't allocated
  if (data == NULL) {
    fprintf(stderr, "Mat:Write error, empty data array\n");
    return 1;
  }

  // print data to screen 
  for (long int i=0; i<nrows; i++) {
    for (long int j=0; j<ncols; j++)
      printf("  %g",data[IDX(i,j,nrows)]);
    printf("\n");
  }
  return 0;
}


// write myself to a file
int Mat::Write(const char *outfile) const {
  // return with failure if data array isn't allocated
  if (data == NULL) {
    fprintf(stderr, "Mat:Write error, empty data array\n");
    return 1;
  }

  // return with failure if 'outfile' is empty
  if (strlen(outfile) < 1) {
    fprintf(stderr, "Mat::Write error, empty outfile\n");
    return 1;
  }

  // open output file
  FILE *fptr = NULL;
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    fprintf(stderr, "Mat::Write error, unable to open %s for writing\n",outfile);
    return 1;
  }

  // print data to file
  for (long int i=0; i<nrows; i++) {
    for (long int j=0; j<ncols; j++)
      fprintf(fptr, "  %g",data[IDX(i,j,nrows)]);
    fprintf(fptr, "\n");
  }

  // close output file and return
  fclose(fptr);
  return 0;
}


// streaming output operator
std::ostream& operator<<(std::ostream &out, Mat &data) {
  // print data to "out" stream and return
  for (long int i=0; i<data.Rows(); i++) {
    for (long int j=0; j<data.Cols(); j++)
      out << "  "<<data(i,j);
    out<<std::endl;
  }
  return out;
}


// access my data -- access value by location in array, 
// gets reference to value so that it can be changed
double& Mat::operator()(long int i, long int j) {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat() error: empty data array\n");
    return data[0];
  }

  // check that (i,j) is a valid entry
  if (i<0 || i>=nrows || j<0 || j>=ncols) {
    fprintf(stderr,"Mat() error: (%li,%li) outside matrix bounds\n",i,j);
    fprintf(stderr,"  returning first entry to avoid seg. fault\n");
    return data[0];
  }
  return data[IDX(i,j,nrows)];
}


// access my data -- access value by location in vector, 
// gets reference to value so that it can be changed
double& Mat::operator()(long int i){
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat() error: empty data array\n");
    return data[0];
  }

  // check that (i) is a valid entry
  if (i<0 || i>=nrows*ncols) {
    fprintf(stderr,"Mat() error: (%li) outside vector bounds\n",i);
    fprintf(stderr,"  returning first entry to avoid seg. fault\n");
    return data[0];
  }
  return data[i];
}



///// Arithmetic operations defined on a given Mat /////


// C = a*C + b*B
int Mat::LinearSum(double a, double b, const Mat &B) {
  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    fprintf(stderr,"Mat::LinearSum error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li, Mat 2 is %li x %li\n", 
	    nrows, ncols, B.nrows, B.ncols);
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || B.data == NULL) {
    fprintf(stderr, "Mat::LinearSum error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = a*data[i] + b*B.data[i];
  
  // return success
  return 0;
}


// C = a*A + b*B
int Mat::LinearSum(double a, const Mat &A, double b, const Mat &B) {
  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols || B.nrows != nrows || B.ncols != ncols) {
    fprintf(stderr,"Mat::LinearSum error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li,  Mat 3 is %li x %li\n",
	    nrows, ncols, A.nrows, A.ncols, B.nrows, B.ncols);
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL || B.data == NULL) {
    fprintf(stderr, "Mat::LinearSum error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = a*A.data[i] + b*B.data[i];
  
  // return success
  return 0;
}


//   C = C+a  (adds scalar a to my data)
int Mat::Add(double a) {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Add error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] += a;
  
  // return success
  return 0;
}


//   C = C.*A (multiplies my data by y, component-wise)
int Mat::Mul(const Mat &A) {
  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    fprintf(stderr,"Mat::Mul error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    nrows, ncols, A.nrows, A.ncols);
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL) {
    fprintf(stderr, "Mat::Mul error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] *= A.data[i];
  
  // return success
  return 0;
}


//   C = a*C  (scales my data by scalar a)
int Mat::Mul(double a) {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Mul error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] *= a;
  
  // return success
  return 0;
}


//   C = A  (copies A into C)
int Mat::Copy(const Mat &A) {
  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    fprintf(stderr,"Mat::Copy error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    nrows, ncols, A.nrows, A.ncols);
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL) {
    fprintf(stderr, "Mat::Copy error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] = A.data[i];
  
  // return success
  return 0;
}


//   C = a  (sets all entries of C to the scalar a)
int Mat::Const(double a) {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Const error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] = a;
  
  // return success
  return 0;
}


// C = C.^p
int Mat::Power(double p) {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Power error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] = pow(data[i], p);
  
  // return success
  return 0;
}


// Cij = |Cij|
int Mat::Abs() {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Abs error: empty data array\n");
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] = fabs(data[i]);
  
  // return success
  return 0;
}


// Cij = Cji
int Mat::Trans() {
  // check that data is not NULL
  if (data == NULL) {
    fprintf(stderr, "Mat::Trans error: empty data array\n");
    return 1;
  }

  // perform operation in place if matrix is square
  if (nrows == ncols) {
    double tmp;
    for (long int i=0; i<nrows; i++)
      for (long int j=0; j<i; j++) {
	tmp = data[IDX(i,j,nrows)];
	data[IDX(i,j,nrows)] = data[IDX(j,i,nrows)];
	data[IDX(j,i,nrows)] = tmp;
      }

  // otherwise we need a new data array to ensure a clean transpose
  } else {

    // create temporary data array, and copy transposed data over 
    double *newdata = new double[nrows*ncols];
    for (long int i=0; i<nrows; i++)
      for (long int j=0; j<ncols; j++)
	newdata[IDX(j,i,ncols)] = data[IDX(i,j,nrows)];

    // copy newdata values into existing data array
    for (long int i=0; i<nrows*ncols; i++)
      data[i] = newdata[i];

    // delete temporary data array
    delete[] newdata;
    
    // swap matrix dimensions
    int tmp = nrows;
    nrows = ncols;
    ncols = tmp;
  }
  
  // return success
  return 0;
}


// fill in a vector with linearly spaced data
int Mat::Linspace(double a, double b) {
  // ensure that this is a row or column vector
  if (nrows != 1 && ncols != 1) {
    fprintf(stderr,"Mat::Linspace error, matrix must be a row/column vector\n");
    fprintf(stderr,"  dimensions are %li x %li\n", nrows, ncols);
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = a + (b-a)/(nrows*ncols-1)*i;
  return 0;
}


// fill in a vector with logarithmically spaced data
int Mat::Logspace(double a, double b) {
  // ensure that this is a row or column vector
  if (nrows != 1 && ncols != 1) {
    fprintf(stderr,"Mat::Logspace error, matrix must be a row/column vector\n");
    fprintf(stderr,"  dimensions are %li x %li\n", nrows, ncols);
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = pow(10.0, a + (b-a)/(nrows*ncols-1)*i);
  return 0;
}


// fill in a vector with uniformly-distributed random numbers in [0,1]
int Mat::Random() {
  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = random() / (pow(2.0,31.0) - 1.0);
  return 0;
}


// fill in a matrix as the identity
int Mat::Eye() {
  // ensure that this is a square matrix
  if (nrows != ncols) {
    fprintf(stderr,"Mat::Eye error, matrix must be square\n");
    fprintf(stderr,"  dimensions are %li x %li\n", nrows, ncols);
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = 0.0;
  for (long int i=0; i<nrows; i++)
    data[IDX(i,i,nrows)] = 1.0;
  return 0;
}


///// New matrix creation operations /////


// C = A+B
Mat& Mat::operator+(const Mat &B) {
  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    fprintf(stderr,"Mat::operator+ error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    nrows, ncols, B.nrows, B.ncols);
    Mat *C = new Mat(0,0);
    return *C;
  }

  // create new Mat for output, and do operation
  Mat *C = new Mat(nrows,ncols);
  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] + B.data[i];

  // return result
  return *C;
}


// C = A-B
Mat& Mat::operator-(const Mat &B) {
  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    fprintf(stderr,"Mat::operator- error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    nrows, ncols, B.nrows, B.ncols);
    Mat *C = new Mat(0,0);
    return *C;
  }

  // create new Mat for output, and do operation
  Mat *C = new Mat(nrows,ncols);
  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] - B.data[i];

  // return result
  return *C;
}


// C = A*B
Mat& Mat::operator*(const Mat &B) {
  // create a pointer to the result
  Mat *C = NULL;

  // determine if either matrix is a scalar
  bool A_scalar = ((nrows==1) && (ncols==1));
  bool B_scalar = ((B.nrows==1) && (B.ncols==1));

  // scalar-times-matrix
  if (A_scalar) {

    // create new Mat for output, and do operation
    C = new Mat(B.nrows,B.ncols);
    for (long int i=0; i<nrows*ncols; i++)  
      C->data[i] = data[0] * B.data[i];
    
  // scalar-times-matrix
  } else if (B_scalar) {
    
    // create new Mat for output, and do operation
    C = new Mat(nrows,ncols);
    for (long int i=0; i<nrows*ncols; i++)  
      C->data[i] = data[i] * B.data[0];
        
  // normal matrix product
  } else {

    // check that array sizes are acceptable
    if (B.nrows != ncols) {
      fprintf(stderr,"Mat::operator* error, inner dimension mismatch\n");
      fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	      nrows, ncols, B.nrows, B.ncols);
      C = new Mat(0,0);
      return *C;
    }
    
    // create new Mat for output, and do operation
    C = new Mat(nrows,B.ncols);
    for (long int i=0; i<nrows; i++) 
      for (long int j=0; j<B.ncols; j++)
	for (long int k=0; k<ncols; k++)
	  C->data[IDX(i,j,nrows)] += data[IDX(i,k,nrows)] * B.data[IDX(k,j,B.nrows)];

  }

  // return result
  return *C;
}


// C = B*a
Mat& Mat::operator*(const double b) {
  // create a pointer to the result
  Mat *C = new Mat(nrows,ncols);

  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] * b;
    
  // return result
  return *C;
}


// C = a*B
Mat& operator*(const double a, const Mat &B) {
  // create a pointer to the result
  Mat *C = new Mat(B.Rows(),B.Cols());
  C->LinearSum(0.0, a, B);
    
  // return result
  return *C;
}


// C = A^T
Mat Mat::T() {
  // create new Mat for output, and do operation
  Mat *C = new Mat(ncols,nrows);
  for (long int j=0; j<ncols; j++)  
    for (long int i=0; i<nrows; i++)  
      C->data[IDX(j,i,ncols)] = data[IDX(i,j,nrows)];

  // return result
  return *C;
}


// column accessor routine (does not alocate its own data, only returns a 
// Mat object whose data array points to an existing Mat column).
Mat Mat::AccessColumn(long int j) {
  // check that requested column exists
  if (j < 0 || j > ncols) {
    fprintf(stderr,"Mat::AccessColumn error, requested column does not exist\n");
    fprintf(stderr,"  j is %li, but matrix only has %li columns\n", j, ncols);
    Mat *C = new Mat(0,0);
  }

  // create new matrix, pointing to column j of this matrix
  Mat *C = new Mat(nrows, 1, &(data[IDX(0,j,nrows)]));

  // return result
  return *C;
}


///// Scalar output operators on matrices /////


// vector 2-norm, square root of the matrix Frobenius norm
double Mat::Norm() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::Norm error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double sum=0.0;
  for (long int i=0; i<nrows*ncols; i++)  sum += data[i]*data[i];
  return sqrt(sum);
}


// vector infinity norm, largest absolute value entry of a matrix
double Mat::MaxNorm() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::MaxNorm error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int i=0; i<nrows*ncols; i++)  
    mx = (mx > fabs(data[i])) ? mx : fabs(data[i]);
  return mx;
}


// row vector one norm, column vector infinity norm, matrix infinity norm
double Mat::InfNorm() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::InfNorm error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int i=0; i<nrows; i++) {
    double sum=0.0;
    for (long int j=0; j<ncols; j++) 
      sum += fabs(data[IDX(i,j,nrows)]);
    mx = (mx > sum) ? mx : sum;
  }
  return mx;
}
 

// row vector infinity norm, column vector one norm, matrix one norm
double Mat::OneNorm() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::OneNorm error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int j=0; j<ncols; j++) {
    double sum=0.0;
    for (long int i=0; i<nrows; i++)
      sum += fabs(data[IDX(i,j,nrows)]);
    mx = (mx > sum) ? mx : sum;
  }
  return mx;
}


// min Cij
double Mat::Min() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::Min error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double mn=data[0];
  for (long int i=0; i<nrows*ncols; i++)  
    mn = (mn < data[i]) ? mn : data[i];
  return mn;
}


// max Cij
double Mat::Max() {
  // check that my data is allocated
  if (data == NULL) {
    fprintf(stderr,"Mat::Max error, data not allocated\n");
    return -1.0;
  }
  
  // perform operation
  double mx=data[0];
  for (long int i=0; i<nrows*ncols; i++)  
    mx = (mx > data[i]) ? mx : data[i];
  return mx;
}


// equivalence-checking operator
bool Mat::operator==(Mat &A)  {
  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    fprintf(stderr,"Mat::operator== error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    nrows, ncols, A.nrows, A.ncols);
    return false;
  }

  // perform operation and return
  bool equal=true;
  for (long int j=0; j<ncols; j++)  
    for (long int i=0; i<nrows; i++)
      equal &= (A(i,j) == data[IDX(i,j,nrows)]);
  return equal;
}


///////// independent Mat routines /////////


// create a row vector of linearly spaced data
Mat Linspace(double a, double b, long int n) {
  // create vector, call in-place Linspace routine, and return
  Mat *x = new Mat(1,n);
  x->Linspace(a,b);
  return *x;
}


// create a row vector of logarithmically spaced data
Mat Logspace(double a, double b, long int n) {
  // create vector, call in-place Logspace routine, and return
  Mat *x = new Mat(1,n);
  x->Logspace(a,b);
  return *x;
}


// create a matrix of uniformly-distributed random data in [0,1]
Mat Random(long int m, long int n) {
  // create matrix, call in-place Random routine, and return
  Mat *R = new Mat(m,n);
  R->Random();
  return *R;
}


// create an n by n identity matrix
Mat Eye(long int n) {
  // create matrix, call in-place Eye routine, and return
  Mat *I = new Mat(n,n);
  I->Eye();
  return *I;
}


// dot-product of x and y
double Dot(Mat &x, Mat &y) {
  // check that array sizes match
  if (y.Rows() != x.Rows() || y.Cols() != x.Cols()) {
    fprintf(stderr,"Dot error, matrix size mismatch\n");
    fprintf(stderr,"  Mat 1 is %li x %li,  Mat 2 is %li x %li\n", 
	    x.Rows(), x.Cols(), y.Rows(), y.Cols());
    return 0.0;
  }
  
  // perform operation and return
  double sum=0.0;
  for (long int j=0; j<x.Cols(); j++)  
    for (long int i=0; i<x.Rows(); i++)  
      sum += x(i,j)*y(i,j);
  return sum;
}


// performs backwards substitution on the linear system U*x = b, filling in the input Mat x
int BackSub(Mat &U, Mat &x, Mat &b) {
  // check that matrix sizes match
  if (U.Rows() != b.Rows() || U.Rows() != U.Cols() || 
      b.Cols() != 1 || x.Rows() != U.Rows() || x.Cols() != 1) {
    fprintf(stderr,"BackSub error, illegal matrix/vector dimensions\n");
    fprintf(stderr,"  Mat is %li x %li,  sol is %li x %li,  rhs is %li x %li\n", 
	    U.Rows(), U.Cols(), x.Rows(), x.Cols(), b.Rows(), b.Cols());
    return 1;
  }
  
  // copy b into x
  x = b;

  // perform column-oriented Backwards Substitution algorithm
  for (long int j=U.Rows()-1; j>=0; j--) {

    // solve for this row of solution
    x(j) = x(j)/U(j,j);

    // update all subsequent rhs
    for (long int i=0; i<j; i++)
      x(i) -= U(i,j)*x(j);

  }

  // return success
  return 0;
}


// performs backwards substitution on the linear system U*x = b, returning x as a new Mat; 
// leaves U and b untouched
Mat BackSub(Mat &U, Mat &b) {
  // check that matrix sizes match
  if (U.Rows() != b.Rows() || U.Rows() != U.Cols() || b.Cols() != 1) {
    fprintf(stderr,"BackSub error, illegal matrix/vector dimensions\n");
    fprintf(stderr,"  Mat is %li x %li,  rhs is %li x %li\n", 
	    U.Rows(), U.Cols(), b.Rows(), b.Cols());
    Mat *x = new Mat(0,0);
    return *x;
  }
  
  // create new Mat for output
  Mat *x = new Mat(U.Rows(),1);

  // call existing BackSub routine for computations
  if (BackSub(U, *x, b) != 0)
    fprintf(stderr,"BackSub Warning: error in BackSub call\n");

  // return result
  return *x;
}


// solves a linear system A*x = b, filling in the input Mat x
int Solve(Mat &A, Mat &x, Mat &b) {

  // create temporary variables
  long int i, j, k, p, n;
  double tmp, Amax;

  // check that matrix sizes match
  if (A.Rows() != b.Rows() || A.Rows() != A.Cols() || 
      b.Cols() != 1 || x.Rows() != A.Rows() || x.Cols() != 1) {
    fprintf(stderr,"Solve error, illegal matrix/vector dimensions\n");
    fprintf(stderr,"  Mat is %li x %li,  sol is %li x %li,  rhs is %li x %li\n", 
	    A.Rows(), A.Cols(), x.Rows(), x.Cols(), b.Rows(), b.Cols());
    return 1;
  }

  // determine maximum absolute entry in A (for singularity check later)
  Amax = A.MaxNorm();

  // perform Gaussian elimination to convert A,b to an upper-triangular system
  n = A.Rows();
  for (k=0; k<n-1; k++) {   // loop over diagonals

    // find the pivot row p
    p=k;
    for (i=k; i<n; i++)  
      if (fabs(A(i,k)) > fabs(A(p,k)))  
	p=i;

    // swap rows in A
    for (j=k; j<n; j++) {
      tmp = A(p,j);
      A(p,j) = A(k,j);
      A(k,j) = tmp;
    }

    // swap rows in b
    tmp = b(p);
    b(p) = b(k);
    b(k) = tmp;

    // check for singular matrix
    if (fabs(A(k,k)) < 1.e-13*Amax) {
      fprintf(stderr,"Solve error: numerically singular matrix!\n");
      return 1;
    }

    // perform elimination on remaining submatrix of A using row k
    for (j=k+1; j<n; j++) 
      for (i=k+1; i<n; i++) 
	A(i,j) = A(i,j) - A(i,k)/A(k,k)*A(k,j);

    // perform elimination on remainder of b using row k
    for (i=k+1; i<n; i++)  b(i) -= A(i,k)/A(k,k)*b(k);

  }

  // check for singularity at end (only need to check final diagonal entry)
  if (fabs(A(n-1,n-1)) < 1.e-13*Amax) {
    fprintf(stderr,"Solve error: numerically singular matrix!\n");
    return 1;
  }

  // perform Backwards Substitution on result
  if (BackSub(A, x, b) != 0) {
    fprintf(stderr,"Solve error in BackSub call\n");
    return 1;
  }

  // return success
  return 0;
}


// solves a linear system A*x = b, returning x as a new Mat
Mat Solve(Mat &A, Mat &b) {
  // check that matrix sizes match
  if (A.Rows() != A.Rows() || A.Rows() != A.Cols() || b.Cols() != 1) {
    fprintf(stderr,"Solve error, illegal matrix/vector dimensions\n");
    fprintf(stderr,"  Mat is %li x %li,  rhs is %li x %li\n", 
	    A.Rows(), A.Cols(), b.Rows(), b.Cols());
    Mat *x = new Mat(0,0);
    return *x;
  }

  // create new Mat for output
  Mat *x = new Mat(A.Rows(),1);

  // call existing Solve routine for computations
  if (Solve(A, *x, b) != 0)
    fprintf(stderr,"Solve Warning: error in Solve call\n");

  // return result
  return *x;
}




// end of file
