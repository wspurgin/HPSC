/* Daniel R. Reynolds
   SMU Mathematics
   Math 3316
   9 August 2013 */

#ifndef MAT_DEFINED__
#define MAT_DEFINED__

// Inclusions
#include <math.h>
#include <iostream>
#include <ostream>

// Macros (column-oriented storage), pi to full double precision
#define IDX(i,j,m) ((m)*(j)+(i))
#define PI 3.141592653589793


// This defines a simple matrix class
class Mat {

 private:

  ///// Contents /////
  long int nrows;
  long int ncols;
  double *data;
  bool own_data;

 public:

  ///// General class routines /////
  // constructor (initializes values to 0.0)
  Mat(long int m, long int n);

  // empty constructor (sets data to point to existing array)
  Mat(long int m, long int n, double *data_array);

  // copy constructor
  Mat(const Mat&);

  // destructor
  ~Mat() { if (own_data) delete[] data; };

  // write myself to stdout
  int Write() const;

  // write myself to a file
  int Write(const char *outfile) const;

  // streaming output routine
  friend std::ostream& operator<<(std::ostream &strm, Mat &data);

  // returns my number of rows
  long int Rows() const { return nrows; };

  // returns my number of columns
  long int Cols() const { return ncols; };

  // access my data -- access value by location in matrix, 
  // gets handle to value so that it can be changed as well
  double& operator()(long int i, long int j);

  // access my data -- access value by location in vector, 
  // gets handle to value so that it can be changed as well
  double& operator()(long int i);


  ///// Arithmetic operations defined on a Mat /////

  // in-place operations (C is the matrix calling the routine)
  int LinearSum(double a, double b, const Mat &B);                // C = a*C + b*B
  int LinearSum(double a, const Mat &A, double b, const Mat &B);  // C = a*A + b*B
  int Add(const Mat &A) { return LinearSum(1.0, 1.0, A); };       // C = C+A
  int Add(double a);                                              // C = C+a
  int Sub(const Mat &A) { return LinearSum(1.0, -1.0, A); };      // C = C-A
  int Sub(double a) { return Add(-a); };                          // C = C-a
  int Mul(const Mat &A);                                          // C = C.*A
  int Mul(double a);                                              // C = a*C
  int Power(double p);                                            // C = C.^p
  int Copy(const Mat &A);                                         // C = A
  int Const(double a);                                            // C = a
  int Abs();                                                      // Cij = |Cij|
  int Trans();                                                    // C = C^T
  int Linspace(double a, double b);                               // C = linear span
  int Logspace(double a, double b);                               // C = log10 span
  int Random();                                                   // Cij is random
  int Eye();                                                      // C=I

  // in-place shortcut operators
  int operator+=(const Mat &A) { return Add(A); };  // C += A
  int operator+=(double a) { return Add(a); };      // C += a
  int operator-=(const Mat &A) { return Sub(A); };  // C -= A
  int operator-=(double a) { return Sub(a); };      // C -= a
  int operator*=(const Mat &A) { return Mul(A); };  // C *= A
  int operator*=(double a) { return Mul(a); };      // C *= a
  int operator^=(double p) { return Power(p); };    // C = C.^p
  int operator=(const Mat &A) { return Copy(A); };  // C = A
  int operator=(double a) { return Const(a); };     // C = a

  // new matrix creation operations (C is the output, A and B are the operands)
  Mat& operator+(const Mat &B);                     // C = A+B
  Mat& operator-(const Mat &B);                     // C = A-B
  Mat& operator*(const Mat &B);                     // C = A*B
  Mat& operator*(const double a);                   // C = B*a
  Mat T();                                          // C = A^T
  Mat AccessColumn(long int j);                     // C = A(:,j)

  // Scalar output operators on matrices
  double Norm();      // sqrt(sum_ij Cij^2) (sqrt matrix Frobenius norm, vector 2-norm)
  double MaxNorm();   // max_ij |Cij|       (vector inf-norm)
  double InfNorm();   // max_i sum_j |Cij|  (row vector 1-norm, col vector inf-norm)
  double OneNorm();   // max_j sum_i |Cij|  (row vector inf-norm, col vector 1-norm)
  bool operator==(Mat &A);

  // other scalar Matrix statistics
  double Min();       // min_ij Cij
  double Max();       // max_ij Cij

};  // end Mat


// extra scalar*matrix routine to support scalar-first order
Mat& operator*(const double a, const Mat &B);                   // C = a*B

// creates a row vector of n linearly spaced values from a through b
Mat Linspace(double a, double b, long int n);

// creates a row vector of n logarithmically spaced values from 10^a through 10^b
Mat Logspace(double a, double b, long int n);

// creates a random n by m matrix (uniform random distribution)
Mat Random(long int m, long int n);

// creates an n by n identity matrix
Mat Eye(long int n);

// dot-product routine
double Dot(Mat &x, Mat &y);

// performs backwards substitution on the linear system U*x = b, filling in the input Mat x
int BackSub(Mat &U, Mat &x, Mat &b);

// performs backwards substitution on the linear system U*x = b, returning x as a new Mat
Mat BackSub(Mat &U, Mat &b);

// solves a linear system A*x = b, filling in the input Mat x
int Solve(Mat &A, Mat &x, Mat &b);

// solves a linear system A*x = b, returning x as a new Mat
Mat Solve(Mat &A, Mat &b);


#endif
