/*
Will Spurgin
10/26/2013
High Performance Scientific Computing
MATH 3316
Project 3 - Interpolation
*/

/*
The following is an adaptation of Daniel R. Reynolds test_lagrange.cpp,
the majority of the written code is his. The only adaptation is the replacement
of the references to lagrange function with the newton_eval function and the
addition of the calculations for the newton coefficients.
*/

#include <iostream>
#include <math.h>
#include "mat.h"

using namespace std;

// function prototypes
inline double f(const double x);
double newton_eval(Mat &x, Mat &y, double z);
double newton_coeffs(Mat& x, Mat& y, Mat& c);


// This routine tests the function lagrange.cpp over the interval [-1,1].
int main(int argc, char* argv[]) {

	  ///////////////
	// first, test with 7 nodes
	int n = 6;                          // set n
	Mat x = Linspace(-1.0, 1.0, n+1);   // set nodes
	Mat y(n+1,1);                       // initialize data
	for (int i=0; i<=n; i++)            // fill data
	  y(i) = f(x(i));

	// set evaluation points z as midpoints between nodes
	double dx = 2.0/n;                  // set node spacing
	Mat z = Linspace(-1.0+dx/2.0, 1.0-dx/2.0, n);
	  
	Mat c(n+1, 1);
	if(newton_coeffs(x, y, c) != 0)
  	{
	  // evaluate the polynomial at the points z, storing in p
	  Mat p(n,1);
	  for (int i=0; i<n; i++) 
	    p(i) = newton_eval(x, c ,z(i));

	  // output errors at each point
	  cout << endl << "interpolant and error using " << n+1 << " nodes:\n";
	  cout << "      z        f(z)               p(z)              error\n";
	  for (int i=0; i<n; i++) 
	    printf("   %6.3f   %16.13f   %16.13f   %g\n",
		   z(i), f(z(i)), p(i), fabs(f(z(i))-p(i)));
	}


	///////////////
  	// repeate test with 15 nodes
  	n = 15;                              // set n
  	Mat x2 = Linspace(-1.0, 1.0, n+1);   // set nodes
  	Mat y2(n+1,1);                       // initialize data
  	for (int i=0; i<=n; i++)             // fill data
    	y2(i) = f(x2(i));

  	// set evaluation points z as midpoints between nodes
  	dx = 2.0/n;                  // set node spacing
  	Mat z2 = Linspace(-1.0+dx/2.0, 1.0-dx/2.0, n);

  	Mat c2(n+1, 1);
  	if(newton_coeffs(x2, y2, c2) != 0)
  	{
	  	// evaluate the polynomial at the points z, storing in p
	  	Mat p2(n,1);
	  	for (int i=0; i<n; i++) 
	    	p2(i) = newton_eval(x2,c2,z2(i));

	  	// output errors at each point
	  	cout << endl << "interpolant and error using " << n+1 << " nodes:\n";
	  	cout << "      z        f(z)               p(z)              error\n";
	  	for (int i=0; i<n; i++) 
	    	printf("   %6.3f   %16.13f   %16.13f   %g\n",
		   	z2(i), f(z2(i)), p2(i), fabs(f(z2(i))-p2(i)));
	}
	return 0;
} // end routine


// function to interpolate
inline double f(const double x) {
  return (exp(x*x));
}

