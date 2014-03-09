/* Daniel R. Reynolds
   SMU Mathematics
   Math 3316
   8 October 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "mat.h"

using namespace std;

// function prototypes
inline double f(const double x);
double lagrange(Mat &x, Mat &y, double z);


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
  
  // evaluate the polynomial at the points z, storing in p
  Mat p(n,1);
  for (int i=0; i<n; i++) 
    p(i) = lagrange(x,y,z(i));

  // output errors at each point
  cout << endl << "interpolant and error using " << n+1 << " nodes:\n";
  cout << "      z        f(z)               p(z)              error\n";
  for (int i=0; i<n; i++) 
    printf("   %6.3f   %16.13f   %16.13f   %g\n",
	   z(i), f(z(i)), p(i), fabs(f(z(i))-p(i)));


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
  
  // evaluate the polynomial at the points z, storing in p
  Mat p2(n,1);
  for (int i=0; i<n; i++) 
    p2(i) = lagrange(x2,y2,z2(i));

  // output errors at each point
  cout << endl << "interpolant and error using " << n+1 << " nodes:\n";
  cout << "      z        f(z)               p(z)              error\n";
  for (int i=0; i<n; i++) 
    printf("   %6.3f   %16.13f   %16.13f   %g\n",
	   z2(i), f(z2(i)), p2(i), fabs(f(z2(i))-p2(i)));

} // end routine


// function to interpolate
inline double f(const double x) {
  return (exp(x*x));
}

