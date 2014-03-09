/* Daniel R. Reynolds
   SMU Mathematics
   Math 3316
   16 September 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;

// function prototypes
double f(const double x);
double bisection(double (*f)(const double), 
		 double a, double b, int maxit, double tol);


// This routine tests the function bisection.cpp on a nonlinear 
// function.  It uses an initial interval [-5,5], and solves for 
// the root to a tolerance of 1e-4. 
int main(int argc, char* argv[]) {

  // set the initial interval, max iteration count, and the tolerance
  double a = -5.0;
  double b = 5.0;
  int maxit = 1000000;
  double tol = 1e-4;

  // call bisection to compute the root, and output result to screen
  double x = bisection(f,a,b,maxit,tol);
  cout << endl << " The approximate root is " << x << endl;

}


// Root-finding residual function for part 1 of project 2.
double f(const double x) {
  return (10.0 + x*x*x - 12.0*cos(x));
}
