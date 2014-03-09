/*
Will Spurgin
11/20/2013
High Performance Scientific Computing
MATH 3316
*/

/*
The following is an adaptation of Daniel R. Reynolds 'test_Gauss2.cpp'.
The majority of the code written below is his with only slight modifications
in order to test the 'composite_int' function. Let credit be given where
credit is due.
*/

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;

// function prototypes
double composite_int(double (*f)(const double), const double a, const double b,
  const int n);

// Integrand
const double c=3.0;
const double d=20.0;
double f(const double x) {
  return (exp(c*x) + sin(d*x));
}


// This routine tests the Gauss-2 method on a simple integral
int main(int argc, char* argv[]) {

  // limits of integration
  double a = -5.0;
  double b = 2.0;

  // true integral value
  double Itrue = 1.0/c*(exp(c*b) - exp(c*a)) - 1.0/d*(cos(d*b)-cos(d*a));
  printf("\n True I = %22.16e\n", Itrue);


  // test the composite_int which is Gauss-3 rule
  cout << "\n Composite-int (Gauss-3) rule:\n";
  cout << "     n              R              relerr    conv rate\n";
  cout << "  ---------------------------------------------------\n";
  int n[] = { 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120 };
  int ntests=10;

  // iterate over n values, computing approximations, error, convergence rate
  double Iapprox, olderr, relerr=0.0;
  for (int i=0; i<ntests; i++) {

    printf("   %6i", n[i]);
    olderr = relerr;
    Iapprox = composite_int(f, a, b, n[i]);
    relerr = fabs(Itrue-Iapprox)/fabs(Itrue);
    if (i == 0) {
      printf("  %22.16e  %7.1e\n", Iapprox, relerr);
    } else {
      printf("  %22.16e  %7.1e   %f\n", Iapprox, relerr, 
         (log(olderr) - log(relerr))/(log(1.0/n[i-1]) - log(1.0/n[i])));
    }
    
  }
  cout << "  ---------------------------------------------------\n";

}