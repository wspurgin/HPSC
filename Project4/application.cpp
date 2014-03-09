/*
Will Spurgin
11/29/2013
High Performance Scientific Computing
MATH 3316
*/

#include <iostream>
#include <cmath>

using namespace std;

//function prototypes

//main carbon concentration function
double carbon(const double x, const double t,const double rtol,
    const double atol);

//Single variable function
double f(const double t);

//root finding function
double fd_newton(double (*f)(const double), double x,
        int maxit, double tol, double alpha);



int main(int argc, char** argv)
{
    double alpha = pow(2.0, -26);
    int maxit = 100;
    double tol = 1e-6;
    cout << "Calculating the time it takes for the carbon concentration"
        << ", which initially is 0.1\% to get " << endl << "to 4\% through "
        << "carburizing using a gas with a carbon concentration"
        << " of 10\% at a distance of 3mm from the metal" << endl;
    double t = fd_newton(f, 120000, maxit, tol, alpha);
    cout << "Approximated time for .04 = C(3e-3,t) or root finding problem "
        << "0 = C(3e-3,t) - .04 is:" << endl << "t = ";
        printf("%6fs\n", t);
    return 0;
}

double f(const double t)
{
    return (carbon(3e-3, t, 1e-14, 1e-15) - .04);
}
