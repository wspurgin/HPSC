/*
Will Spurgin
11/24/2013
High Performance Scietific Computing
MATH 3316
*/

#include <iostream>
#include <cmath>
#ifndef PI
    #define PI 3.141592653589793
#endif

using namespace std;

//prototypes
int adaptive_int(double (*f)(const double), const double a, const double b,
    const double rtol, const double atol, double &R, int &n);

double f_(const double x) { return (exp(-pow(x,2))); }

double erf(const double y, const double rtol, const double atol)
{
    double R; int n;
    int errors = adaptive_int(f_, 0, y, rtol, atol, R, n); 
    if(!errors)
        return (2/sqrt(PI)*R);
    else
    {
        cerr << "The error did not converge with y = " << y << endl;
        return (2/sqrt(PI)*R);
    }
}

double carbon(const double x, const double t,const double rtol,
    const double atol)
{
    double C_s = 0.1;
    double C_0 = 0.001;
    double D = 5e-11;
    double denom = (sqrt(4*D*t));
    if(denom == 0)
        return 0.0;
    double z = (x/denom);
    return (C_s - (C_s - C_0)*erf(z, rtol, atol));
}

