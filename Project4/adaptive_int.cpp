/*
Will Spurgin
11/24/2013
High Performance Scientific Computing
MATH 3316
*/

#include <iostream>
#include <cmath>

using namespace std;

//prototype
double composite_int(double (*f)(const double), const double a, const double b,
  const int n);

int adaptive_int(double (*f)(const double), const double a, const double b,
    const double rtol, const double atol, double &R, int &n)
{

    //set up the interval iterator
    int iterator;
    if(rtol >= 1e-2)
        iterator = 3;
    else if (rtol >= 1e-4)
        iterator = 4;
    else if (rtol >= 1e-6)
        iterator = 10;
    else if (rtol >= 1e-8)
        iterator = 14;
    else
        iterator = 35;

    //primary read
    int maxit = 700;
    int intervals = 1;
    double previous, current;

    n = intervals;

    previous = composite_int(f, a, b, intervals);
    bool converged = false;

    for(int i = 0; i < maxit; i++)
    {
        intervals += iterator;
        current = composite_int(f, a, b, intervals);

        if(abs(current - previous) < (rtol*current + atol))
        {
            converged = true;
            break;
        }
        previous = current;
        n += intervals;
    }
    R = current;

    if(converged)
        return 0; 
    else
        return 1;
}
