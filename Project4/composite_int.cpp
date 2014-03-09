/*
Will Spurgin
11/24/2013
High Performance Scietific Computing
MATH 3316
*/

//inclusions
#include <iostream>

using namespace std;

/*
Composite Integration function will be a Gaussian quadrature
rule with 3 points per subinterval to achieve O(h^6) accuracy.

The function 'fun' below must have the following syntax: y = fun(x)

Usage, F = composite_int(fun, a, b, n)

Where
    fun = integrand
    a = lower limit of integration
    b = upper limit of integration
    n = number of subintervals

    F = value of numerical integral
*/
double composite_int(double (*f)(const double), const double a, const double b,
    const int n)
{
    if (b < a) 
    {
        cerr << "error: illegal interval, b < a" << endl;
        return 0.0;
    }
    if (n < 1) {
        cerr << "error: illegal number of subintervals, n < 1" << endl;
        return 0.0;
    }

    //subinterval width
    double h = (b-a)/n;

    //set weights and nodes
    double x1 = -.774596669241483; //to avoid calling the sqrt function
    double x2 = 0.0;
    double x3 = .774596669241483;
    double w1 = 5.0/9.0;
    double w2 = 8.0/9.0;
    double w3 = 5.0/9.0;

    double F = 0.0;

    double xmid, node1, node2, node3;
    //Iterate through subintervals
    for(int i = 0; i < n; i++)
    {

        // find evaluations points
        xmid  = a + (i+0.5)*h;
        node1 = xmid + 0.5*h*x1;
        node2 = xmid + 0.5*h*x2;
        node3 = xmid + 0.5*h*x3;

        // add approximation to final result
        F += 0.5*h*(w1*f(node1) + w2*f(node2) + w3*f(node3));
    }

    return F;
}
