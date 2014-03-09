/*
Will Spurgin
9/26/2013
HPSC - Nonlinear solvers
MATH 3316
*/

#include <cmath>
#include <iostream>

using namespace std;

double fd_newton(double (*f)(const double), double x,
		int maxit, double tol, double alpha)
{
	double c;
	if(maxit < 1)
	{
		cerr << "Warning: the maximum iterations is set to: " << maxit
			<< " Resetting to 100" << endl;
		maxit = 100;
	}
	for(int i = 0; i < maxit; i++)
	{
		c = x - f(x)/((f(x + alpha) - f(x))/alpha);
		double err = abs(c - x);
		cout << " iter: " << i << " | x = " << x << " | c = " << c 
			<< " | absolute value of f(x) = " << abs(x) << " | err = " << err
			<< endl;
		if(err < tol)
			break;
		x = c;
	}
	return c;	
}
//End of File
