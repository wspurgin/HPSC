/*
Will Spurgin
9/23/2013
HPSC - Nonlinear solvers
MATH 3316
*/

#include <cmath>
#include <iostream>

using namespace std;

double newton(double (*f)(const double), double (*df)(const double), double x,
		int maxit, double tol)
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
		c = x - f(x)/df(x);
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
