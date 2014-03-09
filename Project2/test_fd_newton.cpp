/*
Will Spurgin
9/26/2013
HPSC - Nonlinear solvers
MATH 3316
*/

#include <iostream>

using namespace std;

double f(const double x);
double fd_newton(double (*f)(const double), double x,
		int maxit, double tol, double alpha);

int main(int argc, char** argv)
{
	cout << "Approximating root using Forward-Difference Newton's method" << endl;
	cout << "x = -2, tolerance = 1e-2, alpha = 2e-6" << endl;
	double x = -2;
	double tol = 1e-2;
	int maxit = 20;
	double alpha = 2e-6;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-2, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;
	
	cout << "x = -2, tolerance = 1e-2, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-6, alpha = 2e-6" << endl;
	alpha = 2e-6;
	tol = 1e-6;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-6, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-6, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-10, alpha = 2e-6" << endl;
	alpha = 2e-6;
	tol = 1e-10;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-10, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = -2, tolerance = 1e-10, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-2, alpha = 2e-6" << endl;
	x = 2;
	tol = 1e-2;
	alpha = 2e-6;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-2, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;
	
	cout << "x = 2, tolerance = 1e-2, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-6, alpha = 2e-6" << endl;
	alpha = 2e-6;
	tol = 1e-6;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-6, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-6, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-10, alpha = 2e-6" << endl;
	alpha = 2e-6;
	tol = 1e-10;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-10, alpha = 2e-26" << endl;
	alpha = 2e-26;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	cout << "x = 2, tolerance = 1e-10, alpha = 2e-48" << endl;
	alpha = 2e-48;
	cout << endl << "The approximated root is " << fd_newton(f, x, maxit, tol, alpha)
		<< endl;

	return 0;
}

double f(const double x)
{
	double c = x*x*x - 2*x*x - 3*x;
	return c;
}
