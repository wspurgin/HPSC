/*
Will Spurgin
9/27/2013
HPSC - Nonlinear Solvers
MATH 3316
*/

#include <iostream>

using namespace std;

//prototypes
double f(const double x);
double df(const double x);
double newton(double (*f)(const double), double(*df)(const double), double x,
		int maxit, double tol);

int main(int argc, char** argv)
{
	cout << "Approximating the root with newton's method..." << endl;
	cout << "x = -2, tolerance = 1e-2" << endl;
	double x = -2;
	double tol = 1e-2;
	int maxit = 20;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl;
	cout << "x = -2, tolerance = 1e-6" << endl;
	tol = 1e-6;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl << endl;
	cout << "x = -2, tolerance = 1e-10" << endl;
	tol = 1e-10;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl << endl;
	cout << "x = 2, tolerance = 1e-2" << endl;
	tol = 1e-2;
	x = 2;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl << endl;
	cout << "x = 2, tolerance = 1e-6" << endl;
	tol = 1e-6;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl << endl;
	cout << "x = 2, tolerance = 1e-10" << endl;
	tol = 1e-10;
	cout << endl << "The approximated root is " << newton(f, df, x, maxit, tol)
		<< endl << endl;
}

double f(const double x)
{
	double c = x*x*x - 2*x*x - 3*x;
	return c;
}

double df(const double x)
{
	double c = x*x*3 - 4*x - 3;
	return c;
}
