/*
Will Spurgin
9/30/2013
HPSC - Nonlinear solvers
MATH 3316
*/

#include <cmath>
#include <iostream>

using namespace std;

double f(double V);
double df(double V);
double newton(double (*f)(const double), double (*df)(const double), double x,
		int maxit, double tol);
double fd_newton(double (*f)(const double), double x,
		int maxit, double tol, double alpha);
double bisection(double (*f)(const double), 
		 double a, double b, int maxit, double tol);

int main(int argc, char** argv)
{
	int maxit = 10e3;
	double tol = 10e-5;

	//Part (a) bisection
	double a = 0.5;
	double b = 1.0;
	cout << "The Approximate root for f(v) using Bisection is: " << bisection(
		f, a, b, maxit, tol) << endl << endl;

	//Part (b) newton's method
	double x = 0.75;
	maxit = 30;
	cout << "The Approximate root for f(v) using newton's method is: "
		<< newton(f, df, x, maxit, tol) << endl << endl;

	//Part (c) newton forward-finite difference
	double alpha = 2e-6;
	cout << "The Approximate root for f(v) using newton forward-finite "
		<<" difference alpha = 2e-6 is: "
		<< fd_newton(f, x, maxit, tol, alpha) << endl << endl;

	alpha = 2e-26;
	cout << "The Approximate root for f(v) using newton forward-finite "
		<<" difference alpha = 2e-26 is: "
		<< fd_newton(f, x, maxit, tol, alpha) << endl << endl;

	alpha = 2e-48;
	cout << "The Approximate root for f(v) using newton forward-finite "
		<<" difference alpha = 2e-48 is: "
		<< fd_newton(f, x, maxit, tol, alpha) << endl << endl;
}

double f(double V)
{
	double q = 1.6022e-19;
	double v_oc = 0.8;
	double k_b = 1.3806e-23;
	double T = 315;
	double x = q/(k_b*T);

	return exp(V*x)*(1 + (V*x)) - exp(v_oc*x); 
}

double df(double V)
{
	double q = 1.6022e-19;
	double k_b = 1.3806e-23;
	double T = 315;
	double x = q/(k_b*T);

	return x*exp(V*x)*(1 + V*x) + x*exp(V*x);
}
//End of File
