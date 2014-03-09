/*
Will Spurgin
10/26/2013
High Performance Scientific Computing
Math 3316
Project 3 - interpolation
*/

#include <math.h>
#include <iostream>
#include <chrono>
#include "mat.h"

using namespace std;

// function prototypes
inline double f(const double x);
double newton_eval(Mat &x, Mat &y, double z);
double newton_coeffs(Mat& x, Mat& y, Mat& c);
double lagrange(Mat &x, Mat &y, double z);
void compare(int n, int m);

int main(int argc, char** argv)
{
	int n[] = {10, 25, 50};
	int m[] = {1000, 10000, 100000};
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
			compare(n[i], m[j]);
	}	
	cout << endl << "END OF COMPARISONS" << endl;
	return 0;
}

inline double f(const double x)
{
	return (cos(3*x*x));
}

void compare(int n, int m)
{
	chrono::time_point<chrono::system_clock> stime, ftime;
  	chrono::duration<double> runtime;

	//Tests of n nodes and m evaluation nodes
	Mat x = Linspace(-1.0, 1.0, n+1);
	Mat y(n+1, 1);
	for(int i = 0; i < n+1; i++)
		y(i) = f(x(i));
	Mat z = Linspace(-1.0, 1.0, m+1);

  	Mat p(m, 1);
  	//lagrange
  	stime = chrono::system_clock::now();
	for(int i = 0; i < m; i++)
		p(i) = lagrange(x, y, z(i));
	ftime = chrono::system_clock::now();
	runtime = ftime - stime;
	printf("Using Lagrange Interpolation with %d nodes and %d evaluation nodes,"
		"it took %.4fs to calculate and evaluate the polynomial.\n", n, m,
		 runtime.count());

	//newton
  	stime = chrono::system_clock::now();
  	Mat c(n+1, 1);
  	newton_coeffs(x, y, c);
	for(int i = 0; i < m; i++)
		p(i) = newton_eval(x, c, z(i));
	ftime = chrono::system_clock::now();
	runtime = ftime - stime;
	printf("Using Newton Interpolation with %d nodes and %d evaluation nodes,"
		"it took %.4fs to calculate and evaluate the polynomial.\n", n, m,
		runtime.count());

	cout << endl;
}
