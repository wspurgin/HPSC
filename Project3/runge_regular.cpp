/*
Will Spurgin
10/26/2013
High Performance Scientific Computing
MATH 3316
Project 3 - Interpolation
*/

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include "mat.h"

using namespace std;

inline double f(const double x, const double y);

void run_tests(int n, int m, Mat& a, Mat& b);

double lagrange2D(Mat &x, Mat &y, Mat &z, double a, double b);

int main(int argc, char** arv)
{
	Mat a = Linspace(-5, 5, 201);
	Mat b = Linspace(-5, 5, 201);

	//write a and b to disk
	a.Write("avals.txt");
	b.Write("bvals.txt");
	run_tests(10, 10, a, b);
	run_tests(20, 20, a, b);

	Mat runge(201,201);
	for(int i = 0; i < 201; i++)
		for(int j = 0; j < 201; j++)
			runge(i, j) = f(a(i), b(j));
	runge.Write("runge.txt");
	return 0;
}

void run_tests(int n, int m, Mat& a, Mat& b)
{
	//Matrix creation
	Mat x = Linspace(-5, 5, m+1);
	Mat y = Linspace(-5, 5, n+1);
	Mat z(m+1, n+1);
	for(int i = 0; i < m+1; i++)
		for(int j = 0; j < n+1; j++)
			z(i,j) = f(x(i), y(j));

	Mat p(201, 201);
	for(int i = 0; i < 201; i++)
		for(int j = 0; j < 201; j++)
			p(i,j) = lagrange2D(x, y, z, a(i), b(j));
	string file_out = "p";
	char* buffer = new char[5];
	sprintf(buffer, "%d", m);
	file_out += buffer;
	file_out += "_reg.txt";
	p.Write(file_out.c_str());
}

double f(const double x, const double y)
{
	return 1 / (1 + x*x + y*y);
}