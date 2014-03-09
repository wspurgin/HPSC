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

double russian_nodes(int i, int m, int a, int b);

void run_tests(int n, int m, Mat& a, Mat& b);

double lagrange2D(Mat &x, Mat &y, Mat &z, double a, double b);

int main(int argc, char** arv)
{
	Mat a(201, 1);
	Mat b(201, 1);

	for(int i = 0; i < 201; i++)
		a(i) = b(i) = russian_nodes(i, 200, -5, 5);

	//write a and b to disk
	a.Write("avals_cheb.txt");
	b.Write("bvals_cheb.txt");
	run_tests(10, 10, a, b);
	run_tests(20, 20, a, b);

	return 0;
}

void run_tests(int n, int m, Mat& a, Mat& b)
{
	//Matrix creation
	Mat x(m+1, 1);
	Mat y(n+1, 1);

	//Fill the matrices with chebyshev nodes
	for(int i = 0; i < m+1; i++)
		x(i) = russian_nodes(i, m, -5, 5);

	for(int i = 0; i < m+1; i++)
		y(i) = russian_nodes(i, n, -5, 5);

	//Fill the matrix 'z' with the proper function values
	Mat z(m+1, n+1);
	for(int i = 0; i < m+1; i++)
		for(int j = 0; j < n+1; j++)
			z(i,j) = f(x(i), y(j));

	//Evaluate the lagrange2D interpolant
	Mat p(201, 201);
	for(int i = 0; i < 201; i++)
		for(int j = 0; j < 201; j++)
			p(i,j) = lagrange2D(x, y, z, a(i), b(j));

	//write the interpolant to the proper file
	string file_out = "p";
	char* buffer = new char[5];
	sprintf(buffer, "%d", m);
	file_out += buffer;
	file_out += "_cheb.txt";
	p.Write(file_out.c_str());
}

double f(const double x, const double y)
{
	return 1 / (1 + x*x + y*y);
}

//calculate the proper node for i, m and the interval a to b
double russian_nodes(int i, int m, int a, int b)
{
	double numerator = ((2*i+1)*PI);
	double denominator = (2*m+2);
	double phi = numerator/denominator;
	return (.5*((a+b) + (b-a))) * cos(phi);
}
