/*
Will Spurgin
10/26/2013
High Performance Scientific Computing
MATH 3316
Project 3 - Interpolation
*/

#include <iostream>
#include <math.h>
#include "mat.h"

double lagrange_basis(Mat &x, int i, double z);

double lagrange2D(Mat& x, Mat& y, Mat& z, double a, double b)
{
	double p_value = 0.0;
	for(int i = 0; i < x.Cols()*x.Rows(); i++)
	{
		for(int j = 0; j < y.Cols()*y.Rows(); j++)
			p_value += z(i, j) * lagrange_basis(x, i, a) * lagrange_basis(y, j, b);
	}
	return p_value;
}
