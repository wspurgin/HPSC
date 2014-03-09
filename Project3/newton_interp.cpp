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

using namespace std;

//This version uses Divided Differences to calculate the coefficients
int newton_coeffs(Mat& x, Mat& y, Mat& c)
{
	//Error checking against unequal lengths
	if(x.Cols()*x.Rows() != y.Cols()*y.Rows())
	{
		cerr << "x and y have unequal lengths" << endl;
		return 1;
	}
	//Get the data as an array to speed up random access
	double* xd = x.get_data();
	double* yd = y.get_data();

	//first coefficient will be equal to the first y
	c(0) = yd[0];
	for(int i = 1; i < x.Cols()*x.Rows(); i++)
	{
		//initialize the coefficient to 0
		c(i) = 0;
		for(int j = 0; j <= i; j++)
		{
			/*
			the 'divisor' is made up of the product of all
			x[j] - x[k]s where k != j from 0 <= k <= i
			*/
			double divisor = 1.0;
			for(int k = 0; k <= i; k++)
				if(k != j) divisor *= (xd[j] - xd[k]);
			c(i) += yd[j] / divisor; 
		}
	}
	return 0;
}

double newton_eval(Mat& x, Mat& c, double z)
{
	//Place data in arrays
	double* xd = x.get_data();
	double* cd = c.get_data();

	//initialize sum to zero
	double sum = 0.0;
	//evaluate at z
	for(int i = 0; i < x.Cols()*x.Rows(); i++)
	{
		double product = 1.0;
		for(int j = 0; j < i; j++)
			product *= (z - xd[j]);
		sum += cd[i]*product;
	}
	return sum;
}
