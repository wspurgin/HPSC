/*
* Will Spurgin
* HPSC MATH 3316
* Project 1, Taylor Series and Floating-Point Error
*/

#include <iostream>
#include <math.h>
#include <cmath>
#include "mat.h"

using namespace std;

double nest(Mat& a, double x);
Mat e_factorial(int size);

int main(int argc, char** argv)
{
	// Initialize the vectors and row vector
	Mat z = Linspace(-1.0, 1.0, 201);
	Mat a_8 = e_factorial(8);
	Mat a_12 = e_factorial(12);
	Mat p_8 = Mat(201, 1);
	Mat p_12 = Mat(201, 1);
	Mat etrue = Mat(201, 1);
	Mat err_8 = Mat(201, 1);
	Mat err_12 = Mat(201, 1);
	Mat a_33 = e_factorial(33);
	// Run the calculations for all the values in the seperate vectors
	for(int i = 0; i < 201; i++)
	{
		p_8(i, 0) = nest(a_8, z(0, i));
		p_12(i, 0) = nest(a_12, z(0, i));
		etrue(i, 0) = exp(z(0, i));
		err_8(i, 0) = abs(etrue(i, 0) - p_8(i, 0));
		err_12(i, 0) = abs(etrue(i, 0) - p_12(i, 0));

	}
	// write all the vectors to the disk.
	p_8.Write("p_8.txt");
	p_12.Write("p_12.txt");
	z.Write("z.txt");
	etrue.Write("etrue.txt");
	err_8.Write("err_8.txt");
	err_12.Write("err_12.txt");

	// Initialize the 'true' value of e^-13.2 and the approximation through the taylor polynomial.
	double etrue13 = exp(-13.2);
	double p_33 = nest(a_33, -13.2);
	cout << "e^-13.2 = " << etrue13 << endl << endl;
	// set up a nice little table of values
	cout << "Calculation: 		Value:		Relative Error:" << endl;
	cout << "p33(-13.2) 		" << p_33 << "			" << abs(etrue13 - p_33) / etrue13 << endl;
	p_33 = 1.0 / nest(a_33, 13.2);
	cout << "1/p33(13.2)		" << p_33 << "			" << abs(etrue13 - p_33) / etrue13 << endl;
	etrue13 = exp(13.2);
	p_33 = nest(a_33, 13.2);
	cout << "p33(13.2)		" << p_33 << "				" << abs(etrue13 - p_33) / etrue13 << endl; 

	return 0;
}

// nested multiplication function
double nest(Mat& a, double x)
{
	double p = a(0, a.Cols() - 1);
	for(int i = a.Cols() - 2; i >= 0; i--)
	{
		p = a(0, i) + p*x;
	}
	return p;
}

/*
This function is simple and made specifically for this model
and project. It has no error checking for size, which would be a usful
improvement to make this function better.
*/
Mat e_factorial(int size)
{
	// set up the data to be passed to the Mat constructor	
	double* data = new double[size];
	// The first entry for will be e^0 * x^0 = 1
	data[0] = 1.0;
	// Add the rest of the values for data through nested division
	for(int i = 1; i < size; i++)
		data[i] = (data[i - 1] / (double) i);
	Mat a = Mat(1, size, data);
	return a;
}
