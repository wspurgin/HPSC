/*
*	Will Spurgin
*	HPSC Project MATH 3316
*	Project 1, Taylor Series and Floating-Point Error
*/

#include <iostream>
#include <cmath>
#include "mat.h"

using namespace std;

int main()
{
	// Initialize vectors
	Mat n = Mat(1, 52);
	Mat q = Mat(1, 52);
	Mat q_p = Mat(1, 52);
	// set the constants c1 and c2
	double c1 = (1.0 / pow(3.0, 2.0)) / (2.0 / 3.0);
	double c2 = 9.0 * 1.1 * pow(10.0, -16.0);
	c1 = abs(c1);
	c2 = abs(c2);

	// Calculations for actual and potential errors
	for(int i = 0; i < 52; i++)
	{
		n(i) = (double)i + 1;
		double h = pow(2.0, -n(i));
		double actual_error = abs((log(3 + h) - log(3)) / h);
		actual_error = (actual_error - 1.0/3.0)*3;
		q(i) = -log10(abs(2.0*actual_error));
		double predicted_error = c1 * h + c2 / h;
		q_p(i) = -log10(abs(2.0*predicted_error));	
	}
	//Write vectors to the disk
	n.Write("n.txt");
	q.Write("q.txt");
	q_p.Write("q_p.txt");
	return 0;
}