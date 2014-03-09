/*
Will Spurgin
11/24/2013
High Performance Scientific Computing
MATH 3316
*/

#include <iostream>
#include <cmath>

using namespace std;

//prototypes
int adaptive_int(double (*f)(const double), const double a, const double b,
    const double rtol, const double atol, double &R, int &n);

const double c=3.0;
const double d=20.0;
double f(const double x) {
  return (exp(c*x) + sin(d*x));
}

int main(int argc, char** argv)
{
    double atol = 1e-15;
    double rtol[] = { 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12 };
    double a = -5.0;
    double b = 2.0;

    cout << "rtol               R              n        relerr" << endl
        << "--------------------------------------------------" << endl;
        double Itrue = 1.0/c*(exp(c*b) - exp(c*a)) - 1.0/d*(cos(d*b)-cos(d*a));
        for(int i = 0; i < 6; i++)
        {
            double R = 0.0;
            int n = 0;
            if(adaptive_int(f, a, b, rtol[i], atol, R, n) == 0)
            {
                if(i != 0)
                {
                    double relerr = fabs(Itrue-R)/fabs(Itrue);
                    printf("%7.1e  %22.16e %5i     %7.1e\n", rtol[i], R, n, relerr);
                }
                else
                    printf("%7.1e  %22.16e %5i\n", rtol[i], R, n);
            }
            else
            {
               printf("%7.1e  %22.16e  %5i\n", rtol[i], R, n);
               cout << " (Did not converge) " << endl;
            }
        }
    return 0;
}
