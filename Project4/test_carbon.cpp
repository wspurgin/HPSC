/*
Will Spurgin
11/24/2013
High Performance Scietific Computing
MATH 3316
*/

#include "mat.h"
#include <string>
#include <iostream>

using namespace std;

//prototypes
double carbon(const double x, const double t,const double rtol,
    const double atol);

int main(int argc, char** argv)
{
    Mat x_ = Linspace(0, 4e-3, 400);
    x_.Write("x.txt");
    double* x = x_.get_data();

    Mat t_ = Linspace(0, 172800, 800);
    t_.Write("t.txt");
    double* t = t_.get_data();

    double rtol = 1e-10;
    double atol = 1e-15;

    Mat C(400, 800);
    for(int i = 0; i < 400; i++)
    {
        for(int j = 0; j < 800; j++)
        {
            C(i, j) = carbon(x[i], t[j], rtol, atol);
        }
    }
    C.Write("C.txt");
    cout << "Wrote initial output for C(x, t). Moving to hourly calculations"
         << endl;

    Mat C_hours(400, 1);
    double seconds[] = { 3600, 21600, 43200, 86400, 172800 };
    string output_files[] = { "C_1hour.txt", "C_6hour.txt", "C_12hour.txt",
        "C_24hour.txt", "C_48hour.txt" };
    for(int i = 0; i < 5; i++)
    {
        for(int j = 0; j < 400; j++)
        {
            C_hours(j) = carbon(x[j], seconds[i], rtol, atol);
        }
        C_hours.Write(output_files[i].c_str());
    }

    Mat C_3(800, 1);
    for(int i = 0; i < 800; i++)
        C_3(i) = carbon(3e-3, t[i], rtol, atol);
    C_3.Write("C_3mm.txt");

    cout << "Done." << endl;

    return 0;
}
