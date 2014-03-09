/* Daniel R. Reynolds
   SMU Mathematics
   Math 3316
   11 September 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
//#include <time.h>
#include <chrono>
#include "mat.h"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

// solver prototypes
int BackSubRow1(Mat &U, Mat &x, Mat &b);
int BackSubRow_unsafe(Mat &U, Mat &x, Mat &b);
int BackSubRow(Mat &U, Mat &x, Mat &b);
int BackSubCol_unsafe(Mat &U, Mat &x, Mat &b);
int BackSubCol(Mat &U, Mat &x, Mat &b);
int FwdSubRow(Mat &L, Mat &x, Mat &b);
int FwdSubCol(Mat &L, Mat &x, Mat &b);
int NaiveGaussElimRow1(Mat &A, Mat &b);
int NaiveGaussElimRow(Mat &A, Mat &b);
int NaiveGaussElimCol1(Mat &A, Mat &b);
int NaiveGaussElimCol(Mat &A, Mat &b);
int GaussElimRow_unsafe(Mat &A, Mat &b);
int GaussElimCol_unsafe(Mat &A, Mat &b);
int GaussElimRow(Mat &A, Mat &b);
int GaussElimCol(Mat &A, Mat &b);


// driver routine for testing linear solver algorithms
int main(int argc, char* argv[]) {

  // local variables
  long int sizes[] = {1500, 2000, 4000, 8000};
  int nsizes = 4;
  long int sizes2[] = {300, 400, 500, 600, 700, 800, 1000, 1200};
  int nsizes2 = 8;
  chrono::time_point<chrono::system_clock> stime, ftime;
  chrono::duration<double> runtime;

  // run triangular solver tests
  for (int k=0; k<nsizes; k++) {

    // set the problem size
    long int n = sizes[k];

    // display current problem size
    cout << "\nTesting triangular solvers with " << n << " x " << n << " matrices:\n";

    // allocate the matrices & vectors of this size
    Mat U = Random(n,n);
    Mat L = Random(n,n);
    Mat xtrue = Random(n,1);
    Mat Usave(n,n);
    Mat Lsave(n,n);
    Mat x(n,1);
    Mat b(n,1);
    Mat err(n,1);

    // fix matrices to be appropriate shape; add weight to diagonal
    for (int j=0; j<n; j++) {
      for (int i=j+1; i<n; i++)
        U(i,j) = 0.0;
      U(j,j) += 10.0;
    }
    for (int j=0; j<n; j++) {
      for (int i=0; i<j; i++)
	    L(i,j) = 0.0;
      L(j,j) += 10.0;
    }

    // store matrices for reuse
    Usave = U;
    Lsave = L;
    
    // create RHS from true solution, solve
    b = 0.0;
    b.Axpy(U, xtrue);   // b = b + U*xtrue
    stime = chrono::system_clock::now();
    BackSubRow1(U, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   BackSubRow1:       time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    U = Usave;

    // create RHS from true solution, solve
    b = 0.0;
    b.Axpy(U, xtrue);
    stime = chrono::system_clock::now();
    BackSubRow_unsafe(U, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   BackSubRow_unsafe: time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    U = Usave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(U, xtrue);
    stime = chrono::system_clock::now();
    BackSubRow(U, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   BackSubRow:        time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    U = Usave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(L, xtrue);
    stime = chrono::system_clock::now();
    FwdSubRow(L, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   FwdSubRow:         time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    L = Lsave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(U, xtrue);
    stime = chrono::system_clock::now();
    BackSubCol_unsafe(U, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   BackSubCol_unsafe: time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    U = Usave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(U, xtrue);
    stime = chrono::system_clock::now();
    BackSubCol(U, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   BackSubCol:        time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    U = Usave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(L, xtrue);
    stime = chrono::system_clock::now();
    FwdSubCol(L, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   FwdSubCol:         time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    L = Lsave;

  }


  // run square solver tests
  for (int k=0; k<nsizes2; k++) {

    // set the problem size
    long int n = sizes2[k];

    // display current problem size
    cout << "\nTesting square solvers with " << n << " x " << n << " matrices:\n";

    // allocate the matrices & vectors of this size
    Mat A = Random(n,n);
    Mat xtrue = Random(n,1);
    Mat Asave(n,n);
    Mat x(n,1);
    Mat b(n,1);
    Mat err(n,1);

    // fix matrices to be appropriate shape; add weight to diagonal
    for (int j=0; j<n; j++)
      A(j,j) += 10.0;

    // store matrices for reuse
    Asave = A;
    
    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    NaiveGaussElimRow1(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   NaiveGaussElimRow1:  time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    NaiveGaussElimRow(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   NaiveGaussElimRow:   time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    GaussElimRow_unsafe(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   GaussElimRow_unsafe: time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    GaussElimRow(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   GaussElimRow:        time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    NaiveGaussElimCol1(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   NaiveGaussElimCol1:  time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    NaiveGaussElimCol(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   NaiveGaussElimCol:   time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    GaussElimCol_unsafe(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   GaussElimCol_unsafe: time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

    // create RHS from true solution, solve 
    b = 0.0;
    b.Axpy(A, xtrue);
    stime = chrono::system_clock::now();
    GaussElimCol(A, b);
    BackSubCol(A, x, b);
    ftime = chrono::system_clock::now();
    runtime = ftime-stime;
    err = x-xtrue;
    printf("   GaussElimCol:        time = %.4f,  err = %.2g\n", 
	   runtime.count(), err.Norm());
    A = Asave;

  }

} // end main

