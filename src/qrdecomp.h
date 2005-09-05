// ======================================================================================
// ***** Header for qrdecomp.cpp: Routines for a QR decomposition of a matrix
//                    - grabbed from dqrdc2.f and blas.f
// ======================================================================================
#ifndef _QR_DECOMP_H_
#define _QR_DECOMP_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"

extern "C"{

void 
dqrdc2CPP(double* x, const int* n, const int* p, const double* tol, int* k, double* qraux, int* jpvt);

double
ddotCPP(const int n, double* dx, const int incx, double* dy, const int incy);

void
daxpyCPP(const int n, const double da, double* dx, const int incx, double* dy, const int incy);

void
dscalCPP(const int n, const double da, double* dx, const int incx);

double 
dnrm2CPP(const int n, const double* x, const int incx);

}

#endif
