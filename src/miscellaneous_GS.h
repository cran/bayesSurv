// ========================================================================================================
// ***** Header for miscellaneous_GS.cpp: different helping functions for G-spline applications    *****//
// ========================================================================================================
#ifndef _MISCELLANEOUS_GS_H_
#define _MISCELLANEOUS_GS_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

extern "C"{
void
findClosestKnot(int* index,  const double* knot,  const double* obs,  const int* nknot,  const int* nobs);
}

#endif
