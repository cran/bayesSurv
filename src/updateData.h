// =======================================================
// ***** Header for updateData.cpp: *****
// =======================================================
#ifndef _UPDATE_DATA_H_
#define _UPDATE_DATA_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"

void
updateData(double* YsM,             double* regresResM,
           const double* Y1M,       const double* Y2M,     const int* statusM,
           const int* rM,           double* cumwM,         const double* muM,     const double* invsigma2M,
           const double* Eb0,       const int* kP,         const int* nP,        
           const int* errorTypeP,   const int* randomIntP);


#endif
