// =======================================================
// ***** Header for mixMoments.cpp: *****
// =======================================================
#ifndef _MIX_MOMENTS_H_
#define _MIX_MOMENTS_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

void
mixMoments(double* mixMomentM,  const int* kP,             const double* wM,
           const double* muM,   const double* invsigma2M,  const bool onlySD = false);

void
mixMean(double* mixMeanP,  const int* kP,  const double* wM,  const double* muM);

#endif
