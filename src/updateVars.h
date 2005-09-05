#ifndef _UPDATE_VARS_H_
#define _UPDATE_VARS_H_

#include <R.h>
#include <Rmath.h>

#include "constants.h"
#include "mixMoments.h"

void
updateVars(double* invsigma2M,        double* mixMomentM,      double* Eb0,
           const double* regresResM,  
           const int* kP,             const int* mixtureNM,
           const double* wM,          const double* muM,       const int* rM, 
           const double* zetaP,       const double* etaP,  
           const int* randomIntP,     const int* nP);

#endif
