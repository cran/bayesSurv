// Header for regresResidual_GS.cpp
//
#ifndef _REGRES_RESIDUAL_GS_H_
#define _REGRES_RESIDUAL_GS_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "classBetaGamma.h"

void
regresRes_GS(double* regResid,  const double* YsA,      const BetaGamma* bg,  const double* bA,  
             const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster);

void
regresRes_GS2006(double* regResid,  const double* YsA,      const BetaGamma* bg,  const double* bA,  
                 const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster);

void
linPred_GS(double* linPred,   const BetaGamma* bg,    const double* bA,  
           const double* XA,  const int* nincluster,  const int* nobs,      const int* ncluster);

#endif
