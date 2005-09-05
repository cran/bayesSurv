// ================================================================================================
// ***** Header for cholesky.cpp: Routines for a Cholesky decomposition and related functions *****
// ================================================================================================
#ifndef _CHOLESKY_H_
#define _CHOLESKY_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"

extern "C"{

void
cholesky(double* C,  int* rankC,  const int* nC,  const int* diagI, const double* toler);

void
cholesky2(double* C,  int* rankC,  const int* nC,  const double* toler);

void 
chinv(double* C , const int* nC,  const int* diagI, const int* onlyCholInv);

void 
chinv2(double* C , double* ichol, const int* nC,  const int* diagI);

void
chposDef(const double* C,      double* cholC,      int* rankC,           
         int* attempt,         const int* nC,      const int* diagI,  
         const double* toler,  const double* eps,  const int* nattempt = &AKINT_MAX);

}

#endif
