// =============================================================================================
// ***** Header for random.cpp: set of routines for random sampling and related functions *****
// =============================================================================================
#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"

extern "C"{

void
rltruncGamma(double* x,  
             const double* shape,  const double* rate,   const double* minx,
             const int* n,         const int* callFromR);

void
discreteUniformSampler(int* sampledj,
                       const int* kP,  const int* nP,  const int* callFromR);

void
discreteSampler(int* sampledj,  double* propA,
                const int* kP,  const int* nP,  const int* cumul,  const int* callFromR);

void
discreteSampler2(int* sampledj,  double* propA,
                 const int* kP,  const int* nP,  const int* cumul,  const int* callFromR);

int 
findUniformIndex(const double u,  const int startInd,  const int endInd,  const int k);

int 
findIndex(const double u,  const int startInd,  const int endInd,  const double* ValuesA);

}

#endif
