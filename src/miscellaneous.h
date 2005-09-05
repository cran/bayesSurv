// ===========================================================
// ***** Header for miscellaneous.cpp: Helping functions *****
// ===========================================================
#ifndef _MISCELLANEOUS_H_
#define _MISCELLANEOUS_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "List.h"

void
giveMixtureN(int* mixtureNM, 
             const int* kP,   const int* rM,  const int* nP);

void
giveMixtureN(int* mixtureNM, 
             const int* kP,  const List<int>* invrM);

void
Y2T(double** T,  const double* Y,  const int iter,  const int* nP,  double (*itrans)(const double));

void
giveSigmaAndInvsigma2(double* sigmaM,  double* invsigma2M,  const double* sigma2M,  const int* kP);

void
regresResidual(double* regresResA,
               const double* YsA,    const double* betaA,     const double* bA,
               const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP);

void 
regresPredictor(double* regresPredA,
                const double* betaA,  const double* bA,
                const double* XA,     const int* clusteriA,    const int* randomIntP,   
                const int* indbA,     const int* nP,           const int* nXP,         const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* YsA,   const double* newYsA,
               const int* nP);

void
regresResidual(double* regresResA,
               const double* betaA,   const double* newbetaA,   const int* indnewA,  
               const int* nnewP,      const double* XA,         const int* indbA,         
               const int* nP);

void
regresResidual(double* regresResA,
               const double* bA,      const double* newbA,  const int* indnewA,
               const int* nnewP,      const double* XA,     const int* clusteriA,    const int* randomIntP,   
               const int* indbinXA,   const int* nP,        const int* nXP,          const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP);

void
regresResidual(double* regresResA,
               const double* bA,         const double* bclA,   
               const int* indnewA,       const int* nnewP,     const int* cl,
               const List<int>* indobs,  const double* XA,     const int* randomIntP,
               const int* indbinXA,      const int* nP,        const int* nXP,         const int* nrandomP);

#endif
