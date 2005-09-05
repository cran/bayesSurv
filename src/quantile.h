// =======================================================
// ***** Header for quantile.cpp: *****
// =======================================================
#ifndef _QUANTILE_H_
#define _QUANTILE_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

/***** Functions used in predictive.cpp *******/
void
cumsumQuantile1(double** quant,  double** newval,  const int nquant,  const int nobs,  
                const int iter);

void
cumsumQuantile2(double*** quant,  double*** newval,     const int nquant,  const int nobs,  const int* ngridM,  
                const int iter);

void
meanQuantile1(double** quant,           double** sample,   
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int nquant,      const int sampleSize);

void
meanQuantile2(double*** quant,          double*** sample,  
              const double* quantileA,  const int* indquant1,  const int* indquant2,
              const int nobs,           const int* ngridM,     const int nquant,      const int sampleSize);

/***** Functions used in predictive_GS.cpp *******/
void
resetAverage(double* average,  const int* nobs,  const int* ngrid,  const int* predict);

void
cumsum2average(double* average,  const int* sampleSize,  const int* nobs,  const int* ngrid,  const int* predict);

void
value2quantile(double* value,          double* quantile,  
               const double* probs,    const int* indquant1,  const int* indquant2,  const int* nquant,
               const int* sampleSize,  const int* nobs,       const int* ngrid,      const int* predict,    const int* shift_pointer);

#endif
