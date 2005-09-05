#ifndef _NEWTON_RAPHSON_H_
#define _NEWTON_RAPHSON_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

#ifdef __cplusplus
extern "C" {
#endif

void
newton_raphson(double* x,            double* gx,          double* dgx,  double* ddgx,  
               const double* parmD,  const int* parmI,
	       void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
               int* iter,            const int* maxiter,  const int* max_stephalf,  
               const double* toler,  const double* zero,  int* err);

void
solver_newton_raphson(double* x,            double* gx,          double* dgx,  const double* b,
                      const double* parmD,  const int* parmI,
	              void (*eval2)(const double*,  double*, double*, double*, const double*, const int*, const int&),
                      int* iter,            const int* maxiter,
                      const double* toler,  const double* zero,  int* err);

#ifdef __cplusplus
}
#endif

#endif
