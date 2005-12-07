// ============================================================================================
// ***** Header for marginal_bayesGspline.cpp:                                         ***** //
// ============================================================================================
#ifndef _MARGINAL_BAYES_G_SPLINE_H_
#define _MARGINAL_BAYES_G_SPLINE_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "in_output_GS.h"

extern "C"{

void
marginal_bayesGspline
   (double* average1,   double* average2,      double* value1,       double* value2,
    int* M_now,         const int* onlyAver,   char** dirP,          char** extensP,
    const double* x1,   const double* x2,
    const int* KK,      
    const int* M,       const int* skip,       const int* by,        const int* nwrite,
    const int* nx1,     const int* nx2,        int* errP);
}

void
marginal_evalGspline
   (double* average1,     double* average2,     double* value1,        double* value2,  
    const int* length0,   const int* length1,
    const int* nx1,       const int* nx2,       const double* x1,      const double* x2,
    double** w,           double** mu,          const double* intcpt,  
    const double* sigma,  const double* scale,  double* inv_sigmas,    double* min_half_inv_sig2);

#endif
