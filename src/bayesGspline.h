// ============================================================================================
// ***** Header for bayesGspline.cpp:                                                  ***** //
// ============================================================================================
#ifndef _BAYES_G_SPLINE_CPP_
#define _BAYES_G_SPLINE_CPP_

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
bayesGspline(double* average,          double* value,   int* M_now,          const int* onlyAver,
             char** dirP,              char** extensP,  char** extens_adjP,  
             double* x1,               double* x2,   
             const int* total_length,  const int* M,    const int* skip,     const int* by,         const int* nwrite,
             const int* nx1,           int* nx2,        const int* version,  int* standard,
             int* errP);

}

void
evalGspline(double* average,      double* value,  
            const int* nx1,       const int* nx2,       double** x,
            const int* dim,       const int* k_effect,  
            const double* w,      double** mu,          const double* intcpt,  
            const double* sigma,  const double* scale,  double* min_half_inv_sig2,
            const int* standard,  const double* E_gx,   const double* sd_gx);

void
printReadGspline(const int& iter,       const int& dim,       const int& k_effect,  const double* w,  double** mu,  
                 const double* intcpt,  const double* sigma,  const double* scale);

#endif
