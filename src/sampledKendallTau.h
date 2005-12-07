// ============================================================================================
// ***** Header for bayesGspline.cpp:                                                  ***** //
// ============================================================================================
#ifndef _SAMPLED_KENDALL_TAU_H_
#define _SAMPLED_KENDALL_TAU_H_

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
sampledKendallTau(double* Tau,           int* M_now,
                  char** dirP,           char** extensP,
                  const int* KK,
                  const double* Phi0,    const double* Phi1,   
                  const int* M,          const int* skip,     const int* by,         const int* nwrite,
                  int* errP);

}

void
evalKendallTau(double* value,  const int* dim,  const int* k_effect,  const double* w,  int** ind_mu,  double**** PhiPhi);

#endif
