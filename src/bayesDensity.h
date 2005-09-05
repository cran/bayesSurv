/***** Header for bayesDensity.cpp:                                        *****/
/*** Estimation of the density based on MCMC sample from bayessurvreg1     *****/
/* =========================================================================== */
#ifndef _BAYES_DENSITY_H_
#define _BAYES_DENSITY_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "in_output.h"

extern "C"{

void
bayesDensity(double* aver,        double* staver,        double* centaver,
             double* intercept,   double* scale,         int* Mk,
             char** dirP,       
             const double* grid,  const double* stgrid,  const double* centgrid,
             const int* kmax,     const int* M,          const int* skip,         const int* by,
             const int* ngrid,    const int* nstgrid,    const int* ncentgrid,
             int* errP);

}

#endif
