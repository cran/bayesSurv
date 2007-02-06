/*** Mvtdist3.h ***/

#ifndef _MULTI_VARIATE_DISTRIBUTIONS_THREE_H_
#define _MULTI_VARIATE_DISTRIBUTIONS_THREE_H_

#include <cmath>

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BLAS_LAPACK.h"

namespace Mvtdist3 {

void
rwishartEye3(double *W, double *work, const double *nu, const int *dim);

void
rwishart3(double *W,  double *work,  const double *nu,  double *invS,  const int *dim,  const int &must_decomp);

extern "C"{
  void
  rwishartR3(double *W,  double *work,  const double *nu,  double *invS,  const int *dim,  const int *nrandom);
}

void
rmvnorm2006(double *x,  const double *mu,  const double *L,  const int *nx);

void
rmvnormZero2006(double *x,  const double *L,  const int *nx);

void
rmvnormQ2006(double *x,  const double *mu,  const double *L,  const int *nx);

void
rmvnormQZero2006(double *x,  const double *L,  const int *nx);

void
rmvnormC2006(double *x,  double *b,  const double *L,  const int *nx);

extern "C"{
  void
  rmvnormR2006(double *x,  double *mub,  double *QS,  int *err,  const int *nx,  const int *nrandom,  const int *version);
}

}  /** end of the namespace Mvtdist3 **/

#endif
