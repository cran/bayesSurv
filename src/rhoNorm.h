/*** rhoNorm.h ***/

#ifndef _RHO_NORM_H_
#define _RHO_NORM_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"
#include "AK_BLAS_LAPACK.h"

namespace rhoNorm {

enum rhoNormAlgorithm {_Normal_Around_Mode_, _Langevin_};

const double rhoONE    = 1 - 1e-15;
const double rhoMinONE = -1 + 1e-15;
const double zeps      = 17.61636;        /** = z for rho equal to 1 - 1e-15 **/

const int    _maxiter      = 1;          /** Maximal number of Newton-Raphson iterations when constructing the normal approximation   **/
const double _toler        = 1e-3;       /** Tolerance to detect convergence of the Newton-Raphson                                    **/
const int    _max_stephalf = 10;         /** Maximal number of stephalving steps within Newton-Raphson                                **/


void
z2rho(double *rho,  const double *z);

void
rho2z(double *z,  const double *rho);

void
rho2zError(double *z,  const double *rho);

void
update_pUnif(int *accept,           double *z,             double *rho,          double *work,
             const double *sumu2,   const double *sumv2,   const double *sumuv,  const int *nobs,
             const int *algorithm,  const double *scaleL);

void
ML_est(double *ll,             double *dll,           double *ddll,
       double *z,              double *rho,           int *niter,           int *err,
       const double *sumu2,    const double *sumv2,   const double *sumuv,  const int *nobs,
       const int *maxiter);

void
lposter0(double *ll,           double *rho,          const double *z,  
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs);

void
lposter1(double *ll,           double *dll,          double *rho,          const double *z,  
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs);

void
lposter2(double *ll,           double *dll,          double *ddll,         double *rho,          const double *z,  
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs);

extern "C"{
  void
  mcmc_rhoNorm(int *acceptSample,     double *zSample,      double *rhoSample,     int *iter,
               const double *sumu2,   const double *sumv2,  const double *sumuv,   const int *nobs,  const int *nsimul,
               const int *algorithm,  const double *scaleL);
}

}

#endif


