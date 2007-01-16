// Header for predictive_GS.cpp
#ifndef _PREDICTIVE_GS_H_
#define _PREDICTIVE_GS_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "quantile.h"
#include "classBetaGamma.h"
#include "classRandomEff.h"
#include "classCovMatrix.h"
#include "in_output_GS.h"
#include "regresResidual_GS.h"
#include "structRandomEff32.h"

const double _stres0      = 20.0;                   // dnorm(_stres0, 0, 1, 0) = 5.52e-88, which is zero for me 
const double _zero_weight = 1e-6;                   // limit for weights to ignore mixture components in the marginal distributions

extern "C"{

void
predictive_GS(double *averDens,         double *averS,           double *averHaz,      double *averCumHaz,
              double *valDens,          double *valS,            double *valHaz,       double *valCumHaz,
              double *quantDens,        double *quantS,          double *quantHaz,     double *quantCumHaz,
              const int *dimsP,         const double *X,         const int *obsdims,
              int *M_now,               char **dirP,             char **extensP,       char **extens_adjP,
              const int *GsplI,
              const int *objBetaI,      const double *objBetaD,
              const int *objbI,         const double *objbD,
              const int *b_GsplI,
              const double *gridA,      const double *loggridA,  const int *ngrid,
              double *probsA,           const int *nquant,       int *onlyAver,        const int *predictP,
              const int *M,             const int *skip,         const int *by,        const int *nwrite,
              const int *version,       const int *Onset,        int *errP);

}

void
evalPredFuns(double* averDens,        double* averS,               double* averHaz,            double* averCumHaz,
             double* valDens,         double* valS,                double* valHaz,             double* valCumHaz,
             const int* obsdims,      const int* nobs,             const int* ngrid,           
             const double* gridA,     const double* loggridA,
             const double* linPred,   const int* dim,              const int* Glength,
             double** const w_marg,   double** const mu_sig_marg,  const double* intcpt,
             const double* sigma,     const double* scale,         double* inv_sig_scale,
             const int* predictP,     const double* zero_weight,   const int* shift_pointer);


#endif

