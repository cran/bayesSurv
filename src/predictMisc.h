// =======================================================
// ***** Header for predictSurv.cpp: *****
// =======================================================
#ifndef _PREDICT_MISC_H_
#define _PREDICT_MISC_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "covMatrix.h"
#include "bblocks.h"
#include "random.h"
#include "mvtdist.h"

void
predictSurv(double*** SM,              double*** lambdaM,         double*** LambdaM,
            const int iter,            double** gridM,            double** loggridM,  const double* time0P,
            const double* regresPredM,
            const int* rM,             const double* wM,          const double* muM,  const double* sigmaM,
            const double* Eb0,         const int* kP,             const int* nP,      const int* ngridM,
            const int* errorTypeP,     const int* randomIntP,
            const int* hazardP,        const int* cumhazardP);

void
predictData(double* YsM,             const double* regresPredM,
            int* rM,                 double* cumwM,             const double* muM,     const double* sigmaM,
            const double* Eb0,       const int* kP,             const int* nP,        
            const int* errorTypeP,   const int* randomIntP);

void
predictRandom(double* bM,  
              const double* betaM,     const double* Eb0,    const covMatrix* Dcm,
              const int* nrandomP,     const int* nclusterP,  
              const int* indbinXA,     const int* indUpd);

void
predictET(double** ET,            const double* time0P, const int iter,
          const double* betaM,    const double* wM,     const double*muM,       const double* sigma2M,
          const covMatrix* Dcm,   const double* XA,
          const int* kP,          const int* nP,        const int* nXP,         const int* indbinXA,
          const int* randomIntP,  const int* nrandomP,  const int* errorTypeP);

void
Y2T(double** T,  const double* Y,  const double* time0P, const int iter,  const int* nP,  double (*itrans)(const double));

#endif
