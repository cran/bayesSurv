// Header for randomLogLikelihood.cpp
#ifndef _RANDOM_LOG_LIKELIHOOD_H_
#define _RANDOM_LOG_LIKELIHOOD_H_

#include <R.h>
#include <Rmath.h>

#include "covMatrix.h"
#include "mvtdist.h"

void
randomLogLikelihood(double* randomloglik,    double* randomllcl,
                    const int* clusterUpd,   const int* nclusterUpd,  const int* nclusterP,
                    const double* bM,        const double* betaM,     const covMatrix* Dcm,
                    const double* Eb0,       const int* indbinXA);

void
randomLogLikelihood(double* randomloglik,   double* randomllcl,
                    const int* clusterUpd,  const int* nclusterP, 
                    const double* bclM,     const double* betaM,    const covMatrix* Dcm,
                    const double* Eb0,      const int* indbinXA);

#endif
