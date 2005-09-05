#ifndef _UPDATE_COV_MAT_RANDOM_H_
#define _UPDATE_COV_MAT_RANDOM_H_

#include <R.h>
#include <Rmath.h>

#include "covMatrix.h"
#include "random.h"
#include "mvtdist.h"

void
updateCovMatRandom(covMatrix* Dcm,           const double* bM,             
                   const double* betaM,      const double* Eb0,
                   const int* priorD,        const double* priordf,   const double* priorscaleMat, 
                   const int* indbinXA,      const int* nclusterP,    const int* nP,      
                   const double* tolerChol,  const double* tolerQR);

void
sumSquares(double* sumSq,          const double* bM,       
           const double* betaM,    const double* Eb0,
           const int* indbinXA,    const int* diagI,       
           const int* nclusterP,   const int* nrandomP,    const int* lSS);

#endif
