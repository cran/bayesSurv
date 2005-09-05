#ifndef _UPDATE_MEANS_H_
#define _UPDATE_MEANS_H_

#include <R.h>
#include <Rmath.h>

#include "covMatrix.h"
#include "List.h"
#include "mixMoments.h"

void
updateMeans(double* muM,                double* mixMomentM,
            const double* regresResM,   const double* betaM,        const double* bM,
            const covMatrix* Dcm,
            const int* kP,              const int* mixtureNM,
            const double* wM,           const double* invsigma2M,   const List<int>* invrM,  
            const double* xiInvkappaP,  const double* invkappaP,  
            const int* Eb0dependMix,    const int* randomIntP,      
            const int* nP,              const int* nclusterP,       const int* nrandomP,
            const int* indbinXA);

#endif
