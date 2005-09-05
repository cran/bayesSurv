#ifndef _UPDATE_WEIGHTS_H_
#define _UPDATE_WEIGHTS_H_

#include "mixMoments.h"
#include "templatefun.h"
#include "logLikelihood.h"
#include "randomLogLikelihood.h"

void
updateWeights(double** wM,               double** propwM,            double* mixMomentM,
              double* loglikelhood,      double** loglikobs,         double** proploglikobs,
              double* randomloglik,      double** randomllcl,        double** proprandomllcl,
              const double* regresResM,  const double* YsM,          const double* betaM,        
              const double* bM,          const covMatrix* Dcm,
              const int* kP,             const int* mixtureNM,  
              const double* muM,         const double* invsigma2M,   const int* rM,
              const double* deltaP,
              const int* Eb0dependMix,   const int* randomIntP,      
              const int* nP,             const int* nclusterP,       const int* nrandomP,
              const int* indbinXA,       double (*logdtrans)(double));

#endif
