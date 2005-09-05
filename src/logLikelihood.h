// Header for logLikelihood.cpp
#ifndef _LOG_LIKELIHOOD_H_
#define _LOG_LIKELIHOOD_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "List.h"

void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double));

void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const List<int>* obsUpd,    const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double));

void
clusterlogLikelihood(double* clusterloglik,  const double* loglikobs,  const int* cl,  const List<int>* obsInCluster);

#endif

