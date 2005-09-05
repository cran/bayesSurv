#ifndef _UPDATE_RANDOM_H_
#define _UPDATE_RANDOM_H_

#include <R.h>
#include <Rmath.h>

#include "bblocks.h"
#include "covMatrix.h"
#include "List.h"

#include "mvtdist.h"
#include "MHproposal.h"
#include "miscellaneous.h"
#include "logLikelihood.h"
#include "randomLogLikelihood.h"


void
updateRandom(bblocks* bb,                   double* regresResM,    double* propregresResM,
            double* randomloglik,           double* randomllcl,    double* proprandomllcl,
            double* loglikelhood,           double* loglikobs,     double* proploglikobs,
            const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
            const double* YsM,              const double* XA,      double** ZZt,
            const List<int>* invclusteriA,  const int* indbinXA,   double (*logdtrans)(double),
            const int* randomIntP,          const int* nXP,        const int* nP,
            const int* errorTypeP,          const int* kP,         const int* rM,
            const double* wM,               const double* muM,     const double* invsigma2M,
            const double* tolerChol);

void
GIBBSproposalRandom(double* bM,                     double* regresResM, 
                    double* covpar,                 double* ichicovpar,
                    const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
                    const double* XA,               double** ZZt,        
                    const double* wM,               const double* muM,     const double* invsigma2M,  const int* rM,
                    const int* randomIntP,          const int* nrandomP,       
                    const int* nXP,                 const int* nclusterP,  const int* nP,             const int* errorTypeP,
                    const List<int>* invclusteriA,  const int* indbinXA,   const int* indUpd,
                    const double* tolerChol);

#endif

