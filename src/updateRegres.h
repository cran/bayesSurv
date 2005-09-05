#ifndef _UPDATE_REGRES_H_
#define _UPDATE_REGRES_H_

#include <R.h>
#include <Rmath.h>

#include "covMatrix.h"
#include "MHblocks.h"

#include "miscellaneous.h"
#include "MHproposal.h"
#include "logLikelihood.h"
#include "randomLogLikelihood.h"

void
updateRegres(MHblocks* betaMH,             double** regresResM,   double** propregresResM,   
             double* loglikelhood,         double** loglikobs,    double** proploglikobs,
             double* randomloglik,         double** randomllcl,   double** proprandomllcl,
             const double* YsM,            const double* Eb0,     const double* bM,      
             const covMatrix* Dcm,         const double* XA,      double*** XXt, 
             const int* indbA,             const int* indbinXA,
             double (*logdtrans)(double),  const int* iterTotal,  const int* doadapt,
             const int* randomIntP,        const int* nrandomP,
             const int* nP,                const int* nclusterP,
             const int* errorTypeP,        const int* kP,         const int* rM,
             const double* wM,             const double* muM,     const double* invsigma2M,
             const double* tolerChol);

void
GIBBSproposalFixed(double* betaM,            double* regresResM,
                   double* covpar,           double* ichicovpar,
                   const double* priormean,  const double* priorinvvar,
                   const double* Eb0,        const int* randomIntP,
                   const double* XA,         double** XXtb,              const int* diagIXXtb,
                   const double* wM,         const double* muM,          const double* invsigma2M,  const int* rM,
                   const int* indUpd,        const int* npar,            const int* nupdate,        const int* nP,
                   const int* errorTypeP,    const double* tolerChol);

void
GIBBSproposalMeanRandom(double* betaM,
                        double* covpar,           double* ichicovpar,
                        const double* priormean,  const double* priorinvvar,
                        const double* Eb0,        const double* bM,
                        const covMatrix* Dcm,     const int* diagIXXtb,
                        const int* indbA,         const int* indbinXA,
                        const int* indUpd,        const int* npar,            
                        const int* nupdate,       const int* nrandomP,         const int* nclusterP,
                        const int* errorTypeP,    const double* tolerChol);

#endif


