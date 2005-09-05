// Header for birthDeath.cpp
#ifndef _BIRTH_DEATH_H_
#define _BIRTH_DEATH_H_

#include <R.h>
#include <Rmath.h>

#include "random.h"
#include "templatefun.h"
#include "logLikelihood.h"
#include "randomLogLikelihood.h"
#include "List.h"

void
birthDeath(int* acceptedP,             int* birthP,                   int* kP, 
           double* loglikelhood,       double** loglikobs,            double** proploglikobs,
           double* randomloglik,       double** randomllcl,           double** proprandomllcl,
           double* wM,                 double* muM,                   double* invsigma2M,
           double* mixMomentM,
           int* rM,                    List<int>* invrM,              int* mixtureNM, 
           int* propkP,                          
           double *uM,
           double (*logdu) (const double*, const double*),
           void (*transu) (double*, const double*, const double*),
           void (*invtransu) (double*, const double*, const double*),
           double (*logJtransu) (const double*, const double*, const double*),
           const double* regresResM,   const double* YsM,
           const double* bM,           const double* betaM,           const covMatrix* Dcm,
           const int* kmaxP,              
           const double* piBirthM,     const double* logpiBirthM,     const double* logpiDeathM,
           const double* deltaP,       const double* xiP,             const double* invkappaP, 
           const double* sqrtkappaP,   const double* halfl2pikappaP,  const double* zetaP,
           const double* etaP,         const double* lgammazetaP,     const double* llambdaP,
           const double* priorParmu,   double* transParmu,          
           const int* priorForkP,      const int* Eb0dependMix,       const int* randomIntP,
           const int* indbinXA,        const int* nP,                 const int* nclusterP,
           double (*logdtrans) (double));

void
proposeDeath(int& jdeath,        double* vM,
             const int& numE,    const int* emptyCompM,
             const double* wM,   const double* muM,        const double* invsigma2M);

double
logdtransBirthDeath(const double* vM,          const double* uM,  
                    const double* priorParmu,  const double* transParmu,          
                    double (*logdu) (const double*, const double*),
                    double (*logJtransu) (const double*, const double*, const double*),
                    const bool* RichardsonGreenP);

double
logPostRatioJacobianBirthDeath(const int* shortkP,       const double* vM,
                               const int* nP,
                               const double* deltaP,
                               const double* xiP,        const double* invkappaP,  const double* halfl2pikappaP,
                               const double* zetaP,      const double* etaP,       const double* lgammazetaP,
                               const double* llambdaP,   const int* priorForkP,    const bool* RichardsonGreenP);

void
moveParamsBirthDeath(int& jdeath,                   
                     double* wM,          double* muM,        double* invsigma2M,
                     int* rM,             List<int>* invrM,   int* mixtureNM,
                     const int* propkP,   const double* vM,   const int* birthP);

int
numEmpty(int* emptyCompM, const int *kP,  const int* mixtureNM);

#endif
