// Header for splitCombine.cpp
#ifndef _SPLIT_COMBINE_H_
#define _SPLIT_COMBINE_H_

#include <R.h>
#include <Rmath.h>

#include "constants.h"
#include "random.h"
#include "templatefun.h"
#include "logLikelihood.h"
#include "List.h"

void
splitCombine(int* acceptedP,                int* splitP,                 
             double* loglikelhood,          double** loglikobs,          double** proploglikobs,
             int* kP,                             
             double** wM,                   double** muM,                double** invsigma2M,
             int* rM,                       List<int>** invrM,           int** mixtureNM,
             int* propkP,                         
             double** propwM,               double** propmuM,            double** propinvsigma2M,
             int* proprM,                   List<int>** propinvrM,       int** propmixtureNM,
             double* uM,              
             double (*logdu) (const double*, const double*),
             void (*transu) (double*, const double*, const double*),
             void (*invtransu) (double*, const double*, const double*),
             double (*logJtransu) (const double*, const double*, const double*),
             const double* regresResM,      const double* YsM,            const double* Eb0,
             const int* kmaxP,              const int* randomIntP,
             const double* piSplitM,        const double* logpiSplitM,    const double* logpiCombineM,
             const double* deltaP,          const double* xiP,            const double* invkappaP,                     
             const double* halfl2pikappaP,  const double* zetaP,          const double* etaP,
             const double* lgammazetaP,     const double* llambdaP,       
             const double* priorParmu,      const double* transParmu,
             const int* priorForkP,         double (*logdtrans) (double), const int* nP);

void
proposeSplit(int* acceptedP,
             double* propwM,      double* propmuM,     double* propinvsigma2M,
             const double* wM,    const double* muM,   const double* invsigma2M,
             const double* vM,    const int jsplit,    const int* kP);

void
proposeCombine(int* acceptedP,
               double* vM,
               double* propwM,    double* propmuM,     double* propinvsigma2M,
               const double* wM,  const double* muM,   const double* invsigma2M,
               const int jsplit,  const int* propkP);

double
allocSplit(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
           const int* rM,             const List<int>* invrM,  const  int* mixtureNM,
           const double* propwM,      const double* propmuM,   const double* propinvsigma2M,
           const int jsplit,          const int* kP,           
           const double* regresResM,  const double* Eb0,       const int* randomIntP);

double
allocCombine(int* proprM,               List<int>* propinvrM,    int* propmixtureNM,
             const int* rM,             const List<int>* invrM,  const int* mixtureNM,
             const double* wM,          const double* muM,       const double* invsigma2M,
             const int jsplit,          const int* propkP,
             const double* regresResM,  const double* Eb0,       const int* randomIntP);

double 
logPostRatioSplitCombine(const int jsplit,              const int* shortkP,
                         const double* longwM,          const double* shortwM,     
                         const double* longmuM,         const double* shortmuM,
                         const double* longinvsigma2M,  const double* shortinvsigma2M,
                         const int* longmixtureNM,      const int* shortmixtureNM,
                         const double* deltaP,          const double* xiP,              const double* invkappaP,                   
                         const double* halfl2pikappaP,  const double* zetaP,            const double* etaP,                        
                         const double* lgammazetaP,     const double* llambdaP,         const int* priorForkP);

double
logJacobianSplitCombine(const double w,
                        const double mu0,           const double mu1,
                        const double invsigma20,    const double invsigma21,
                        const double invsigma2,     const double* vM);

#endif
