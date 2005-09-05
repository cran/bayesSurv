// ***** Header for arrays2Mat.cpp: Functions to create various matrices from input arrays *****
//
#ifndef _ARRAYS_2_MAT_H_
#define _ARRAYS_2_MAT_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "List.h"
#include "miscellaneous.h"

void
create2pointer(double**& point,    const int lpoint,  const int llayer2,  
               const double init,  const bool doinit);

void
create3pointer(double***& point,   const int lpoint,  const int* llayer2,  const int llayer3,  
               const double init,  const bool doinit);

void
create3pointer(double***& point,   const int lpoint,  const int llayer2,  const int* llayer3,  
               const double init,  const bool doinit);

void
createDataShort(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA,
                const double* XA,
                double** ZZt,           int* diagIZZt,             int* indbinXA,     
                const int* nP,          const int* nclusterP,
                const int* nXP,         const int* nfixedP,        const int* nrandomP,
                const int* randomIntP,  const int* indbA);

void 
createData(int* nwithinA,          int* clusteriA,            List<int>* invclusteriA, 
           int* statusA,           double*& Y1A,              double*& Y2A,           
           double** ZZt,           int* diagIZZt,             int* indbinXA,     
           double*** XXt,          int** diagIXXt,      
           const double* XA,       double* YA,
           const int* nP,          const int* nclusterP,      const int* nYP,
           const int* nXP,         const int* nfixedP,        const int* nrandomP,
           const int* randomIntP,  const int* indbA,      
           const int& nBlockbeta,  const int* nInBlockbeta,   int** indBlockbeta,        const int* typeUpdbeta);

void 
createParam(const int* nP,                 const int* kmaxP,             const double* mixtureA,        
            double* wM,                    double* muM,                  double* invsigma2M,
            int* rM,                       List<int>* invrM,             int* mixtureNM,  
            double* propwM,                double* propmuM,              double* propinvsigma2M,
            int* proprM,                   List<int>* propinvrM,         int* propmixtureNM);

void
createPriors(const int *kmaxP,     const double* priorParD,
             double* piSplitM,       
             double* logpiSplitM,  double* logpiCombineM,
             double* piBirthM,       
             double* logpiBirthM,  double* logpiDeathM);

#endif
