#ifndef _UPDATE_ALLOC_H_
#define _UPDATE_ALLOC_H_

#include <R.h>
#include <Rmath.h>

#include "List.h"
#include "random.h"

void
updateAlloc(int* rM,                List<int>* invrM,          int* mixtureNM, 
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP);

void
updateAlloc(int* rM,                
            const double* wM,       const double* muM,         const double* invsigma2M,
            const int* kP,          const double* regresResM,  const double* Eb0,
            const int* randomIntP,  const int* nP);

#endif
