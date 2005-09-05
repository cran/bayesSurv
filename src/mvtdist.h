// ======================================================================================
// ***** Header for mvtdist.cpp: Routines for dealing with the multivariate distributions
// ======================================================================================
#ifndef _MVT_DIST_H_
#define _MVT_DIST_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "cholesky.h"
//#include "classes.h"
//#include "random.h"

extern "C"{

void
LxMxtL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI);

void
tLxMxL(double* LAL,     const double* L,   const double* A,
       const int* nA,   const int* diagI);

void
Mxa(double* Ma,     const double* a,  const double* A,  const int* inda,  
    const int* na,  const int* nA,    const int* diagI);

void
Mxa2(double* Ma,     const double* a,  const double* A,  const int* indA,  
     const int* na,  const int* nA,    const int* diagI);

void
Wxa(double* Wa,   const double* a,   const double* A,  const int* indr,  const int* indc,
    const int* na,  const int* nA,  const int* nrow,   const int* diagI);

void
axMxa(double* aMa,    const double* a,  const double* A,  const int* inda,  
      const int* na,  const int* nA,    const int* diagI);

void
dmvtnorm(double* dens,      const double* x,      const double* mean,  const double* vari,  
         const int* indx,   const int* indxrepl, 
         const int* nx,     const int* nmean,     const int* nxrepl,    
         const int* nP,     const int* diagI,     const int* logP);

void
rmvtnorm(double* x,        const double* mean,   const double* L,  
         const int* indx,  const int* indxrepl,
         const int* nx,    const int* nmean,     const int* nxrepl,      
         const int* nP,    const int* diagI,     const int* callFromR);

void
rmvtnorm2(double* x,        const double* mean,   const double* iLi,  
          const int* indx,  const int* indxrepl,
          const int* nx,    const int* nmean,     const int* nxrepl,      
          const int* nP,    const int* diagI,     const int* callFromR);

void
rmvtiunif(double* x,        const double* mean,   const double* halfRange,  
          const int* indx,  const int* indxrepl,  
          const int* nx,    const int* nmean,     const int* nxrepl,
          const int* nP,    const int* callFromR);

void
rwishart(double* w,
         const int* p,      const double* nu,   const double* L,
         const int* diagI,  const int* nP,      const int* callFromR);

void
rinvwishart(double* w,
            const int* p,      const double* nu,   const double* L,
            const int* diagI,  const int* nP,      const int* callFromR);

void
rwishart2(double* w,
          const int* p,      const double* nu,   const double* iLi,
          const int* diagI,  const int* nP,      const int* callFromR);


}

#endif
