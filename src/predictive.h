// Header for predictive.cpp: Sampling of predictive survivor, hazard, cummulative hazard curves and predictive survivor times
//
#ifndef _PREDICTIVE_H_
#define _PREDICTIVE_H_

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "List.h"
#include "in_output.h"
#include "arrays2Mat.h"
#include "covMatrix.h"
#include "bblocks.h"
#include "miscellaneous.h"
#include "mixMoments.h"
#include "random.h"
#include "predictMisc.h"
#include "quantile.h"

extern "C" {

void
predictive(const int* errorTypeP,  const char** dirP,        int* dimsP,          int* dims2P, 
           const double* XA,       const int* indbA,         double* quantileA,   double* gridA,
           const int* priorParI,   const double* priorParD,  const int* nsimulP,  const int* skipP,          
           const int* byP,         int* onlyAver,            int* predictP,       int* storeP,
           const double* tolersP,  int* errP);

}

#endif
