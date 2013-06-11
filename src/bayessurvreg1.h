// Header for bayessurvreg1.cpp
// ============================
#ifndef _BAYESSURVREG_ONE_H_
#define _BAYESSURVREG_ONE_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"

#include "MHblocks.h"
#include "bblocks.h"
#include "covMatrix.h"
#include "List.h"

#include "arrays2Mat.h"
#include "miscellaneous.h"
#include "in_output.h"
#include "mixMoments.h"
#include "logLikelihood.h"
#include "randomLogLikelihood.h"

#include "splitCombine.h"
#include "birthDeath.h"
#include "updateWeights.h"
#include "updateMeans.h"
#include "updateVars.h"
#include "updateEta.h"
#include "updateAlloc.h"
#include "updateData.h"
#include "updateRegres.h"
#include "updateRandom.h"
#include "updateCovMatRandom.h"

#include "moodyRing.h"
#include "propVector.h"


// =========================================================
// Log of the Jacobian of the transformation
//   * for different transformations
// =========================================================
inline double
logdlogtrans(const double y){ return -y;}

inline double
logdidenttrans(const double y){ return 0;}

extern "C"{

void 
bayessurvreg1(char** dirP,
              int* dimsP,
              double* YA,           double* XA,            int* indbA,
              int* iterM,           double* loglikM,       double* mixtureM,   double* mixMomentM,
              double* betaM,        double* bM,            double* DM,
              int* rM,              double* YsM,           double* otherpM,
              double* uM,
              int* priorParI,       double* priorParD,
              int* revJumpParI,     double* revJumpParD,
              int* priorBetaI,      double* priorBetaD,
              int* priorbI,         double* priorbD,
              int* nsimulP,         int* storeP,           double* tolersP,
              int* errP);

}

#endif
