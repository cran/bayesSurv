#ifndef _UPDATE_ETA_H_
#define _UPDATE_ETA_H_

#include <R.h>
#include <Rmath.h>

#include "constants.h"

void
updateEta(double* etaP,
          const int* kP,        const double* invsigma2M,
          const double* zetaP,  const double* gP,           const double* hP);

#endif
