// =======================================================
// ***** moodyRing.h: Routines for a moody ring *****
// =======================================================
#ifndef _MOODY_RING_H_
#define _MOODY_RING_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"

extern "C"{

void
moodyRing(double* uA,              double* moodP,             
          const double* timeDepP,  const double* componentDepP,     
          const int* nuP,          const int* corrP,
          const int* callFromR);

void
corr_moodyRing(double* uA,              double* moodA,                double* initmoodP,
               const double* timeDepP,  const double* componentDepP,  const int* nuP,
               const int* nP,           const int* callFromR);

void
indep_moodyRing(double* uA,               double* inituA,
                const double* timeDepP,   const int* nuP,
                const int* nP,            const int* callFromR);


}

#endif

