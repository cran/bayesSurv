#ifndef _UPDATE_SENS_SPEC_MISCLASSIFICATION_H_
#define _UPDATE_SENS_SPEC_MISCLASSIFICATION_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "AK_Error.h"

const double ZERO_SENS_SPEC = 1e-6;

void
update_sens_spec_misclassification(double* sens,
                                   double* spec,
                                   const double* priorPar,
                                   const int* n00,
                                   const int* n10,
                                   const int* n01,
                                   const int* n11,
                                   const int* nExaminer,
                                   const int* nFactor);

#endif
