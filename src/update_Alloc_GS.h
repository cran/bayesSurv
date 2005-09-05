// =======================================================================================================
// ***** Header for update_Alloc_GS.cpp: update labels of component pertinences                   ***** //
// =======================================================================================================
#ifndef _UPDATE_ALLOC_GS_H_
#define _UPDATE_ALLOC_GS_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "random.h"
#include "Gspline.h"

void
update_Alloc_GS(int* rM,            int* mixtureNM,            double** mu,    double* loglik,  double* logpr,
                const Gspline* gg,  const double* regresResM,  const int* nP,   
                int* iwork,         double* dwork);

#endif
