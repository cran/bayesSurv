// Header for bayesBisurvreg.cpp
#ifndef _BAYES_BI_SURV_REG_H_
#define _BAYES_BI_SURV_REG_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "openFile.h"

#include "constants.h"
#include "classBetaGamma.h"
#include "regresResidual_GS.h"
#include "in_output_GS.h"
#include "update_Data_GS.h"
#include "update_Alloc_GS.h"


#ifdef __cplusplus  
extern "C"{
#endif
  
void
bayesBisurvreg(char** dirP,             const int* dimsP,         const double* X1,     const double* X2,
               const double* y1_left,   const double* y1_right,   const int* status1,
               const double* t2_left,   const double* t2_right,   const int* status2,               
               double* Y1,              double* Y2,
               int* r1,                 int* r2,                  const int* specif,
               int* GsplineI1,          double* GsplineD1,
               int* GsplineI2,          double* GsplineD2,
               int* priorBeta1I,        double* priorBeta1D,
               int* priorBeta2I,        double* priorBeta2D,
               int* iterM,              int* nsimulP,             int* storeP,
               const int* mainSimul,    int* errP);

#ifdef __cplusplus    
}    /*** end of extern "C" ***/
#endif

#endif
