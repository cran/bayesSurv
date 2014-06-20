#ifndef _UPDATE_DATA_GS_H_
#define _UPDATE_DATA_GS_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "constants.h"
#include "AK_Error.h"
#include "Gspline.h"

const double _ZERO_TIME_     = 1e-10;
const double _LOG_ZERO_TIME_ = log(_ZERO_TIME_);

// =======================================================================================================
// ***** Header for update_Data_GS.cpp: update possibly censored data                             ***** //
// =======================================================================================================
void
update_Data_GS(double* YsM,
               const double*  y_left,  
               const double*  y_right,   
               const int*     status,
               const int*     rM,         
               const Gspline* gg,
               const int*     nP,         
               const int*     n_censored);

void
update_Data_GS_regres(double* YsM,           
                      double* regresResM,
                      const double*  y_left,  
                      const double*  y_right,   
                      const int*     status,
                      const int*     rM,         
                      const Gspline* gg,       
                      const int* nP);

void
update_Data_GS_doubly(double* Yevent,        
                      double* regresResM,
                      const double*  Yonset, 
  	              const double*  t_left,  
                      const double*  t_right,    
                      const int*     status,
                      const int*     rM,         
                      const Gspline* gg,       
                      const int*     nP);

void
update_Data_GS_regres_misclass(double* YsM,           
                               double* regresResM,
                               int*    n00,
                               int*    n10,
                               int*    n01,
                               int*    n11,
                               double* iPML,
                               double* dwork,
                               const double*  sens,
                               const double*  spec,
                               const double*  logvtime,
                               const int*     status,
                               const int*     nExaminer,
                               const int*     nFactor,
                               const int*     nvisit,
                               const int*     maxnvisit,
                               const int*     Examiner,
                               const int*     Factor,
                               const int*     rM,         
                               const Gspline* gg,       
                               const int*     nP);

extern "C"{

void
iPML_misclass_GJK(double* iPML,                  
                  double* dwork,
                  const double* etaM,
                  const double* sens,
                  const double* spec,
                  const double* logvtime,
                  const int*    status,
                  const int*    nExaminer,
                  const int*    nFactor,
                  const int*    nvisit,
                  const int*    maxnvisit,
                  const int*    Examiner,
                  const int*    Factor,
                  const int*    gg_K,
                  const double* gg_gamma,
                  const double* gg_delta,         
                  const double* gg_sigma,        
                  const double* gg_intcpt,        
                  const double* gg_scale,
                  const double* gg_w,        
                  const int*    nP);

}         // end of extern "C"

#endif
