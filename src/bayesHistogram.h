// ***** Header for bayesHistrogram.cpp: Main routine to be called from R to smooth a (bi)-variate density        ***** //
//                                       using G-splines                                                          ***** //
#ifndef _BAYES_HISTOGRAM_H
#define _BAYES_HISTOGRAM_H

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "openFile.h"
#include "Gspline.h"
#include "in_output_GS.h"
#include "templatefun_GS.h"
#include "update_Alloc_GS.h"
#include "update_Data_GS.h"

extern "C" {

void
bayesHistogram(char** dirP,           int* dimsP,
               double* y_left,        double* y_right,      int* status,
               int* rM,               double* YsM,
               int* iterM,
               const int* specif,     int* GsplineI,        double* GsplineD,
               int* nsimulP,          int* storeP,
               const int* mainSimul,  int* errP);

}

#endif
