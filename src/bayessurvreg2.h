#ifndef _BAYES_SURV_REG_TWO_H_
#define _BAYES_SURV_REG_TWO_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "openFile.h"

#include "classCovMatrixExtend.h"
#include "classBetaGammaExtend.h"
#include "classRandomEff.h"

#include "regresResidual_GS.h"
#include "in_output_GS.h"
#include "update_Data_GS.h"
#include "update_Alloc_GS.h"
#include "update_sens_spec_misclassification.h"

//#include "AK_BLAS_LAPACK.h"
#include "rhoNorm.h"
#include "structRandomEff32.h"

const double _toler_chol_bayessurvreg2 = 1e-10;

extern "C"{

void
bayessurvreg2(char** dirP,             
              const int*    dimsP,         
              const double* X1,     
              const double* X2,
              const double* y1_left,   
              const double* y1_right,   
              const int*    status1,
              const double* t2_left,   
              const double* t2_right,   
              const int*    status2,
              double* iPML,               
              double* Y1,              
              double* Y2,
              int* r1,                 
              int* r2,                  
              const int* specif,    
              int* r_b1,               
              int* r_b2,                
              const int* specif_b,
              int*    GsplineI1,          
              double* GsplineD1,
              int*    GsplineI2,          
              double* GsplineD2,
              int*    priorBeta1I,        
              double* priorBeta1D,
              int*    priorBeta2I,        
              double* priorBeta2D,
              int*    priorb1I,           
              double* priorb1D,
              int*    priorb2I,           
              double* priorb2D,
              int*    priorCovMat1I,      
              double* priorCovMat1D,
              int*    priorCovMat2I,      
              double* priorCovMat2D,              
              double* rhob,            
              int*    rho_accept,          
              const int*    rhobI,     
              const double* rhobD,
              double*       mcsensspec,
              const double* mclogvtime,
              const int*    mcstatus,
              const int*    mcparI,
              const double* mcparD,
              int* iterM,              
              int* nsimulP,             
              int* storeP,
              const int* version,      
              const int* mainSimul,     
              int* errP);
}

void
adjust_intcpt(Gspline*   g_eps,  
              Gspline*   g_b,  
              RandomEff* bb);

void
print_iter_info(int& writeAll,  
                int& backs,  
                const int& iter,  
                const int& nwrite,  
                const int& lastIter);

#endif
