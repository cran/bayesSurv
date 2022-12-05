/*
** This file causes the entry points of my .C routines to be preloaded.
**
** Added on 29/03/2017 at the request of R CMD check.
**
** It adds one more layer of protection by declaring the number of arguments,
** and perhaps a tiny bit of speed.
*/
#include <R.h>
// #include <Rinternals.h>   // Causes some conflict of its macro LENGTH(x) with
                             // another macro of the same name which is present
                             // in one of the C++ headers
                             // included by openFile.h.
                             // But it seems that Rinternals.h are not needed
                             // for bayesSurv compilation.
#include <R_ext/Rdynload.h>

#include "bayesBisurvreg.h"
#include "bayesDensity.h"
#include "bayesGspline.h"
#include "bayesHistogram.h"
#include "bayessurvreg1.h"
#include "bayessurvreg2.h"
#include "cholesky.h"
#include "update_Data_GS.h"
#include "miscellaneous_GS.h"
#include "checkImputeData.h"
#include "marginal_bayesGspline.h"
#include "predictive_GS.h"
#include "predictive.h"
#include "Mvtdist3.h"
#include "sampledKendallTau.h"

static const R_CMethodDef Centries[] = {
    {"C_bayesBisurvreg",        (DL_FUNC) &bayesBisurvreg,         28},    // bayesBisurvreg
    {"C_bayesDensity",          (DL_FUNC) &bayesDensity,           18},    // bayesDensity
    {"C_bayesGspline",          (DL_FUNC) &bayesGspline,           19},    // bayesGspline
    {"C_bayesHistogram",        (DL_FUNC) &bayesHistogram,         15},    // bayesHistogram
    {"C_bayessurvreg1",         (DL_FUNC) &bayessurvreg1,          28},    // bayessurvreg1
    {"C_bayessurvreg2",         (DL_FUNC) &bayessurvreg2,          50},    // bayessurvreg2
    {"C_cholesky",              (DL_FUNC) &cholesky,                5},    // cholesky  
    {"C_iPML_misclass_GJK",     (DL_FUNC) &iPML_misclass_GJK,      21},    // update_Data_GS
    {"C_findClosestKnot",       (DL_FUNC) &findClosestKnot,         5},    // miscellaneous_GS
    {"C_midimputeData",         (DL_FUNC) &midimputeData,           6},    // checkImputeData    
    {"C_midimputeDataDoubly",   (DL_FUNC) &midimputeDataDoubly,    10},    // checkImputeData
    {"C_marginal_bayesGspline", (DL_FUNC) &marginal_bayesGspline,  18},    // marginal_bayesGspline
    {"C_predictive_GS",         (DL_FUNC) &predictive_GS,          39},    // predictive_GS
    {"C_predictive",            (DL_FUNC) &predictive,             18},    // predictive
    {"C_rmvnormR2006",          (DL_FUNC) &Mvtdist3::rmvnormR2006,  7},    // Mvtdist3
    {"C_rwishartR3",            (DL_FUNC) &Mvtdist3::rwishartR3,    6},    // Mvtdist3
    {"C_sampledKendallTau",     (DL_FUNC) &sampledKendallTau,      12},    // sampledKendallTau
    {NULL, NULL, 0}    
};

extern "C" void R_init_bayesSurv(DllInfo *dll){
    R_registerRoutines(dll, Centries, NULL,  NULL,     NULL);
    /*                      .C        .Call  .Fortran  .External         */
    
    /* The following line makes only those routines defined above
       available to outside packages, i.e., internal C++ things
       are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE);
 
    /*
    ** This line makes them only be available via the symbols above
    **  i.e., .C("bayesBisurvreg", ) won't work but .C(C_bayesBisurvreg, )  will
    */
    R_forceSymbols(dll, TRUE);
}
    
