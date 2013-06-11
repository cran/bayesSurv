//
// Function to update sensitivities and specificities in a misclassification model
//
// 31/05/2013:  implemented (in Santiago de Chile)
//
// =================================================================================
#include "update_sens_spec_misclassification.h"


// ***** update_sens_spec_misclassification *******************************
//
//  sens[nFactor * nExaminer]  INPUT:  whatsever
//                            OUTPUT:  updated sensitivities
//
//  spec[nFactor * nExaminer]  INPUT:  whatsever
//                            OUTPUT:  updated specificities
//
//  priorPar[4]                a.sens, b.sens, a.spec, b.spec: parameters of beta priors for sensitivities and specificities
//  n00[nExaminer * nFactor]   numbers of (0-0) correctly classified events for each examiner:factor
//  n10[nExaminer * nFactor]   numbers of (Classification = 1 | True = 0) incorrectly classified events for each examiner:factor
//  n01[nExaminer * nFactor]   numbers of (Classification = 0 | True = 1) incorrectly classified events for each examiner:factor
//  n11[nExaminer * nFactor]   numbers of (1-1) correctly classified events for each examiner:factor
//
void
update_sens_spec_misclassification(double* sens,
                                   double* spec,
                                   const double* priorPar,
                                   const int* n00,
                                   const int* n10,
                                   const int* n01,
                                   const int* n11,
                                   const int* nExaminer,
                                   const int* nFactor)
{
  double a, b, F_lower, u;

  const double *a_sens = priorPar;
  const double *b_sens = a_sens + 1;
  const double *a_spec = b_sens + 1;
  const double *b_spec = a_spec + 1;

  double *sensP = sens;
  double *specP = spec;
  const int *n00P = n00;
  const int *n10P = n10;
  const int *n01P = n01;
  const int *n11P = n11;

  for (int i = 0; i < *nFactor * *nExaminer; i++){
    /*** Update sensitivity ***/
    a = *a_sens + *n11P;
    b = *b_sens + *n01P;    
    if (*specP < ZERO_SENS_SPEC) *sensP = 1.0;
    else{
      F_lower = pbeta(1 - *specP, a, b, 1, 0);    // F(1 - spec), F = CDF of Beta(a, b)
      u = runif(0, 1);
      *sensP = qbeta(u + (1 - u) * F_lower, a, b, 1, 0);
    }    

    /*** Update specificity ***/
    a = *a_spec + *n00P;
    b = *b_spec + *n10P;         
    if (*sensP < ZERO_SENS_SPEC) *specP = 1.0;
    else{
      F_lower = pbeta(1 - *sensP, a, b, 1, 0);    // F(1 - sens), F = CDF of Beta(a, b)
      u = runif(0, 1);
      *specP = qbeta(u + (1 - u) * F_lower, a, b, 1, 0);
    }    

    /*** Shift the pointers ***/
    sensP++;
    specP++;
    n00P++;
    n10P++;
    n01P++;
    n11P++;
  }

  return;
}


                                   
                                   
