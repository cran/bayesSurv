// Function to compute a log-likelihood in the AFT model
//  with censored data replaced by imputed exact data,
//  i.e. all observations in this function are exactly observed.

// Various structures for the error term are allowed.

// 06/11/2003: start working on it
// 13/03/2004: changed handling of a random intercept

#include "bayessurvreg.h"

using namespace std;

// ======================================================================
//
// ***** logLikelihood *****
//
// * Function with two prototypes
//
// ======================================================================
//
// VERSION 1: Compute the value of the log-likelihood from scratch
//
// ======================================================================
//
// OUTPUT PARAMETERS:
//
// loglikelhood ...... value of the log-likelihood
// loglikobs ......... an array of personal contributions to the log-likelihood 
//
// INPUT PARAMETERS:
//
// nP ................ number of observations
// regresResM ........ regression residuals (y - x'beta - z'b)            (*nP)
// YsM ............... (augmented) data                                   (*nP)
// kP ................ number of mixture components             
// rM ................ vector of pertainence indicators                   (*nP)
//                     (ignored if *errorTypeP = 1)
// wM ................ vector of mixture weights                          (at least *kP)
//                     (ignored if *errorTypeP = 0)
// muM ............... vector of component means                          (at least *kP)
// invsigma2M ........ vector of component inverse-variances              (at least *kP)
// Eb0 ............... mean of the random intercept                       (1)
//                     (not always needed)
// errorTypeP ........ indicator of assumed error density (i.e. rM is then ignored) 
//                     or a classical mixture where it is assumed that each residual belongs to one mixture component
// randomIntP ........ 0/1 indicating whether random intercept is in the model
//                     if 1, it causes with some error structures (which are not of zero mean)
//                     that the mean of the random intercept must be explicitely added 
//                     to regression residuals before using them to compute a log-likelihood
// (*logdtrans) ...... -log of the Jacobian of the inverse transformation
//
void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double))
{
  int obs, j;
  bool minInfty = false;
  double intcptadd;
  double* halfLogInvSigmaMinLogS2PI;

  if (!(*nP)){
    *loglikelhood = 0.0;    
    return;
  }

  *loglikelhood = 0.0;              // reset it completely

  switch (*errorTypeP){
  case Mixture:
    halfLogInvSigmaMinLogS2PI = new double[*kP];
    for (j = 0; j < *kP; j++) halfLogInvSigmaMinLogS2PI[j] = 0.5 * log(invsigma2M[j]) - LOG_SQRT_2PI;

      // Personal contributions to the log-likelihood
    if (*randomIntP) intcptadd = *Eb0;
    else             intcptadd = 0.0;
    for (obs = 0; obs < *nP; obs++) {
      loglikobs[obs] = logdtrans(YsM[obs]);
      loglikobs[obs] = halfLogInvSigmaMinLogS2PI[rM[obs]];
      loglikobs[obs] -= 0.5 * ((regresResM[obs] - muM[rM[obs]] + intcptadd) * (regresResM[obs] - muM[rM[obs]] + intcptadd)) * invsigma2M[rM[obs]];
      if (loglikobs[obs] <= -FLT_MAX){
        loglikobs[obs] = -FLT_MAX;     // out of the probability scale
        minInfty = true;
      }
      *loglikelhood += loglikobs[obs];
    }

    if (minInfty == true) *loglikelhood = -FLT_MAX;
    delete [] halfLogInvSigmaMinLogS2PI;
    break;

  case Spline:
    returnR("C++ Error: Not yet implemented part (Spline) of function logLikelihood called.", 1);
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function logLikelihood called.", 1);
    break;
      
  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function logLikelihood.", 1);
  }

  return;
}


// ===============================================================================
//
// VERSION 2: Update only log-likelihood contributions of desired observations
//            and appropriately also the total log-likelihood
//
// ===============================================================================
//
// INPUT PARAMETERS:
//
// obsUpd ............ indeces of observations for which the log-likelihood contributions are to be updated
// 
void
logLikelihood(double* loglikelhood,       double* loglikobs,
              const List<int>* obsUpd,    const int* nP,
              const double* regresResM,   const double* YsM,
              const int* kP,              const int* rM,                
              const double* wM,           const double* muM,      const double* invsigma2M,
              const double* Eb0,          const int* errorTypeP,  const int* randomIntP,
              double (*logdtrans)(const double))
{
  int obs, obb, j;
  bool minInfty = false;
  double intcptadd;
  double* halfLogInvSigmaMinLogS2PI;

  int nobsUpd = obsUpd->length();

  if (!(*nP)){
    *loglikelhood = 0.0;    
    return;
  }

  switch (*errorTypeP){
  case Mixture:
    halfLogInvSigmaMinLogS2PI = new double[*kP];
    for (j = 0; j < *kP; j++) halfLogInvSigmaMinLogS2PI[j] = 0.5 * log(invsigma2M[j]) - LOG_SQRT_2PI;

    if (*randomIntP) intcptadd = *Eb0;
    else             intcptadd = 0.0;
    for (obb = 0; obb < nobsUpd; obb++){
      obs = (*obsUpd)[obb];
      *loglikelhood -= loglikobs[obs];              // subtract the old value
      loglikobs[obs] = logdtrans(YsM[obs]);
      loglikobs[obs] = halfLogInvSigmaMinLogS2PI[rM[obs]];
      loglikobs[obs] -= 0.5 * ((regresResM[obs] - muM[rM[obs]] + intcptadd) * (regresResM[obs] - muM[rM[obs]] + intcptadd)) * invsigma2M[rM[obs]];
      if (loglikobs[obs] <= -FLT_MAX){
        loglikobs[obs] = -FLT_MAX;     // out of the probability scale
        minInfty = true;
      } 
      *loglikelhood += loglikobs[obs];             // add the new value
    }

    if (minInfty == true) *loglikelhood = -FLT_MAX;
    delete [] halfLogInvSigmaMinLogS2PI;
    break;

  case Spline:
    returnR("C++ Error: Not yet implemented part (Spline) of function logLikelihood called.", 1);
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function logLikelihood called.", 1);
    break;
     
  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function logLikelihood.", 1);
  }

  return;

}   // end of function logLikelihood


// ======================================================================
//
// ***** clusterlogLikelihood *****
//
// Compute a value of the log-likelihood in a specific cluster only
//  (this function only sums personal log-likelihood contributions)
//
// =======================================================================
//
// loglikobs ......... an array of personal contributions to the log-likelihood 
// cl ................ index of a cluster
// obsInCluster ...... indeces of observations within that cluster
// 
void
clusterlogLikelihood(double* clusterloglik,  const double* loglikobs,  const int* cl,  const List<int>* obsInCluster)
{
  int i, obs;

  int nobs = obsInCluster->length();

  *clusterloglik = 0.0;
  for (i = 0; i < nobs; i++){
    obs = (*obsInCluster)[i];
    if (loglikobs[obs] <= -FLT_MAX){
      *clusterloglik = -FLT_MAX;
      return;
    }
    *clusterloglik += loglikobs[obs]; 
  }

  return;
}

