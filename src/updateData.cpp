// Functions to update the latent data (if they are censored)
//   and to sample predictive data
// * Gibbs step is used

// 24/11/2003: start woking on it
// 06/02/2004: robustness towards out of probability scale data
// 13/03/2004: allow a general mean of the random intercept

#include "bayessurvreg.h"

using namespace std;


// ======================================================================
//
// PARAMETERS:
//
// YsM ........ on INPUT: current vector of (imputed) log(event times)                (nP x 1)
//              on OUTPUT: updated vector
// regresResM . on INPUT: current vectoe of regression residuals (y - x'beta - z'b))  (nP x 1)
//              on OUTPUT: updated vector
//
// INPUT PARAMETERS:
//
// Y1M ........ lower limit of interval censored obs.    (nP x 1)
// Y2M ........ upper limit of interval censored obs.    (nP x 1)  
// statusM .... censoring status                         (nP x 1)
// rM ......... vector of pertainence indicators         (nP x 1)
//              (when *errorTypeP = Spline, this is a working space)
// cumwM ...... cumulative mixture weights               (nP x 1)
//              (it is used if errorTypeP = Spline)
//              (it is not 'const' but it is not modified by this function)
// muM ........ vector of component means                (at least kP x 1)             
// invsigma2M . vector of component inverse-variances    (at least kP x 1)
// Eb0 ........ mean of the random intercept
// kP ......... current number of mixture components     (1)
// nP ......... number of observations
// errorTypeP . indicator whether the error density is assumed to be a spline 
//              (i.e. mixtureNM and rM are then ignored
//              or a classical mixture where it is assumed that each residual belongs to one mixture component
//              or something else
// randomIntP . 0/1 indicating whether a random intercept is in the model
//
void
updateData(double* YsM,             double* regresResM,
           const double* Y1M,       const double* Y2M,     const int* statusM,
           int* rM,                 double* cumwM,         const double* muM,     const double* invsigma2M,
           const double* Eb0,       const int* kP,         const int* nP,        
           const int* errorTypeP,   const int* randomIntP)
{
  int obs, j;
  double PhiL, PhiU;  
  double u, PhiInv, stres, intcptadd;

  double* invsigmaM = new double[*kP];
  double* sigmaM = new double[*kP];
  for (j = 0; j < *kP; j++){
    invsigmaM[j] = sqrt(invsigma2M[j]);
    sigmaM[j] = 1/invsigmaM[j];
  }

  switch (*errorTypeP){
  case Mixture:
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      switch (statusM[obs]){
      case 1:    // precisely observed observation
        YsM[obs] = Y1M[obs];
        break;

      case 0:    // right censored observation
        regresResM[obs] -= YsM[obs];
        stres = (Y1M[obs] + regresResM[obs] - muM[rM[obs]] + intcptadd) * invsigmaM[rM[obs]];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        if (PhiL >= 1 - NORM_ZERO)        // censored time irrealistic large (out of the prob. scale)
          YsM[obs] = Y1M[obs];
        else
          if (PhiL <= NORM_ZERO){         // censoring time equal to "zero", generate an exact time from N(mean, variance), 
                                            //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] = -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv; 
          }
          else{
            u = runif(0, 1) * (1 - PhiL) + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] =  -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv;
          }              
        regresResM[obs] += YsM[obs];
        break;

      case 2:    // left censored observation
        regresResM[obs] -= YsM[obs];
        stres = (Y1M[obs] + regresResM[obs] - muM[rM[obs]] + intcptadd) * invsigmaM[rM[obs]];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        if (PhiU <= NORM_ZERO)           // left censoring time irrealistic low (equal to "zero")
          YsM[obs] = Y1M[obs];
        else
          if (PhiU >= 1 - NORM_ZERO){    // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                           //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] = -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiU;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] =  -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv;
          }
        regresResM[obs] += YsM[obs];     
        break;

      case 3:    // interval censored observation
        regresResM[obs] -= YsM[obs];
        stres = (Y1M[obs] + regresResM[obs] - muM[rM[obs]] + intcptadd) * invsigmaM[rM[obs]];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        stres = (Y2M[obs] + regresResM[obs] - muM[rM[obs]] + intcptadd) * invsigmaM[rM[obs]];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        PhiInv = PhiU - PhiL;
        if (PhiInv <= NORM_ZERO){       // too narrow interval, or the interval out of the probability scale
                                        //   (both limits in "zero" probability region)
                                        //   generate something inbetween
          u = runif(0, 1);
          YsM[obs] = Y1M[obs] + u*(Y2M[obs] - Y1M[obs]); 
        }
        else
          if (PhiInv >= 1 - NORM_ZERO){ // too large interval, practically (-infty, +infty), generate an exact time from N(mean, variance)
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] = -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiInv + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            YsM[obs] =  -regresResM[obs] + muM[rM[obs]] - intcptadd + sigmaM[rM[obs]]*PhiInv;
          }
        regresResM[obs] += YsM[obs];     
        break;
      }   // end of switch 
    }   // end of the loop over observations
    break;

  case Spline:
    returnR("C++ Error: Not yet implemented part (Spline) of function updateData called.", 1);
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function updateData called.", 1);
    break;
      
  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function updateData.", 1);
  }
    
  delete [] invsigmaM;
  delete [] sigmaM;

  return;
}   // end of function updateData 




