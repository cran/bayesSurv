// Functions to update the latent data (if they are censored)
//   and to sample predictive data
// * Gibbs step is used

// 24/11/2003: start woking on it
// 06/02/2004: robustness towards out of probability scale data
// 13/03/2004: allow a general mean of the random intercept
// 09/12/2004: check for generated NaN added

#include "updateData.h"

// ======================================================================
//
// PARAMETERS:
//
// YsM ........ on INPUT: current vector of (imputed) log(event times)                (nP x 1)
//              on OUTPUT: updated vector
// regresResM . on INPUT: current vector of regression residuals (y - x'beta - z'b))  (nP x 1)
//              on OUTPUT: updated vector
//
// INPUT PARAMETERS:
//
// Y1M ........ lower limit of interval censored obs.    (nP x 1)
// Y2M ........ upper limit of interval censored obs.    (nP x 1)  
// statusM .... censoring status                         (nP x 1)
// rM ......... vector of pertainence indicators         (nP x 1)
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
           const int* rM,           double* cumwM,         const double* muM,     const double* invsigma2M,
           const double* Eb0,       const int* kP,         const int* nP,        
           const int* errorTypeP,   const int* randomIntP)
{
  int obs, j;
  double PhiL = 0;
  double PhiU = 0;  
  double u = 0;
  double PhiInv = 0;
  double stres = 0;
  double intcptadd = 0;

  double* invsigmaM = new double[*kP];
  double* sigmaM = new double[*kP];
  for (j = 0; j < *kP; j++){
    invsigmaM[j] = sqrt(invsigma2M[j]);
    sigmaM[j] = 1/invsigmaM[j];
  }

  double* y_obs = YsM;
  double* regRes = regresResM;
  const double* y1 = Y1M;
  const double* y2 = Y2M;
  const int* status = statusM;
  const int* rp = rM;

  switch (*errorTypeP){
  case Mixture:
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      switch (*status){
      case 1:    // precisely observed observation
        *y_obs = *y1;
        break;

      case 0:    // right censored observation
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - muM[*rp] + intcptadd) * invsigmaM[*rp];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        if (PhiL >= 1 - NORM_ZERO){        // censored time irrealistic large (out of the prob. scale)
          *y_obs = *y1;
        }
        else{
          if (PhiL <= NORM_ZERO){         // censoring time equal to "zero", generate an exact time from N(mean, variance), 
                                            //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv; 
          }
          else{
            u = runif(0, 1) * (1 - PhiL) + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_PosInf){       // u was equal to 1, additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs =  -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv;
            }
          }  
        }            
        *regRes += (*y_obs);
        break;

      case 2:    // left censored observation
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - muM[*rp] + intcptadd) * invsigmaM[*rp];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        if (PhiU <= NORM_ZERO){           // left censoring time irrealistic low (equal to "zero")
          *y_obs = *y1;
        }
        else{
          if (PhiU >= 1 - NORM_ZERO){    // left censoring time equal to "infty", generate an exact time from N(mean, variance), 
                                           //   i.e. from the full  not-truncated distribution
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiU;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (PhiInv == R_NegInf){      // u was equal to 0,  additional check added 16/12/2004
              *y_obs = *y1;
            }
            else{
              *y_obs =  -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv;
            }
          }
        }
        *regRes += *y_obs;     
        break;

      case 3:    // interval censored observation
        *regRes -= *y_obs;
        stres = (*y1 + (*regRes) - muM[*rp] + intcptadd) * invsigmaM[*rp];
        PhiL = pnorm(stres, 0, 1, 1, 0);
        stres = (*y2 + (*regRes) - muM[*rp] + intcptadd) * invsigmaM[*rp];
        PhiU = pnorm(stres, 0, 1, 1, 0);
        PhiInv = PhiU - PhiL;
        if (PhiInv <= NORM_ZERO){       // too narrow interval, or the interval out of the probability scale
                                        //   (both limits in "zero" probability region)
                                        //   generate something inbetween
          u = runif(0, 1);
          *y_obs = *y1 + u*((*y2) - (*y1)); 
        }
        else{
          if (PhiInv >= 1 - NORM_ZERO){ // too large interval, practically (-infty, +infty), generate an exact time from N(mean, variance)
            u = runif(0, 1);
            PhiInv = qnorm(u, 0, 1, 1, 0);
            *y_obs = -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv; 
          }
          else{
            u = runif(0, 1) * PhiInv + PhiL;
            PhiInv = qnorm(u, 0, 1, 1, 0);
            if (!R_finite(PhiInv)){    // u was either zero or one,  additional check added 16/12/2004
              u = runif(0, 1);
              *y_obs = *y1 + u*((*y2) - (*y1)); 
            }
            else{
              *y_obs =  -(*regRes) + muM[*rp] - intcptadd + sigmaM[*rp]*PhiInv;
            }
          }
	}
        *regRes += *y_obs;     
        break;
      }   // end of switch 

      /*** This section performs just additional check to prevent simulation working with NaNs ***/
      if (!R_finite(*y_obs) || !R_finite(*regRes)){
        int condit;
        REprintf("\nY[%d]=%e,  regRes[%d]=%e, r[%d]=%d,  status[%d]=%d,  stres=%e", 
                  obs, *y_obs, obs, *regRes, obs, *rp, obs, *status, stres);
        REprintf("\nk=%d", *kP);
        REprintf("\nmean="); for (j=0; j < *kP; j++) REprintf("%e, ", muM[j]);
        REprintf("\ninvvar="); for (j=0; j < *kP; j++) REprintf("%e, ", invsigma2M[j]);
        REprintf("\nsigma="); for (j=0; j < *kP; j++) REprintf("%e, ", sigmaM[j]);
        REprintf("\ninvsigma="); for (j=0; j < *kP; j++) REprintf("%e, ", invsigmaM[j]);
        REprintf("\nu=%3.20e,  PhiL=%3.20e,  PhiU=%3.20e,  PhiInv=%3.20e", u, PhiL, PhiU, PhiInv);
        REprintf("\nNORM_ZERO=%3.20e,  1-NORM_ZERO=%3.20e", NORM_ZERO, 1-NORM_ZERO);
        switch (*status){
        case 0:
          condit = 1*(PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiL <= NORM_ZERO);
          REprintf("\nPhiL <= NORM_ZERO: %d", condit);
          break;
        case 2:
          condit = 1*(PhiU >= 1 - NORM_ZERO);
          REprintf("\nPhiU >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU <= NORM_ZERO);
          REprintf("\nPhiU <= NORM_ZERO: %d", condit);
          break;
        case 3:
          condit = 1*(PhiU-PhiL >= 1 - NORM_ZERO);
          REprintf("\nPhiU-PhiL >= 1 - NORM_ZERO: %d", condit);
          condit = 1*(PhiU-PhiL <= NORM_ZERO);
          REprintf("\nPhiU-PhiL <= NORM_ZERO: %d", condit);
          break;
        }        
        REprintf("\n");
        delete [] invsigmaM;
        delete [] sigmaM;
        throw returnR("Trap in updateData: NaN generated.", 1);
      }

      y_obs++;
      regRes++;
      y1++;
      y2++;
      status++;
      rp++;
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




