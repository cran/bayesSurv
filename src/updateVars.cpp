// Function to update the mixture inverse variances
//   Gibbs step is used here.

// 13/01/2004: start woking on it
// 14/03/2004: general mean of the random intercept allowed

#include "bayessurvreg.h"

using namespace std;

// =======================================
//
// INPUT PARAMETERS:
//
// Eb0 ............ mean of the random intercept                                (1) 
//                  (it is not changed by this function but it cannot be 'const'
//                   since often Eb0 = mixMomentM + 0 and part of mixMomentM is recalculated)
// regresResM ..... regression residuals (y - x'beta - z'b)                     (nP x 1)
// kP ............. current number of mixture components
// mixtureNM ...... numbers of observations belonging to each mixture component (at least kP x 1)
// muM ............ curent mixture means                                        (at least kP x 1)
// rM ............. current component pertinences                               (nP x 1)
// zetaP .......... prior hyperparameter  (shape of inverse-gamma)
// etaP ........... hyperparameter (scale of inverse-gamma, rate of gamma)
// randomIntP ..... 0/1 indicating whether a rando intercept is in the model
// nP ............. number of observations
//
void
updateVars(double* invsigma2M,        double* mixMomentM,      double* Eb0,
           const double* regresResM,  
           const int* kP,             const int* mixtureNM,
           const double* wM,          const double* muM,       const int* rM, 
           const double* zetaP,       const double* etaP,  
           const int* randomIntP,     const int* nP)
{
  int obs, j;
  double intcptadd;

//double* resid = new double[*nP];
//double* pars = new double[*kP *2];

  if (*randomIntP) intcptadd = *Eb0;
  else             intcptadd = 0.0;  

  // Compute the rate of the proposal distribution
  double* proposalShape = new double[*kP];
  double* proposalScale = new double[*kP];
  for (j = 0; j < *kP; j++){
    proposalShape[j] = *zetaP;
    proposalScale[j] = 0.0;
  }
  for (obs = 0; obs < *nP; obs++){
//resid[obs] = regresResM[obs] - muM[rM[obs]];
    proposalScale[rM[obs]] += ((regresResM[obs] - muM[rM[obs]] + intcptadd) * (regresResM[obs] - muM[rM[obs]] + intcptadd));
  }

//writeToFile(proposalScale, 1, *kP, "/home/arnost/temp", "/sumsqbijVar.sim", 'a');
//writeToFile(resid, 1, *nP, "/home/arnost/temp", "/residbijVar.sim", 'a');
//writeToFile(mixtureNM, 1, *kP, "/home/arnost/temp", "/NbijVar.sim", 'a');


  // Adjust further the proposal parameters and sample
  for (j = 0; j < *kP; j++){
    proposalScale[j] *= 0.5;
    proposalScale[j] += (*etaP);                   // this is now the rate of the proposal distrib.
    proposalScale[j] = 1/proposalScale[j];         // this is now its scale
    proposalShape[j] += (0.5 * mixtureNM[j]);
    if (proposalScale[j] <= SCALE_ZERO) proposalScale[j] = SCALE_ZERO;
    invsigma2M[j] = rgamma(proposalShape[j], proposalScale[j]);    
//pars[2*j] = proposalShape[j];
//pars[2*j + 1] = 1/proposalScale[j];
  }
//writeToFile(pars, 1, *kP*2, "/home/arnost/temp", "/parsbijVar.sim", 'a');

//delete [] pars;
//delete [] resid;

    // Recalculate SD of the error
  mixMoments(mixMomentM, kP, wM, muM, invsigma2M, true);

  delete [] proposalShape;
  delete [] proposalScale;
  return;
}   // end of function updateVars


