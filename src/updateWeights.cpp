// Function to update the mixture weights
//   * Eb0 does not depend on mixture: Gibbs step is used to update weights
//   * Eb0 depends on mixture: Metropolis step is used where proposal is sampled as in the case
//     when Eb0 does not depend on mixture => proposal ratio is then ratio of log-likelihood and random log-likelihood

// 26/11/2003: start woking on it
// 15/03/2004: dependence of Eb0 on mixture allowed

#include "bayessurvreg.h"

using namespace std;

// ======================================================================
//
// INPUT PARAMETERS:
// 
// kP ......... current number of mixture components
// mixtureNM .. number of observations belonging to each mixture component
// deltaP ..... prior hyperparameter
//
void
updateWeights(double** wM,               double** propwM,            double* mixMomentM,
              double* loglikelhood,      double** loglikobs,         double** proploglikobs,
              double* randomloglik,      double** randomllcl,        double** proprandomllcl,
              const double* regresResM,  const double* YsM,          const double* betaM,        
              const double* bM,          const covMatrix* Dcm,
              const int* kP,             const int* mixtureNM,  
              const double* muM,         const double* invsigma2M,   const int* rM,
              const double* deltaP,
              const int* Eb0dependMix,   const int* randomIntP,      
              const int* nP,             const int* nclusterP,       const int* nrandomP,
              const int* indbinXA,       double (*logdtrans)(double))
{
  int errorType = Mixture;
  int j;
  double wsum = 0.0;

  if (*kP == 1) return;             // update is not necessary

  if ((*Eb0dependMix) && (*randomIntP)){
    double propEb0 = 0.0;
    double proploglik, proprandomloglik;
    double unif, Paccept;
    for (j = 0; j < *kP; j++){
      (*propwM)[j] = rgamma(*deltaP + double(mixtureNM[j]), 1.0);
      wsum += (*propwM)[j];
    }
    for (j = 0; j < *kP; j++){
      (*propwM)[j] /= wsum;
      propEb0 += (*propwM)[j]*muM[j];
    }
    logLikelihood(&proploglik, *proploglikobs, nP, regresResM, YsM, kP, rM, *propwM, muM, invsigma2M, &propEb0, 
                  &errorType, randomIntP, logdtrans);
    randomLogLikelihood(&proprandomloglik, *proprandomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaM, Dcm, &propEb0, indbinXA);
    Paccept = exp(proploglik + proprandomloglik - (*loglikelhood) - (*randomloglik));
    if (Paccept < 1.0){
      unif = runif(0, 1);
      if (unif > Paccept){
        return;
      }  
    }
    changePointers(wM, propwM);
    *loglikelhood = proploglik;
    *randomloglik = proprandomloglik;
    changePointers(loglikobs, proploglikobs);
    changePointers(randomllcl, proprandomllcl);
    mixMomentM[0] = propEb0;
    mixMoments(mixMomentM, kP, *wM, muM, invsigma2M, true);
  }

    // Eb0 does not depend on mixture or no random intercept in the model
  else{
    for (j = 0; j < *kP; j++){
      (*wM)[j] = rgamma(*deltaP + double(mixtureNM[j]), 1.0);
      wsum += (*wM)[j];
    }
    for (j = 0; j < *kP; j++) (*wM)[j] /= wsum;

    mixMoments(mixMomentM, kP, *wM, muM, invsigma2M, false);
  }

  return;
}


