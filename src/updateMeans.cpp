// Function to update the mixture means
//   The components are updated using a Metropolis-Hastings step
//   which is almost equal to the Gibbs move.
//   The mixture means are sampled from their full conditionals (apart the factorial)
//   and accepted provided the ordering is unchanged.

// 26/11/2003: start woking on it
// 14/03/2004: general mean of the random intercept allowed

// Assumption: input means (muM) are sorted in increasing order

// Remarks: * with mean of random intercept being a mean of the mixture, this mean must
//            be recalculated after each mixture mean is updated
//          * moreover, the full conditional distribution is somewhat more complex
//            since also distribution of random effects must be taken into account

#include "updateMeans.h"

using namespace std;

// =======================================
//
// mixMomentM ..... mean and standard deviation of the error distribution
//
// INPUT PARAMETERS:
//
// regresResM ..... regression residuals (y - x'beta - z'b)                     (nP x 1)
// kP ............. current number of mixture components
// mixtureNM ...... numbers of observations belonging to each mixture component (at least kP x 1)
// invsigma2M ..... current mixture inverse variances                           (at least kP x 1)
// rM ............. current component pertinences                               (nP x 1)
// xiInvkappaP .... prior hyperparameter (mean * inverse variance)
// invkappaP ...... prior hyperparameter (inverse variance)               
// nP ............. number of observations
//
void
updateMeans(double* muM,                double* mixMomentM,
            const double* regresResM,   const double* betaM,        const double* bM,
            const covMatrix* Dcm,
            const int* kP,              const int* mixtureNM,
            const double* wM,           const double* invsigma2M,   const List<int>* invrM,  
            const double* xiInvkappaP,  const double* invkappaP,  
            const int* Eb0dependMix,    const int* randomIntP,      
            const int* nP,              const int* nclusterP,       const int* nrandomP,
            const int* indbinXA)
{
  int i, j, obs;
  double proposalMean, proposalSD, proposalmu;
  double intcptadd = 0.0;                              // alpha[j] in my notes
  double sumb0 = 0.0;                                  // sum(b[i,0])
  double sumW = 0.0;                                   // sum(V[0,1]*(b[i] - gamma))
  bool accept;

  if ((*Eb0dependMix) && (*randomIntP)){
    double* sumbmingamma = new double[*nrandomP - 1];          // sum(b[i] - gamma)
    for (j = 1; j < *nrandomP; j++){
      sumbmingamma[j] = -(*nclusterP)*betaM[indbinXA[j]];
    }
    for (i = 0; i < *nclusterP; i++){
      sumb0 += bM[(*nrandomP)*i];
      for (j = 1; j < *nrandomP; j++){
        sumbmingamma[j] += bM[(*nrandomP)*i + j];
      }
    }
    for (j = 1; j < *nrandomP; j++){                          // sum(V[0,1] * (b[i] - gamma))
      sumW += Dcm->icovm[j] * sumbmingamma[j - 1];
    }
    delete [] sumbmingamma;
  }

  for (j = 0; j < *kP; j++){
    accept = false;

      // Compute the mean and variance of the proposal distributions
    proposalMean = 0.0;
    if ((*Eb0dependMix) && (*randomIntP)) intcptadd = mixMomentM[0] - wM[j]*muM[j];
    for (i = 0; i < invrM[j].length(); i++){
      obs = invrM[j][i];
      proposalMean += regresResM[obs] + intcptadd;
    }
    if ((*Eb0dependMix) && (*randomIntP)){
      proposalSD = 1 / ((*invkappaP) + mixtureNM[j]*invsigma2M[j]*(1 - wM[j])*(1 - wM[j]) +
                   (*nclusterP)*wM[j]*wM[j]*Dcm->icovm[0]);
      proposalMean *= (1 - wM[j])*invsigma2M[j];
      proposalMean += (*xiInvkappaP) + wM[j]*Dcm->icovm[0]*(sumb0 - (*nclusterP)*intcptadd) +
                      wM[j]*sumW;
    }
    else{
      proposalSD = 1 / ((*invkappaP) + mixtureNM[j] * invsigma2M[j]);
      proposalMean *= invsigma2M[j];
      proposalMean += (*xiInvkappaP); 
    }
    proposalMean *= proposalSD;
    proposalSD = sqrt(proposalSD);

      // Sample and check the adjacency condition
      // *kP == 1
      // =========
    if (*kP == 1){
      proposalmu = rnorm(proposalMean, proposalSD);
      muM[0] = proposalmu;
      mixMoments(mixMomentM, kP, wM, muM, invsigma2M, false);
      return;
    }

      // *kP >= 2
      // =========
    // j = 0 (only one neighbour on the right)
    if (j == 0){
      proposalmu = rnorm(proposalMean, proposalSD);    
      if (proposalmu < muM[1]) accept = true;
    }
    else{
      // j = 1, ..., *kP - 2 (neighbours from both sides)
      if (j < *kP - 1){
        proposalmu = rnorm(proposalMean, proposalSD);
        if (proposalmu > muM[j - 1] && proposalmu < muM[j + 1]) accept = true;
      }
      // j = *kP - 1 (only one neighbour on the left)      
      else{
        proposalmu = rnorm(proposalMean, proposalSD);    
        if (proposalmu > muM[*kP - 2]) accept = true;
      }
    }

    // Recalculate mean of the error if necessary
    if (accept){
      mixMomentM[0] += wM[j]*(proposalmu - muM[j]);
      muM[j] = proposalmu;
    }
  }      // end of the loop over mixture components

    // Recalculate SD of the error
  mixMoments(mixMomentM, kP, wM, muM, invsigma2M, true);

  return;
}   // end of function updateMeans


