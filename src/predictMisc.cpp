// Various functions used to compute either predictive quantities or other
//   quantities from the resulting Markov chain

//
// 07/03/2004: function 'predictSurv' 
// 18/03/2004: function 'predictData' 
// 24/03/2004: 'predictSurv' rewritten using 3-dimensional arrays 
// 25/03/2004: function 'predictRandom' 
// 26/03/2004: function 'predictET'
//

#include "bayessurvreg.h"

using namespace std;

// =====================================================================
// ******* predictSurv *******
//
// Function to compute  
//  * predictive survivor functions
//  * predictive hazards and cumulative hazards
//  * prediction for 
//  for various combination of covariates
//  (these combinations are implicitely given by the value of regresResM)
// Survivor functions etc. may be evaluated in different grids for different observations
//
//  * marginalization over the vector of random effects is done via Monte Carlo integration
//    (i.e. by using sampled values of b)
//
// PARAMETERS:
//
// SM[nP][ngridM[obs]][niter] .......... predicted survivor functions evaluated in a grid of values 
//                                       given by gridM and loggridM
// lambdaM[nP][ngridM[obs]][niter] ..... predicted hazard functions evaluated in a grid of values
// LambdaM[nP][ngridM[obs]][niter] ..... predictive cumulative hazard function evaluated in a grid of values
//
// iter ........ index of current iteration 
// gridM ....... grids of values where SM and other quantities should be evaluated
//                     (nP arrays, ith array of length ngridM[i]
//                * this function does not change it
// loggridM .... logarithm of grids of values where SM and other quantities should be evaluated
//                     (nP arrays, ith array of length ngridM[i]
// time0P ...... starting time of the survival model
//                * this function does not change it
// regresPredM . regression predictors (x'beta + z'b) (all)                   (nP x 1)
// rM .......... component pertinences (if you do not not integrate over them)
// wM .......... component weights                                            (kP x 1)
//               (used if errorType == Mixture, Spline)
// muM ......... component means                                              (kP x 1)
// sigmaM ...... component standard deviations                                (kP x 1)
// Eb0 ......... mean of the random intercept
// kP .......... number of mixture components                                 (1)
// nP .......... number of observations (all)                                 (1)
// ngridM ...... lengths of grids for each predictive survivor function       (nP x 1)
// errorTypeP .. type of the error density 
//               (0 = Mixture, 1 = Spline, 2 = PolyaTree)
// randomIntP .. 0/1 indicating whether a random intercept is in the model
// hazardP ..... 0/1, indicator whether the hazard functions are to be computed
// cumhazardP .. 0/1, indicator whether the cumulative hazard functions are to be computed
//
void
predictSurv(double*** SM,              double*** lambdaM,         double*** LambdaM,
            const int iter,            double** gridM,            double** loggridM,  const double* time0P,
            const double* regresPredM,
            const int* rM,             const double* wM,          const double* muM,  const double* sigmaM,
            const double* Eb0,         const int* kP,             const int* nP,      const int* ngridM,
            const int* errorTypeP,     const int* randomIntP,
            const int* hazardP,        const int* cumhazardP)
{
  int obs, i, j;
  double logt_eta;      // log(t) - beta'x - b'z
  double semidensity;   // density not divided by t
  double intcptadd;

  switch (*errorTypeP){
  case Mixture:
  case Spline:
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      for (i = 0; i < ngridM[obs]; i++){
        logt_eta = loggridM[obs][i] - regresPredM[obs];
        SM[obs][i][iter] = 0.0;
        for (j = 0; j < *kP; j++){
          if (wM[j] <= 0.0) continue;
          SM[obs][i][iter] += wM[j] * pnorm(logt_eta, muM[j] - intcptadd, sigmaM[j], 0, 0);
        }
        if (*hazardP){
          semidensity = 0.0;
          for (j = 0; j < *kP; j++){
            if (wM[j] <= 0.0) continue;
            semidensity += wM[j] * dnorm(logt_eta, muM[j] - intcptadd, sigmaM[j], 0);
          }
          if (SM[obs][i][iter] <= 0.0) lambdaM[obs][i][iter] = FLT_MAX;
          else                         lambdaM[obs][i][iter] = (1 / (gridM[obs][i] - (*time0P))) * (semidensity / SM[obs][i][iter]);
          if (lambdaM[obs][i][iter] >= FLT_MAX) lambdaM[obs][i][iter] = FLT_MAX;
        }
        if (*cumhazardP){
          if (SM[obs][i][iter] <= 0.0) LambdaM[obs][i][iter] = FLT_MAX;
          else                         LambdaM[obs][i][iter] = -log(SM[obs][i][iter]);
          if (LambdaM[obs][i][iter] <= 0) LambdaM[obs][i][iter] *= (-1);           // this happens if SM is 1
          if (LambdaM[obs][i][iter] >= FLT_MAX) LambdaM[obs][i][iter] = FLT_MAX;
        }
      }    // end of the loop over grid values
    }    // end of the loop over observations
    break;

  case WhoKnows:    // version when rM is given and it is not averaged over it
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      for (i = 0; i < ngridM[obs]; i++){
        logt_eta = loggridM[obs][i] - regresPredM[obs];
	SM[obs][i][iter] = pnorm(logt_eta, muM[rM[obs]] - intcptadd, sigmaM[rM[obs]], 0, 0);
        if (*hazardP){
          semidensity = dnorm(logt_eta, muM[rM[obs]] - intcptadd, sigmaM[rM[obs]], 0);
          if (SM[obs][i][iter] <= 0.0) lambdaM[obs][i][iter] = FLT_MAX;
          else                         lambdaM[obs][i][iter] = (1 / (gridM[obs][i] - (*time0P))) * (semidensity / SM[obs][i][iter]);
          if (lambdaM[obs][i][iter] >= FLT_MAX) lambdaM[obs][i][iter] = FLT_MAX;
        }
        if (*cumhazardP){
          LambdaM[obs][i][iter] = -pnorm(logt_eta, muM[rM[obs]] - intcptadd, sigmaM[rM[obs]], 0, 1);
          if (LambdaM[obs][i][iter] >= FLT_MAX) LambdaM[obs][i][iter] = FLT_MAX;
        }
      }    // end of the loop over grid values
    }    // end of the loop over observations
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function predictSurv called.", 1);
    break;
      
  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function predictSurv.", 1);

  }

  return;
}    // end of the function predictSurv


// ===================================================================================================
// ****** predictData *******
//
// Function used to compute predictive event times
// * this corresponds to sampling from a full conditional
//   given T > 0
//
// PARAMETERS:
// 
// YsM ......... on INPUT: current vector of sampled log(event times)                  (nP x 1)
//               on OUTPUT: newly sampled log(event times)
// INPUT PARAMETERS:
//
// regresPredM .. current vector of regression predictors (x'beta - z'b)     (nP x 1)
// rM .......... vector of pertainence indicators         (nP x 1)
//               (when *errorTypeP = Spline, Mixture, this is a working space)
// cumwM ....... vector of cumulative component weights   (at least kP x 1)
//               (it is not 'const' but it is not modified by this function)
// muM ......... vector of component means                (at least kP x 1)
// sigmaM ...... vector of component standard deviations  (at least kP x 1)
// Eb0 ......... mean of the random intercept
// kP .......... current number of mixture components     (1)
// nP .......... number of predictive event times that are sampled
// errorTypeP . indicator whether the error density is assumed to be a spline 
//              (i.e. mixtureNM and rM are then ignored
//              or a classical mixture where it is assumed that each residual belongs to one mixture component
//              (i.e. cumwM is then ignored)
//              or something else
// randomIntP . 0/1 indicating whether a random intercept is in the model
//
void
predictData(double* YsM,             const double* regresPredM,
            int* rM,                 double* cumwM,             const double* muM,     const double* sigmaM,
            const double* Eb0,       const int* kP,             const int* nP,
            const int* errorTypeP,   const int* randomIntP)
{
  int obs, j;
  double intcptadd;

  switch (*errorTypeP){
  case Mixture:
  case Spline:
    discreteSampler(rM, cumwM, kP, nP, &ONE_INT, &ZERO_INT);
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      YsM[obs] = rnorm(regresPredM[obs] + muM[rM[obs]] - intcptadd, sigmaM[rM[obs]]);
    }
    break;

  case WhoKnows:    // version when rM is given in advance
    intcptadd = (*randomIntP ? (*Eb0) : 0.0);
    for (obs = 0; obs < *nP; obs++){
      YsM[obs] = rnorm(regresPredM[obs] + muM[rM[obs]] - intcptadd, sigmaM[rM[obs]]);
    }
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function predictData called.", 1);
    break;

  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function predictData.", 1);
  }

  return;
}

// ================================================================================
// ******* Y2T *******
//
// Compute survivor times from their transformation
//   (usually T = exp(Y), i.e. itrans = exp)
//
// T[*nP][iteration] ...... predictive survivor times
// Y[*nP] ................. predictive log(survivor times)
// iter ................... index of current iteration
//
void
Y2T(double** T,  const double* Y,  const double* time0P, const int iter,  const int* nP,  double (*itrans)(const double))
{
  int obs;
  for (obs = 0; obs < *nP; obs++){
    T[obs][iter] = (*itrans)(Y[obs]) + (*time0P);
  }

  return;
}


// ===================================================================================================
// ****** predictRandom *******
//
// Function to sample values of random effects given their means and variance matrix
// (means and variance matrix are given by the Markov chain)
// * mean is obtained from betaM
// * variance matrix is given by Dcm
//
void
predictRandom(double* bM,
              const double* betaM,     const double* Eb0,    const covMatrix* Dcm,
              const int* nrandomP,     const int* nclusterP,
              const int* indbinXA,     const int* indUpd)
{
  if (*nrandomP < 1) return;

  int j, cl;
    // Mean of the normal distribution to sample from
  double* gamma = new double[*nrandomP];
  for (j = 0; j < *nrandomP; j++) gamma[j] = (indbinXA[j] < 0 ? (*Eb0) : betaM[indbinXA[j]]);

    // Loop over clusters, sample
  for (cl = 0; cl < *nclusterP; cl++){
    rmvtnorm2(bM + (*nrandomP)*cl, gamma, Dcm->ichicovm, &ZERO_INT, indUpd, nrandomP, nrandomP, nrandomP, 
              &ONE_INT, Dcm->diagI, &ZERO_INT);
  }

  delete [] gamma;
  return;
}


// ===================================================================================================
// ****** predictET ******
//
// Function to compute expectation of the survivor time for given combinations of covariates
//  based on sampled Markov chain
//
// * marginalization over a vector of random effects b is done analytically
// * mean of random intercept does not have to be added to regresMean since it would have
//   to be immediately subtracted from each mixture mea
//   -> that is why, this function does not have to know *Eb0 
//
// PARAMETERS:
//
// ET[nP][niter] ....... expectation of the survivor times
// iter ................ index of current iteration
// betaM[nXP] .......... regression parameters
// wM[kP]
// muM[kP]
// sigma2M[kP]
// Dcm ................. covariance matrix of random effects
// XA[nP x nXP] ........ design matrix
// kP
// nP
// indbinXA[nrandomP] .. where to find a design for random effects
//
void
predictET(double** ET,            const double* time0P, const int iter,
          const double* betaM,    const double* wM,     const double*muM,       const double* sigma2M,
          const covMatrix* Dcm,   const double* XA,
          const int* kP,          const int* nP,        const int* nXP,         const int* indbinXA,
          const int* randomIntP,  const int* nrandomP,  const int* errorTypeP)
{

  int obs, j;
  double helpd;
  double* zvec = new double [*nrandomP];
  double zDz;

  switch (*errorTypeP){
  case Mixture:
  case Spline:
    for (obs = 0; obs < *nP; obs++){

      // Compute exp(beta'x + gamma'z + 0.5*z'Dz)
      helpd = 0;
      for (j = 0; j < *nXP; j++) helpd += XA[(*nP)*j + obs] * betaM[j];
      if (*nrandomP){
        if (*randomIntP) zvec[0] = 1;
        for (j = (*randomIntP); j < *nrandomP; j++) zvec[j] = XA[(*nP)*indbinXA[j] + obs];
        axMxa(&zDz, zvec, Dcm->covm, &ZERO_INT, nrandomP, nrandomP, Dcm->diagI);
        helpd += (0.5 * zDz);
      }
      ET[obs][iter] = exp(helpd);

      // Compute sum_{j=1}^k w_j * exp(mu_j + 0.5*sigma^2_j)
      helpd = wM[0] * exp(muM[0] + 0.5*sigma2M[0]);
      for (j = 1; j < *kP; j++) helpd += wM[j] * exp(muM[j] + 0.5*sigma2M[j]);

      // Compute fimal ET
      ET[obs][iter] *= helpd;
      ET[obs][iter] += (*time0P);
    }
    break;

  case PolyaTree:
    returnR("C++ Error: Not yet implemented part (PolyaTree) of function predictET called.", 1);
    break;

  default:
    returnR("C++ Error: Unknown errorType appeared in a call to function predictET.", 1);
  }

  delete [] zvec;
  return;
}
