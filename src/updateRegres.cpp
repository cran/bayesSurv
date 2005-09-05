// Functions to update the regression parameters,
//  i.e. the fixed effects and the means of the random effects.
//
// It is allowed to update some parameters in blocks.
//
// Either Gibbs move or standard Metropolis-Hastings or adaptive Metropolis (AM) algorithm are used.
//   * the proposal for Metropolis-like algorithms is sampled from a mixture of (multivariate) normal and 
//     a product of univariate uniform distributions, all of them centered at a current point
//   * so that the whole block of betas is updated either using the uniform distribution
//     or using the normal distribution

// 30/01/2004: start working on it
// 09/03/2004: 'GIBBSproposalFixed' function added and Gibbs move for fixed effects allowed
//

// INPUT PARAMETERS
// 
// loglikelhood ...... log-likelihood of the observations
// randomloglik ...... log-likelihood of random effects
// randomllcl ........ pointer to an array of values of the log-likelihood of random effects, separately for each cluster
// proprandomllcl .... pointer to a proposal array of randomllcl
// doadapt ........... 0/1, 1 if this iteration is stored (determined by the thinning interval)
//                     if doadapt == 1 and adaptive algorithm is used, then quantities needed for it are
//                     recalculated
// errorTypeP......... type of the error distribution
//
#include "updateRegres.h"

using namespace std;

void
updateRegres(MHblocks* betaMH,             double** regresResM,   double** propregresResM,   
             double* loglikelhood,         double** loglikobs,    double** proploglikobs,
             double* randomloglik,         double** randomllcl,   double** proprandomllcl,
             const double* YsM,            const double* Eb0,     const double* bM,      
             const covMatrix* Dcm,         const double* XA,      double*** XXt, 
             const int* indbA,             const int* indbinXA,
             double (*logdtrans)(double),  const int* iterTotal,  const int* doadapt,
             const int* randomIntP,        const int* nrandomP,
             const int* nP,                const int* nclusterP,
             const int* errorTypeP,        const int* kP,         const int* rM,
             const double* wM,             const double* muM,     const double* invsigma2M,
             const double* tolerChol)
{
//  if (*doadapt){
//    double* temp = new double[betaMH->nBlock()];
//    for (int block = 0; block < betaMH->nBlock(); block++) temp[block] = betaMH->covpar[block][0];
//    writeToFile(temp, 1, betaMH->nBlock(), "/home/arnost/temp", "/propcov.sim", 'a');    
//    writeToFile(betaMH->meanpar, 1, betaMH->nBlock(), "/home/arnost/temp", "/means.sim", 'a');    
//    delete [] temp;
//  }

  int nParam = betaMH->nParam();

  // All priors are assumed to be normal. Include here this part of code to allow for possible future changes
  // ========================================================================================================
  double (*dprior) (double, double, double, int);
  dprior = dnorm;

  int j, accepted, block;
  double logdpriorNew, loglikNew, randomloglikNew, Paccept, u;

  for (block = 0; block < betaMH->nBlock(); block++){
    switch (betaMH->typeUpd[block]){
    case RandomWalk:
    case AdaptiveM:

      accepted = 1;

        // Propose a new value
        // ====================
      switch (betaMH->typeUpd[block]){
      case RandomWalk:             // standard MH step
        MHproposal(betaMH->proppar, betaMH->chcovpar[block], betaMH->par, betaMH->indBlock[block], 
                   &nParam, betaMH->nInBlock + block, betaMH->diagI[block], 
                   betaMH->halfRangeUnif, betaMH->weightUnif + block);      
        break;

      case AdaptiveM:             // adaptive Metropolis algorithm
	AMproposal(betaMH->proppar, betaMH->chcovpar[block], betaMH->covpar[block], betaMH->par, betaMH->indBlock[block],
                   &nParam, betaMH->nInBlock + block, betaMH->diagI[block],
                   betaMH->halfRangeUnif, betaMH->weightUnif + block,
                   betaMH->eps + block, betaMH->sdNum + block, tolerChol);
        break;
      }

        // Compute the acceptance ratio and accept/reject
        // ===============================================
          // log-prior density evaluated at the proposal
      logdpriorNew = 0.0;
      for (j = 0; j < betaMH->nInBlock[block]; j++){
        logdpriorNew += dprior(betaMH->proppar[betaMH->indBlock[block][j]], 
                            betaMH->priorMean[betaMH->indBlock[block][j]], betaMH->priorSD[betaMH->indBlock[block][j]], 1);
      }

          // needed log-posterior parts evaluated at the proposal
      loglikNew = *loglikelhood;
      randomloglikNew = *randomloglik;

      if (betaMH->nFixedB[block] > 0){     // new regresRes and new log-likelihood must be computed
        for (j = 0; j < *nP; j++) (*propregresResM)[j] = (*regresResM)[j];
        regresResidual(*propregresResM, betaMH->par, betaMH->proppar, betaMH->indBlock[block], betaMH->nInBlock + block, XA, indbA, nP);
        logLikelihood(&loglikNew, *proploglikobs, nP, *propregresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, 
                      errorTypeP, randomIntP, logdtrans);
      }

      if (betaMH->nRandomB[block] > 0){    // new value of the density of random effects must be computed
        randomLogLikelihood(&randomloglikNew, *proprandomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaMH->proppar, Dcm, Eb0, indbinXA); 
      }
    
          // log(acceptance probability) and acceptance probability 
      Paccept = loglikNew - (*loglikelhood);
      Paccept += randomloglikNew - (*randomloglik);
      Paccept += logdpriorNew - betaMH->logdprior[block];
      Paccept = exp(Paccept);

      if (Paccept < 1.0){
        u = runif(0, 1);
        if (u > Paccept){
          accepted = 0;
        }  
      }


        // Perform changes in the case of acceptance and clean in the case of rejectance
        // ==============================================================================
      if (accepted){
        betaMH->sumAccept[block]++;
        *loglikelhood = loglikNew;
        *randomloglik = randomloglikNew;
        if (betaMH->nFixedB[block] > 0){
          changePointers(regresResM, propregresResM); 
          changePointers(loglikobs, proploglikobs);
        }
        if (betaMH->nRandomB[block] > 0){
          changePointers(randomllcl, proprandomllcl); 
        }
        for (j = 0; j < betaMH->nInBlock[block]; j++){
          betaMH->par[betaMH->indBlock[block][j]] = betaMH->proppar[betaMH->indBlock[block][j]];
        }
      }
      else{
        for (j = 0; j < betaMH->nInBlock[block]; j++){
          betaMH->proppar[betaMH->indBlock[block][j]] = betaMH->par[betaMH->indBlock[block][j]];
        }
      }

        // Adapt the proposal covariance matrix and sampled mean (if needed)
        // ==================================================================
      if (betaMH->typeUpd[block] == AdaptiveM && (*doadapt)){     // (AM) adaptive step
        AMadapt(betaMH->covpar[block], betaMH->meanpar, betaMH->par, betaMH->indBlock[block], betaMH->nInBlock + block,
                betaMH->diagI[block], iterTotal, betaMH->eps + block, betaMH->sdNum + block);
      }
      break;

    case Gibbs:                // Gibbs step
        // Sample new value of beta in a given block and update also regresResM
        // =====================================================================
      if (betaMH->nFixedB[block] > 0 && betaMH->nRandomB[block] == 0){
        GIBBSproposalFixed(betaMH->par, *regresResM, betaMH->covpar[block], betaMH->chcovpar[block], 
                           betaMH->priorMean, betaMH->priorInvVar, Eb0, randomIntP,
                           XA, XXt[block], betaMH->diagI[block], wM, muM, invsigma2M, rM, 
                           betaMH->indBlock[block], &nParam, betaMH->nInBlock + block, nP, errorTypeP, tolerChol);
      }
      else{
        if (betaMH->nRandomB[block] > 0 && betaMH->nFixedB[block] == 0){
          GIBBSproposalMeanRandom(betaMH->par, betaMH->covpar[block], betaMH->chcovpar[block], 
                                  betaMH->priorMean, betaMH->priorInvVar, Eb0, bM, Dcm, betaMH->diagI[block], indbA, indbinXA, 
                                  betaMH->indBlock[block], &nParam, betaMH->nInBlock + block, 
                                  nrandomP, nclusterP, errorTypeP, tolerChol);
        }
        else
          throw returnR("C++ Error(updateRegres): Gibbs move for blocks with both fixed effects and means of random effects not implemented.", 1);
      }

        // Update quantities that change
        // =============================
      betaMH->sumAccept[block]++;
      if (betaMH->nFixedB[block] > 0){
        logLikelihood(loglikelhood, *loglikobs, nP, *regresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);
      }
      if (betaMH->nRandomB[block] > 0){
        randomLogLikelihood(randomloglik, *randomllcl, &ZERO_INT, nclusterP, nclusterP, bM, betaMH->par, Dcm, Eb0, indbinXA); 
      }
      break;

    default:
      throw returnR("C++ Error(updateRegres): Unknown or unimplemented update type for regression parameters.", 1);
    }  // end of the first switch;

  }  // end of the loop over blocks

  return;
}   // end of the function updateRegres


//
// ***** GIBBSproposalFixed *****
//
// Update of a sub-vector of regression parameters, all corresponding to fixed effects
//  using a Gibbs move
//
// Update also regression residuals
//
// betaM ........ vector of all regression parameters
//       on OUTPUT: appropriate ones are replaced by newly sampled values
// regresResM ... regression residuals for all observations
//       on OUTPUT: apropriately updated
// covpar ....... working space to compute a covariance matrix of full conditional
//       on OUTPUT: this matrix is here stored
// ichicovpar ... working space to compute an inverse of the Cholesky decomposition of the inverse covariance matrix of the full conditional
//       on OUTPUT: this matrix is here stored
// priormean .... vector of prior means for all betas
// priorinvvar .. vector of prior inverse-variances for all betas
// XA ........... full design matrix (for all observations)
// XXtb ......... array of lower triangles of matrices x*x' for given block and for all observations
// diagIXXtb .... indeces of diagonal elements in XXtb 
// indUpd ....... indeces of betas that are to be updated
// npar ......... length of vector betaM
// nupdate ...... number of betas that are to be updated
// 
void
GIBBSproposalFixed(double* betaM,            double* regresResM,
                   double* covpar,           double* ichicovpar,
                   const double* priormean,  const double* priorinvvar,
                   const double* Eb0,        const int* randomIntP,
                   const double* XA,         double** XXtb,              const int* diagIXXtb,
                   const double* wM,         const double* muM,          const double* invsigma2M,  const int* rM,
                   const int* indUpd,        const int* npar,            const int* nupdate,        const int* nP,
                   const int* errorTypeP,    const double* tolerChol)
{

  double* propMean;
  double* propMeanTemp;
  double intcptadd;  
  int lcov;
  int i, j, obs, rank;
  double helpd;

  switch (*errorTypeP){
  case Mixture:
    lcov = (*nupdate * (*nupdate + 1)) / 2;

    intcptadd = (*randomIntP ? (*Eb0) : 0.0);

    propMean = new double[*nupdate];
    propMeanTemp = new double[*nupdate];

      // Mean of full conditional distribution, part 1 (Psi^{-1}*nu) AND
      // Inverse variance of full conditional distribution, part 1 (Psi^{-1}) (store it in covpar)
    for (j = 0; j < *nupdate; j++){
      propMeanTemp[j] = priorinvvar[indUpd[j]] * priormean[indUpd[j]];
      covpar[diagIXXtb[j]] = priorinvvar[indUpd[j]];
      for (i = j + 1; i < *nupdate; i++){
        covpar[diagIXXtb[j] + i - j] = 0.0;
      }
    }

      // Loop over observations
    for (obs = 0; obs < *nP; obs++){

        // Compute y - (mu - Eb0) - x(-M)'beta(-M) - z'b, store it in regresResM
      for (j = 0; j < *nupdate; j++){
        regresResM[obs] += XA[indUpd[j]*(*nP) + obs] * betaM[indUpd[j]];
      }
      regresResM[obs] -= (muM[rM[obs]] - intcptadd);
  
        // Inverse variance of full conditional distribution, part 2 (+ 1/sigma2 * x(M) * x(M)')
      for (j = 0; j < lcov; j++){
        covpar[j] += invsigma2M[rM[obs]] * XXtb[obs][j];
      }

        // Mean of full conditional distribution, part 2 (+ 1/sigma2 * x(M) * v)
      helpd = invsigma2M[rM[obs]] * regresResM[obs];
      for (j = 0; j < *nupdate; j++){
        propMeanTemp[j] += helpd * XA[indUpd[j]*(*nP) + obs];
      }
    }  // end of the loop over observations

      // Cholesky decomposition of the inverse variance of full conditional distrib.
    cholesky(covpar, &rank, nupdate, diagIXXtb, tolerChol);
 
      // Variance of the full conditional distribution
      //  and the inverse of the Cholesky decomposition of the inverse variance
    chinv2(covpar, ichicovpar, nupdate, diagIXXtb);

      // Mean of full conditional distribution, part 3 (* var(beta(M)|...))
    Mxa(propMean, propMeanTemp, covpar, &ZERO_INT, nupdate, nupdate, diagIXXtb);

      // Sample
    rmvtnorm2(betaM, propMean, ichicovpar, &ZERO_INT, indUpd, npar, nupdate, nupdate, &ONE_INT, diagIXXtb, &ZERO_INT);

      // Update regresResM
    for (obs = 0; obs < *nP; obs++){
      for (j = 0; j < *nupdate; j++){
        regresResM[obs] -= XA[indUpd[j]*(*nP) + obs] * betaM[indUpd[j]];
      }
      regresResM[obs] += (muM[rM[obs]] - intcptadd);
    }

    delete [] propMean;
    delete [] propMeanTemp;
    break;

  case Spline:
    throw returnR("C++ Error(GIBBSproposalFixed): Gibbs update of fixed effects has not been implemented for a spline version.", 99);
    break;

  case PolyaTree:
    throw returnR("C++ Error(GIBBSproposalFixed): Gibbs update of fixed effects has not been implemented for a Polya tree version.", 99);
    break;

  default:
    throw returnR("C++ Error(GIBBSproposalFixed): Unknown errorType appeared in a call to GIBBSproposalFixed function.", 99);
  }    // end of switch (*errorTypeP)

  return;
}    // end of the function GIBBSproposalFixed


//
// ***** GIBBSproposalMeanRandom *****
//
// Update of a sub-vector of regression parameters, all corresponding to means of random effects
//  using a Gibbs move
//
// REMARK: both spline and mixture versions have same Gibbs move in this case
//
// betaM ........ vector of all regression parameters
//       on OUTPUT: appropriate ones are replaced by newly sampled values
// covpar ....... working space to compute a covariance matrix of full conditional
//       on OUTPUT: this matrix is here stored
// ichicovpar ... working space to compute an inverse of the Cholesky decomposition of the inverse covariance matrix of the full conditional
//       on OUTPUT: this matrix is here stored
// priormean .... vector of prior means for all betas
// priorinvvar .. vector of prior inverse-variances for all betas
// Eb0 .......... mean of the random intercept
// bM ........... values of random effects (all)
// Dcm .......... covariance matrix of random effects
// diagIXXtb .... indeces of diagonal elements in XXtb for given block (it is used to work with covpar and ichicovpar)
// indbA ........ indeces of betas to which random effects correspond
// indbinXA ..... indeces of random effects to which betas correspond
// indUpd ....... indeces of betas that are to be updated
// npar ......... length of vector betaM
// nupdate ...... number of betas (gammas) that are to be updated
// nrandomP ..... total number of random effects
// nclusterP .... number of clusters
// errorTypeP ... type of the error distribution
// tolerChol .... tolerance for the Cholesky decomposition
// 
void
GIBBSproposalMeanRandom(double* betaM,
                        double* covpar,           double* ichicovpar,
                        const double* priormean,  const double* priorinvvar,
                        const double* Eb0,        const double* bM,
                        const covMatrix* Dcm,     const int* diagIXXtb,
                        const int* indbA,         const int* indbinXA,
                        const int* indUpd,        const int* npar,            
                        const int* nupdate,       const int* nrandomP,         const int* nclusterP,
                        const int* errorTypeP,    const double* tolerChol)
{
  double* propMean;
  double* propMeanTemp;
  double* sumgammab;
  double* sumbM;
  int* indRandomKeep;
  int* indRandomUpdate;
  int i, j, ii, jj, cl, rank;
  bool match;
  double helpd;

  switch (*errorTypeP){
  case Mixture:
  case Spline:
    propMean = new double[*nupdate];
    propMeanTemp = new double[*nupdate];
    sumgammab = new double[*nrandomP - *nupdate];
    sumbM = new double[*nupdate];
    indRandomKeep = new int[*nrandomP - *nupdate];
    indRandomUpdate = new int[*nupdate];

      // Indeces of random effects whose means are not updated and whose means are updated
      // (these indeces correspond to the covariance matrix D)
    ii = 0;
    for (j = 0; j < *nrandomP; j++){
      match = false;
      for (i = 0; i < *nupdate && (!match); i++){
        if (indbinXA[j] == indUpd[i]) match = true;
      }
      if (!match){
        indRandomKeep[ii] = j;
        ii++;
      }
    }
    for (i = 0; i < *nupdate; i++){
      indRandomUpdate[i] = indbA[indUpd[i]];
    }

      // Inverse variance of full conditional distribution (Psi^{-1} + N*V_M) (store it in covpar) AND
      // mean of the full conditional distribution, part 1 (Psi^{-1}*nu)
    for (j = 0; j < *nupdate; j++){
      jj = indbA[indUpd[j]];
      if (jj < 0) throw returnR("C++ Error(GIBBSproposalMeanRandom): programming error somewhere, contact the author.", 99);
      covpar[diagIXXtb[j]] = priorinvvar[indUpd[j]] + (*nclusterP) * (Dcm->icovm[Dcm->diagI[jj]]);
      for (i = j + 1; i < *nupdate; i++){
        ii = indbA[indUpd[i]];
        if (ii > jj) covpar[diagIXXtb[j] + i - j] = (*nclusterP) * (Dcm->icovm[Dcm->diagI[jj] + ii - jj]);
        else         covpar[diagIXXtb[j] + i - j] = (*nclusterP) * (Dcm->icovm[Dcm->diagI[ii] + jj - ii]);
      }
      propMeanTemp[j] = priorinvvar[indUpd[j]] * priormean[indUpd[j]];
    }    

      // Cholesky decomposition of the inverse variance of full conditional distrib.
    cholesky(covpar, &rank, nupdate, diagIXXtb, tolerChol);
 
      // Variance of the full conditional distribution
      //  and the inverse of the Cholesky decomposition of the inverse variance
    chinv2(covpar, ichicovpar, nupdate, diagIXXtb);

      // Mean of the full conditional distribution, part 2 (+ V_M*\sum b_M - W*\sum(gamma_{-M} - b_{-M}))
      //  a) \sum b_M (store it in sumbM)
      //  b) += V_M * \sum b_M (store it first in propMean)
      //  c) \sum (gamma_{-M} - b_{-M}) (store it in sumgammab)
      //  d) -= W * \sum(gamma_{-M} - b_{-M}) (store it first in propMean)
    for (j = 0; j < *nupdate; j++){
      sumbM[j] = 0.0;
      for(cl = 0; cl < *nclusterP; cl++){
        sumbM[j] += bM[cl*(*nrandomP) + indbA[indUpd[j]]];
      }
    }
    Mxa2(propMean, sumbM, Dcm->icovm, indRandomUpdate, nupdate, nrandomP, Dcm->diagI);    
    for (j = 0; j < *nupdate; j++) propMeanTemp[j] += propMean[j];
    jj = *nrandomP - *nupdate;
    if (jj > 0){
      for (j = 0; j < jj; j++){
        sumgammab[j] = 0.0;
        helpd = (indbinXA[indRandomKeep[j]] < 0) ? (*Eb0) : betaM[indbinXA[indRandomKeep[j]]];
        for(cl = 0; cl < *nclusterP; cl++){
          sumgammab[j] += (helpd - bM[cl*(*nrandomP) + indRandomKeep[j]]); 
        }
      }
      Wxa(propMean, sumgammab, Dcm->icovm, indRandomUpdate, indRandomKeep, &jj, nrandomP, nupdate, Dcm->diagI);
      for (j = 0; j < *nupdate; j++) propMeanTemp[j] -= propMean[j];
    }

      // Mean of full conditional distribution, part 3 (* var(gamma(M)|...))
    Mxa(propMean, propMeanTemp, covpar, &ZERO_INT, nupdate, nupdate, diagIXXtb);

      // Sample
    rmvtnorm2(betaM, propMean, ichicovpar, &ZERO_INT, indUpd, npar, nupdate, nupdate, &ONE_INT, diagIXXtb, &ZERO_INT);

    delete [] propMean;
    delete [] propMeanTemp;
    delete [] sumbM;
    delete [] sumgammab;
    delete [] indRandomKeep;
    delete [] indRandomUpdate;
    break;

  case PolyaTree:
    throw returnR("C++ Error(GIBBSproposalMeanRandom): Gibbs update of means of random effects has not been implemented for a Polya tree version.", 99);
    break;

  default:
    throw returnR("C++ Error(GIBBSproposalMeanRandom): Unknown errorType appeared in a call to GIBBSproposalFixed function.", 99);

  }    // end of switch (*errorTypeP)

  return;
}    // end of the function GIBBSproposalMeanRandom               
    



