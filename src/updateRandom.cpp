// Functions to update random effects
//
//  * either a Metropolis-Hastings step with a normal proposal is used
//  * or a Gibbs step
//
//  * in the case of (MH) step, random effects within one cluster can be divided
//    into blocks that are then separately updated
//  * in case of Gibbs step, all random effects within one cluster are updated
//    using just one full conditional distribution
//
// 04/02/2004: start working on it
// 12/02/2004: first working version
// 14/03/2004: added possibility of a general mean of a random intercept
//
#include "bayessurvreg.h"

using namespace std;

// 
// ***** updateRandom *****
// 
// PARAMETERS:
// 
// bb .................. object where random effects are stored
// regresResM .......... regression residuals                                             (*nP)
//                       (they are updated)
// propregresResM ...... working space, this function does not change pointers            (*nP)
// randomloglik ........ value of a random log-likelihood                                 (1)
//                       (it is updated)
// randomllcl .......... values of a random log-likelihood separately for each cluster    (*nclusterP)
//                       (they are updated)
// proprandomllcl ...... working space, this function does not change pointers            (*nclusterP)
// loglikelhood ........ value of a log-likelihood                                        (1)
//                       (it is updated)
// loglikobs ........... values of log-likelihood contributions, separately               (*nP)
//                       for each observation                       
//                       (they are updated)
// proploglikobs ....... working space, this function does not change pointers            (*nP)
// betaM ............... values of all beta parameters                                    (*nXP)
// Dcm ................. covariance matrix of random effects                              
// YsM ................. (augmented) log(event times)                                     (*nP)
// XA .................. design matrix                                                    (*nP x *nXP)
// ZZt ................. array of lower triangles of matrices z*z' for each observation   (array of length *nP)
//                                     (each component of this array is of length 0.5*(*nrandomP)*(1+(*nrandomP)))
// invclusteriA ........ array of lists with indeces of observations belonging            (array of length *nclusterP)
//                       to each cluster
// indbinXA ............ indeces of columns in X corresponding to random effects          (*nrandomP)
// (*logdtrans) ........ Jacobian of a data transformation
// randomIntP .......... 0/1 indicating whether a random intercept is present             (1)
// nXP ................. number of columns in the design matrix                           (1)
// nP .................. number of observations                                           (1)
// errorTypeP .......... indicating whether the error density is spline or a mixture  
//                       or something else                                                (1)
// kP .................. current number of mixture/spline components                      (1)
// rM .................. pertinences of observations to the mixture                       (*nP)
//                       (ignored if splineP = 1)
// wM .................. current mixture/spline weights                                   (*kmaxP)
// muM ................. current mixture/spline means                                     (*kmaxP)
// invsigma2M .......... current mixture/spline inverse variances                         (*kmaxP)
// tolerChol ........... tolerance for a Cholesky decomposition                           (1)
//
void
updateRandom(bblocks* bb,                   double* regresResM,    double* propregresResM,
            double* randomloglik,           double* randomllcl,    double* proprandomllcl,
            double* loglikelhood,           double* loglikobs,     double* proploglikobs,
            const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
            const double* YsM,              const double* XA,      double** ZZt,
            const List<int>* invclusteriA,  const int* indbinXA,   double (*logdtrans)(double),
            const int* randomIntP,          const int* nXP,        const int* nP,
            const int* errorTypeP,          const int* kP,         const int* rM,
            const double* wM,               const double* muM,     const double* invsigma2M,
            const double* tolerChol)
{
  int i, j, cl, block, nobs, obs;

  int nRandom = bb->nRandom();
  int nCluster = bb->nCluster();
  int nBlock = bb->nBlock();

  int accepted;
  double Paccept, u;

  double* propb = new double[nRandom];
  double proploglik, proprandomloglik, cl_loglik, propcl_loglik;
  bool minInfty;

  switch (bb->typeUpd[0]){
    case RandomWalk:          // standard MH step or adaptive Metropolis step
    case AdaptiveM:
      for (cl = 0; cl < nCluster; cl++){
        nobs = (invclusteriA[cl]).length();

	  // Update a particular block
        for (block = 0; block < nBlock; block++){
          accepted = 1;

  	    // Copy all values of b of a given cluster to propb
            // Copy the value of randomllcl of a given cluster to proprandomllcl
            // Copy the values of loglikobs of a given cluster to proploglikobs
            // Copy the values of regresResM of a given cluster to propregresResM
            // Compute the value of the log-likelihood in the cluster only
          for (j = 0; j < nRandom; j++){
            propb[j] = bb->bM[nRandom*cl + j];
          }

          proprandomllcl[cl] = randomllcl[cl];
          proprandomloglik = *randomloglik;

          for (i = 0; i < nobs; i++){
            obs = invclusteriA[cl][i];
            proploglikobs[obs] = loglikobs[obs];
            propregresResM[obs] = regresResM[obs];
          }
          proploglik = *loglikelhood;
          clusterlogLikelihood(&cl_loglik, loglikobs, &cl, invclusteriA + cl);
         
            // Propose a new value
          switch (bb->typeUpd[0]){
 	    case RandomWalk:    // standard MH step
	       MHproposal(propb, bb->chcovpar[block], bb->bM + nRandom*cl, bb->indBlock[block], 
	                  &nRandom, bb->nInBlock + block, bb->diagI[block], 
	                  bb->halfRangeUnif, bb->weightUnif + block);      
              break;

 	    case AdaptiveM:    // adaptive Metropolis step
              throw returnR("C++ Error: Adaptive Metropolis algorithm not implemented for update of random effects.", 1);
          }

            // Compute new regres. residuals for observations from a given cluster      
          regresResidual(propregresResM, bb->bM, propb, bb->indBlock[block], bb->nInBlock + block, &cl, invclusteriA + cl,
                         XA, randomIntP, indbinXA, nP, nXP, &nRandom);

       	    // Compute a new value of a random log-likelihood for a given cluster
          randomLogLikelihood(&proprandomloglik, proprandomllcl, &cl, &nCluster, propb, betaM, Dcm, Eb0, indbinXA);

  	    // Compute a new value of a log-likelihood for observations from a cluster
            // and the value of the log-likelihood in the cluster
	  logLikelihood(&proploglik, proploglikobs, invclusteriA + cl, nP, propregresResM, YsM,
                        kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);
          clusterlogLikelihood(&propcl_loglik, proploglikobs, &cl, invclusteriA + cl);

           // log(acceptance probability) and acceptance probability 
          Paccept = proprandomllcl[cl] - randomllcl[cl];
          Paccept += propcl_loglik - cl_loglik;
          Paccept = exp(Paccept);

          if (Paccept < 1.0){
            u = runif(0, 1);
            if (u > Paccept){
              accepted = 0;
            }  
          }

  	   // Perform changes if accepted
          if (accepted){
            (bb->sumAccept[nBlock*cl + block])++;
            for (i = 0; i < bb->nInBlock[block]; i++){
              j = bb->indBlock[block][i];
              bb->bM[nRandom*cl + j] = propb[j];
            }
            randomllcl[cl] = proprandomllcl[cl];
            *randomloglik = proprandomloglik;
            for (i = 0; i < nobs; i++){
              obs = invclusteriA[cl][i];
              loglikobs[obs] = proploglikobs[obs];
              regresResM[obs] = propregresResM[obs];
            }
            *loglikelhood = proploglik;
          }     // end of if (accepted)
        }    // end of the loop over blocks
      }   // end of the loop over clusters
      break;

    case Gibbs:          // Gibbs step
        // Sample new values of b and update also regresResM
      GIBBSproposalRandom(bb->bM, regresResM, bb->covpar[0], bb->chcovpar[0], betaM, Eb0, Dcm, XA, ZZt,
                          wM, muM, invsigma2M, rM, randomIntP, &nRandom, nXP, &nCluster, nP, errorTypeP,
                          invclusteriA, indbinXA, bb->indBlock[0], tolerChol);

        // Update quantities that change
      for (i = 0; i < nCluster; i++) bb->sumAccept[i] += 1;             // Gibbs move was of course accepted for all clusters  
      logLikelihood(loglikelhood, loglikobs, nP, regresResM, YsM, kP, rM, wM, muM, invsigma2M, Eb0, errorTypeP, randomIntP, logdtrans);
      randomLogLikelihood(randomloglik, randomllcl, &ZERO_INT, &nCluster, &nCluster, bb->bM, betaM, Dcm, Eb0, indbinXA);
      break;

    default:
      throw returnR("C++ Error: Unknown or unimplemented update type.", 1);
  }

  delete [] propb;
  return;
}

//
// ***** GIBBSproposalRandom *****
//
// Update of the whole vector of random effects within all clusters
//   using a Gibbs move.
//
// Update also regression residuals.
//
//
// This function provides directly the loop over clusters since it is more
//   effective (some parts of the proposal distribution are common for all clusters).
//
// 
// bM ............ random effects of all clusters 
// regresResM .... regression residuals for all observations
// covpar ........ working space of length Dcm->larray
//     on OUTPUT: a covariance matrix of the full conditional is here stored
// ichicovpar .... working space of length Dcm->larray
//     on OUTPUT: an inverse of the Cholesky decomposition of the inverse covariance matrix
//                of the full conditional is here stored
// XA ............ full design matrix (for all observations)
// ZZt ........... array of lower traingles of matrices z*z' for all observations
//
// invclusteriA .. an array of lists of indeces of observations belonging to a specific cluster
//
// indUpd ........ indeces of b's that should be updated, with this version of the function
//                 * this must be an array {0, 1, ..., *nrandomP - 1}
//                 * it will be obtained as bb->indBlock[0] 
//
void
GIBBSproposalRandom(double* bM,                     double* regresResM, 
                    double* covpar,                 double* ichicovpar,
                    const double* betaM,            const double* Eb0,     const covMatrix* Dcm,
                    const double* XA,               double** ZZt,        
                    const double* wM,               const double* muM,     const double* invsigma2M,  const int* rM,
                    const int* randomIntP,          const int* nrandomP,       
                    const int* nXP,                 const int* nclusterP,  const int* nP,             const int* errorTypeP,
                    const List<int>* invclusteriA,  const int* indbinXA,   const int* indUpd,
                    const double* tolerChol)
{
  double* Digamma;
  double* propMean;
  double* propMeanTemp;
  double intcptadd;
  int lD;
  int cl, i, j, nobs, obs, rank;
  double helpd;

  switch (*errorTypeP){
  case Mixture:
    lD = Dcm->larray();

    intcptadd = (*randomIntP ? (*Eb0) : 0.0);

    Digamma = new double[*nrandomP];
    propMean = new double[*nrandomP];
    propMeanTemp = new double[*nrandomP];

      // Mean of full conditional distribution, part 1 (D^{-1}gamma)
    if (*randomIntP && *nrandomP == 1) 
      Digamma[0] = (Dcm->icovm[0]) * (*Eb0);
    else{
      if (*randomIntP && *nrandomP > 1){     // gamma[0] = (*Eb0), but it is not contained in betaM vector
          // first column (= row) of D^{-1} times gamma (the first element of the sum is sometimes zero)
        Digamma[0] = Dcm->icovm[0] * (*Eb0);
        for (i = 1; i < *nrandomP; i++) Digamma[0] += Dcm->icovm[i] * betaM[indbinXA[i]];

          // all other columns of D^{-1} times gamma 
        for (j = 1; j < *nrandomP; j++){
          Digamma[j] = Dcm->icovm[Dcm->diagI[j]] * betaM[indbinXA[j]];
          for (i = j + 1; i < *nrandomP; i++) Digamma[j] += Dcm->icovm[Dcm->diagI[j] + i - j] * betaM[indbinXA[i]];
          for (i = 1; i < j; i++)             Digamma[j] += Dcm->icovm[Dcm->diagI[i] + j - i] * betaM[indbinXA[i]];
          Digamma[j] += Dcm->icovm[Dcm->diagI[0] + j] * (*Eb0);
        }
      }
      else{                                  // no random intercept, gamma = subvector of betaM
        Mxa(Digamma, betaM, Dcm->icovm, indbinXA, nXP, nrandomP, Dcm->diagI);
      }
    }
      // Loop over clusters
    for (cl = 0; cl < *nclusterP; cl++){
      nobs = invclusteriA[cl].length();

        // Mean of full conditional distribution, part 1 (D^{-1}gamma)
      for (j = 0; j < *nrandomP; j++) propMeanTemp[j] = Digamma[j];
 
        // Inverse variance of full conditional distribution, part 1 (D^{-1}) (store it in covpar)
      for (j = 0; j < lD; j++) covpar[j] = Dcm->icovm[j];

        // Loop over observations in a given cluster
      for (i = 0; i < nobs; i++){
        obs = invclusteriA[cl][i];

          // Compute y - mu - x'beta, store it in regresResM
        if (*randomIntP) regresResM[obs] += bM[(*nrandomP)*cl];
        for (j = *randomIntP; j < *nrandomP; j++){
          regresResM[obs] += XA[(*nP)*indbinXA[j] + obs] * bM[(*nrandomP)*cl + j];
        }
        regresResM[obs] -= (muM[rM[obs]] - intcptadd);    

          // Inverse variance of full conditional distribution, part 2 (+ 1/sigma2 * z * z')
        for (j = 0; j < lD; j++){
          covpar[j] += invsigma2M[rM[obs]] * ZZt[obs][j];
        }
          // Mean of full conditional distribution, part 2 (+ 1/sigma2 * z * v)
        helpd = invsigma2M[rM[obs]] * regresResM[obs];
        if (*randomIntP) propMeanTemp[0] += helpd;
        for (j = *randomIntP; j < *nrandomP; j++){
          propMeanTemp[j] += helpd * XA[(*nP)*indbinXA[j] + obs];
        }     

      }  // end of the loop over observations in a given cluster

        // Cholesky decomposition of the inverse variance of full conditional distrib.
      cholesky(covpar, &rank, nrandomP, Dcm->diagI, tolerChol);

        // Variance of the full conditional distribution
        //  and the inverse of the Cholesky decomp. of the inverse variance
      chinv2(covpar, ichicovpar, nrandomP, Dcm->diagI);

        // Mean of full conditional distribution, part 3 (* var(b|...))
      Mxa(propMean, propMeanTemp, covpar, &ZERO_INT, nrandomP, nrandomP, Dcm->diagI);

        // Sample
      rmvtnorm2(bM + (*nrandomP)*cl, propMean, ichicovpar, &ZERO_INT, indUpd, nrandomP, nrandomP, nrandomP, 
                &ONE_INT, Dcm->diagI, &ZERO_INT);

      //    cout << "Cluster " << cl;
      //    cout << "  propMean = " << propMean[0] << "  propVar = " << covpar[0] << "  ichicovpar" << ichicovpar[0];
      //    cout << "  sampled b = " << bM[(*nrandomP*cl)] << endl;

        // Update regresResM
      for (i = 0; i < nobs; i++){
        obs = invclusteriA[cl][i];
        if (*randomIntP) regresResM[obs] -= bM[(*nrandomP)*cl];
        for (j = *randomIntP; j < *nrandomP; j++){
          regresResM[obs] -= XA[(*nP)*indbinXA[j] + obs] * bM[(*nrandomP)*cl + j];
        }
        regresResM[obs] += (muM[rM[obs]] - intcptadd);    
      }

    }   // end of the loop over clusters

    delete [] propMean;
    delete [] propMeanTemp;
    delete [] Digamma;
    break;

  case Spline:
    throw returnR("C++ Error: Gibbs update of random effects has not been implemented for a spline version.", 99);
    break;

  case PolyaTree:
    throw returnR("C++ Error: Gibbs update of random effects has not been implemented for a Polya tree version.", 99);
    break;

  default:
    throw returnR("C++ Error: Unknown errorType appeared in a call to GIBBSproposalRandom function.", 99);
  }    // end of switch (*errorTypeP)
  
  return;
}
