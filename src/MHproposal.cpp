// Functions to propose a value for the Metropolis-Hastings algorithm
//  * different methods are used

// 28/01/2004: start working on it

#include "MHproposal.h"

using namespace std;


// ==========================================================================
// *** AMproposal ***
// ======================
// Proposal done using an adaptive Metropolis algorithm
//   ref. Haario, Saksman, Tamminen (2001), Bernoulli, vol. 7, pp. 223-242
//
// AM proposal is sampled from a (multivariate) normal distribution,
//   centered at a current value
//   with a covariance matrix equal to the empirical covariance matrix of all up to now
//   sampled values
// It is further mixed with a uniform proposal, also centered at a current value
//
// OUTPUT PARAMETERS:
//
// proppar ........ a vector of parameters where desired components are replaced by the proposal
//                  -> array of length npar
// chcovpar ....... Cholesky decomposition of last used proposal covariance matrix
//                  or of a positive definite matrix which is as close as possible to the covariance matrix
//                  (only a lower triangle)
//                   -> array of length (nupdate * (nupdate + 1))/2
//
// INPUT PARAMETERS:
//
// covpar ......... updated covariance matrix of up to now sampled parameters
//                  (<=> covariance matrix for desired components only)
//                  (only a lower triangle is stored)
//                   -> array of length (nupdate * (nupdate + 1))/2
// par ............ a vector of current parameters values
//                  -> array of length npar
// indUpd ......... indeces of parameters that are to be updated 
//                  -> array of length nupdate
// npar ........... length of the arrays 'proppar', 'meanpar', 'par', 'halfRangeUnif'
// nupdate ........ length of the array indUpd,
//                  number of parameters that are to be updated,
//                  it determines the dimension of the matrices 'covpar' and 'icovpar'
// diagI .......... indeces of diagonal elements of the covariance matrices in the storage arras
//                  -> array of length nupdate
// halfRangeUnif .. specification of the range for the uniform part of the proposal,       
//                  which is always Unif(current beta +- halfRangeUnif)
//                  -> array of length npar
// weightUnif ..... weight of the uniform parts of the proposal in the mixture            
//                  -> scalar
// eps ............ epsilon added always to the diagonal of the proposal covariance matrix 
//                  -> scalar
// sdNum .......... an s[d] number used in the (AM) algorithm                   
//                  -> scalar
// tolerChol ...... tolerance for the Cholesky decomposition
//
void
AMproposal(double* proppar,              double* chcovpar,     
           const double* covpar,         const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif,   
           const double* eps,            const double* sdNum,        const double* tolerChol)
{
  double helpd;

    // Cholesky decomposition of matrix covpar (try to make it positive definite if it's not)
    // (some initial matrix is used if *iter = 1)
    // store it in icovpar
    // ====================================================
  helpd = (*sdNum)*(*eps);        // what to add to the diagonal in the case of non-positive definitness
  int rank = 0;
  int attempt = 0;
  chposDef(covpar, chcovpar, &rank, &attempt, nupdate, diagI, tolerChol, &helpd);
  if (rank < *nupdate)     // I do not expect that this ever happens
    throw returnR("C++ Error: (AM) algorithm failed for one of the blocks.", 99);


  // Give the proposal
  // ==================
    // Decide whether the proposal will be sampled from normal or uniform
    //   and sample from the appropriate distribution
  helpd = runif(0, 1);
  if (helpd >= *weightUnif) rmvtnorm(proppar, par, chcovpar, indUpd, indUpd, npar, npar, nupdate, &ONE_INT, diagI, &ZERO_INT);
  else                      rmvtiunif(proppar, par, halfRangeUnif, indUpd, indUpd, npar, npar, nupdate, &ONE_INT, &ZERO_INT);

  return;
}    // end of function AMproposal


// ========================================================================================
// *** AMadapt ***
// =========================
// Function to adapt a covariance matrix for the proposal to be used in the next iteration
//   and the mean of up to now sampled values
//
// meanpar ........ mean of up to now sampled values with desired components updated
//                   -> array of length npar
//
// iter ......... number of the iteration 
//
void
AMadapt(double* covpar,     double* meanpar,     const double* par,  
        const int* indUpd,  const int* nupdate,  const int* diagI,   const int* iter,
        const double* eps,  const double* sdNum)  
{

  int i, j, k;
  double helpd;
  double* oldmean = new double[*nupdate];

      // Update the mean 
  for (i = 0; i < *nupdate; i++){
    oldmean[i] = meanpar[indUpd[i]];               // store here now the old mean
    meanpar[indUpd[i]] *= (*iter);
    meanpar[indUpd[i]] += par[indUpd[i]];
    meanpar[indUpd[i]] /= (*iter + 1);         
  }

      // Update the covariance matrix of the proposal
  for (j = 0; j < *nupdate; j++){                   // over columns
    for (i = j; i < *nupdate; i++){                 // over rows
      k = diagI[j] + i - j;                         // index of the element in the array        
      covpar[k] *= (double(*iter - 1)/double(*iter));
      helpd = (*iter) * oldmean[i] * oldmean[j] - (*iter + 1) * meanpar[indUpd[i]] * meanpar[indUpd[j]]
              + par[indUpd[i]] * par[indUpd[j]] + (i == j ? *eps : 0.0);
      covpar[k] += ((*sdNum)/(*iter)) * helpd;
    }
  }

  delete [] oldmean;
  return;
}


// ==========================================================================
// *** MHproposal ***
// ======================
//
// PARAMETERS different from the previous function:
//
// chcovpar ...... Cholesky decomposition of the proposal covariance matrix for
//                 updated components of the vector par stored as a lower triangle
//                 in an array
void
MHproposal(double* proppar,              
           const double* chcovpar,       const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif)
{
    // Decide whether the proposal will be sampled from normal or uniform
    //   and sample from the appropriate distribution
  double helpd;
  helpd = runif(0, 1);
  if (helpd >= *weightUnif) rmvtnorm(proppar, par, chcovpar, indUpd, indUpd, npar, npar, nupdate, &ONE_INT, diagI, &ZERO_INT);
  else                      rmvtiunif(proppar, par, halfRangeUnif, indUpd, indUpd, npar, npar, nupdate, &ONE_INT, &ZERO_INT);

  return;

}    // end of function MHproposal
