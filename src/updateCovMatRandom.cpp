// Function to update a covariance matrix of random effects
//
// * a Gibbs move is used
//
// 12/02/2004: start working on it
// 13/03/2004: allow a general mean of the random intercept
// 02/08/2004: uniform prior for std. deviation of univariate random effects added
//
#include "bayessurvreg.h"

using namespace std;

//
// ***** updateCovMatRandom *****
//
// priorD ......... indicator of the prior distribution for matrix D
//                  0 = InvWishart, 1 = SDUniform
// priordf ........ prior degrees of freedom for inverse-Wishart if priorD = InvWishart
//                  1/B^2 if priorD = SDUniform, i.e. Unif(0, B) prior for std. dev.(b)
// priorscaleMat .. prior scale matrix for inverse-Wishart if priorD = InvWishart
//                  B if priorD = SDUniform, i.e. Unif(0, B) prior for std. dev.(b)
//
void
updateCovMatRandom(covMatrix* Dcm,           const double* bM,             
                   const double* betaM,      const double* Eb0,
                   const int* priorD,        const double* priordf,   const double* priorscaleMat, 
                   const int* indbinXA,      const int* nclusterP,    const int* nP,      
                   const double* tolerChol,  const double* tolerQR)
{
  int i, rank;

  int nRandom = Dcm->nrow();
  int lSS = Dcm->larray();
  double df, shape, rate;

    // Compute sum of squares of random effects, store it in Dcm->covm
  sumSquares(Dcm->covm, bM, betaM, Eb0, indbinXA, Dcm->diagI, nclusterP, &nRandom, &lSS);

  switch (*priorD){
  case InvWishart:
      // Degrees of freedom for a full conditional
    df = *priordf + (*nclusterP);

      // Compute a scale matrix of the inverse-Wishart distribution which is a full conditional,
      //   store it in Dcm->covm
    for (i = 0; i < lSS; i++){
      Dcm->covm[i] += priorscaleMat[i];
    } 

      // Cholesky decomposition of a scale matrix
      //   and an inverse of this Cholesky decomposition
    cholesky(Dcm->covm, &rank, &nRandom, Dcm->diagI, tolerChol);
    if (rank < nRandom) throw returnR("C++ Error: Scale matrix for update of D is not positive definite.", 1);
    chinv(Dcm->covm, &nRandom, Dcm->diagI, &ONE_INT);

      // Sample new inverse D matrix
    rwishart2(Dcm->icovm, &nRandom, &df, Dcm->covm, Dcm->diagI, &ONE_INT, &ZERO_INT);
    break;

  case SDUniform:
      // Shape and rate of the full conditional distribution
    shape = 0.5*((*nclusterP) - 1);  
    rate = 0.5*(Dcm->covm[0]);
    
      // Sample new inverse var(b), i.e. new inverse D matrix
    rltruncGamma(Dcm->icovm, &shape, &rate, priordf, &ONE_INT, &ZERO_INT);     
    break;
  }

    // Compute D, and iLi where iLi = L^{-1} and L*L' = D^{-1}
  Dcm->updateFromInv(tolerChol, tolerQR);
  return;
}    // end of the function updateCovMatRandom


//
// ***** sumSquares *****
//
// Compute a vectorized sum of squares
// i.e. sum (b_i - gamma)*(b_i - gamma)'
//
// Return its lower triangle
//
//
// lSS ....... length of the array sumSq, it should be equal to (*nrandomP * (*nrandomP + 1))/2
//
void
sumSquares(double* sumSq,          const double* bM,       
           const double* betaM,    const double* Eb0,
           const int* indbinXA,    const int* diagI,       
           const int* nclusterP,   const int* nrandomP,    const int* lSS)
{
  int i, j, cl;
  double* bMinGamma = new double[*nrandomP];

  for (i = 0; i < *lSS; i++){
    sumSq[i] = 0.0;
  }

  for (cl = 0; cl < *nclusterP; cl++){
    for (j = 0; j < *nrandomP; j++){
      bMinGamma[j] = bM[(*nrandomP)*cl + j] - (indbinXA[j] < 0 ? (*Eb0) : betaM[indbinXA[j]]);
    }
    for (j = 0; j < *nrandomP; j++){
      for (i = j; i < *nrandomP; i++){
        sumSq[diagI[j] + i - j] += bMinGamma[i] * bMinGamma[j];
      }
    }
  }

  delete [] bMinGamma;
  return;
}    // end of function sumSquares
