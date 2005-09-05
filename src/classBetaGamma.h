// Header for classBetaGamma.cpp
//
#ifndef _CLASS_BETA_GAMMA_H_
#define _CLASS_BETA_GAMMA_H_

#include <R.h>
#include <Rmath.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <cmath>

#include "AK_Error.h"
#include "constants.h"
#include "Gspline.h"
#include "classCovMatrix.h"
#include "cholesky.h"
#include "mvtdist.h"

const double _toler_chol_BetaGamma = 1e-10;

class BetaGamma
{
  protected:
    /*** Dimensions ***/
    int _nbeta;                // [1] number of all regression parameters
    int _nFixed;               // [1] number of fixed effects
    int _ngamma;               // [1] number of means of random effects that are not fixed (i.e. without possible random intercept)
    int _randomIntcpt;         // [1] 0/1 is random intercept?
    int _nRandom;              // [1] _ngamma + 1*(is.random.intercept)

    /*** Values of regression parameters ***/
    double* _beta;             // [_nbeta] regression parameters (both beta and gamma)

    /*** Indexing arrays ***/
    int* _indbA;               // [_nbeta] indicators of fixed and random effects
                               //           _indbA[i] = -1 iff _beta[i] is a fixed effect
                               //           _indbA[i] = j  iff _beta[i] is a mean of the jth random effect
                               //
                               //  i.e. _indbA[i] \in {-1, 1,2,...,nrandom-1} if there is random intercept in the model
                               //       _indbA[i] \in {-1, 0,1,2,...,nrandom-1} if there is NO random intercept in the model
    int* _indFixed;            // [_nFixed] indices of fixed effects in the total vector _beta
    int* _indgamma;            // [_ngamma] indices of gamma (means of random effects except random intcpt) parameters 
                               //           in the total vector _beta
    int* _indbinXA;            // [_nRandom] indeces of columns of X which correspond to a given random effect, -1 for random intercept

    /*** Prior parameters ***/
    double* _priorMean;        // [_nbeta] prior means for each parameter
    double* _priorSD;          // [_nbeta] prior standard deviations for each parameter 
    double* _priorInvVar;      // [_nbeta] with prior inverse-variances for each parameter

    /*** Stuff for update of fixed effects ***/
    int _lcovFixed;             // [1]          number of unique elements in _covFixed and _chcovFixed;
    double* _meanFixed;         // [_nFixed]    working place for the mean of the full conditional for fixed effects
    double* _meanFixedTemp;     // [_nFixed]    working place to compute the mean of the full conditional for fixed effects
    double* _covFixed;          // [_lcovFixed] working place for covariance matrix of the full conditional for fixed effects
    double* _ichicovFixed;      // [_lcovFixed] working place for Cholesky decomposition of the covariance matrix 
                                //              of the full conditional for fixed eff. and its inversions
    int* _diagIFixed;           // [_nFixed]    diagonal indeces of covparFixed and chcovparFixed

    /*** Stuff for update of means of random effects ***/
    int _lcovgamma;             // [1]           number of unique elements in _covgamma and _chcovgamma 
    double* _meangamma;         // [_ngamma]     working place for the mean of the full conditional for means of random effects
    double* _meangammaTemp;     // [_ngamma]     working place to compute the mean of the full conditional for means of random eff.
    double* _covgamma;          // [_lcovgamma]  working place for covariance matrix of the full conditional for means of random effects
    double* _ichicovgamma;      // [_lcovgamma]  working place for Cholesky decomposition of the covariance matrix 
                                //               of the full conditional for means of random effects and its inversion
    int* _diagIgamma;           // [_ngamma]    diagonal indeces of covpargamma and chcovpargamma
    double* _sumbM;             // [_ngamma]             sum of b's that correspond to gamma's
    double* _sumgammab;         // [_nRandom - _ngamma]  sum of b's that do not correspond to gamma's (here, only random intercepts)
    int* _indRandomUpdate;      // [_ngamma]             indices of random effects whose means are updated
    int* _indRandomKeep;        // [_nRandom - _ngamma]  indices of random effect whose means are not updates (here, either 0 or nothing)

 public:
   BetaGamma();
   BetaGamma(const int* parmI,  const double* parmD);
   BetaGamma(const BetaGamma& bg);
   BetaGamma& operator=(const BetaGamma& bg);
   ~BetaGamma();

   void BetaGamma2initArray(int* parmI,  double* parmD) const;

   inline int
   nbeta() const { return _nbeta;}

   inline int
   nFixed() const { return _nFixed;}

   inline int
   ngamma() const { return _ngamma;}

   inline int
   randomIntcpt() const { return _randomIntcpt;}

   inline int
   nRandom() const { return _nRandom;}

   inline double*
   betaP() const { return _beta;}

   inline double
   beta(const int& i) const
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::beta(i).", 1);
    return _beta[i];
   }

   inline void
   new_beta(const int& i, const double& newbeta)
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::new_beta(i, newbeta).", 1);
    _beta[i] = newbeta;
    return;
   }

   inline double
   priorMean(const int& i) const
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::priorMean(i).", 1);
    return _priorMean[i];
   }

   inline double
   priorSD(const int& i) const
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::priorSD(i).", 1);
    return _priorSD[i];
   }

   inline double
   priorInvVar(const int& i) const
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::priorInvVar(i).", 1);
    return _priorInvVar[i];
   }

   inline int
   lcovFixed() const{ return _lcovFixed;}

   inline int
   lcovgamma() const{ return _lcovgamma;}

   inline int
   diagIFixed(const int& i) const
   {
    if (i >= _lcovFixed || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::diagIFixed(i).", 1);
    return _diagIFixed[i];
   }  

   inline int
   diagIgamma(const int& i) const
   {
    if (i >= _lcovgamma || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::diagIgamma(i).", 1);
    return _diagIgamma[i];
   }  

   inline int
   indbA(const int& i) const
   {
    if (i >= _nbeta || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::indbA(i).", 1);
    return _indbA[i];
   }  

   inline int
   indFixed(const int& i) const
   {
    if (i >= _nFixed || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::indFixed(i).", 1);
    return _indFixed[i];
   }  

   inline int
   indgamma(const int& i) const
   {
    if (i >= _ngamma || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::indgamma(i).", 1);
    return _indgamma[i];
   }  

   inline int*
   indbinXP() const { return _indbinXA; }

   inline int
   indbinXA(const int& i) const
   {
    if (i >= _nRandom || i < 0) throw returnR("C++ Error: Incorrect i in BetaGamma::indbinXA(i).", 1);
    return _indbinXA[i];
   }  

   void
   GIBBSfixed(double* regresResM,  const int* nP,
              const double* XA,    const double* XXtb,
              const Gspline* gg,   double** const mu,       const int* rM);

};    // end of class BetaGamma

#endif
