// Class to store random effects together with the prior information on them
//   and methods to update them
// * it should replace 'bblocks' class in all future version sof bayessurvreg
//
#ifndef _CLASS_RANDOM_EFF_H_
#define _CLASS_RANDOM_EFF_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "constants.h"
#include "mvtdist.h"
#include "random.h"
#include "cholesky.h"

#include "Gspline.h"
#include "classBetaGamma.h"
#include "classCovMatrix.h"

const double _Eb0_ = 0.0;                // mean of the random intercept

class RandomEff
{
  private:
    int _nRandom;       // dimension of each random effect
    int _nCluster;      // number of clusters
    int _lbMarray;      // _nRandom*_nCluster
    int _larray;        // (_nRandom * (_nRandom + 1))/2
    int* _nwithinCl;    // [_nCluster] number of observations within each cluster

    int _type_prior;    // type of the prior distribution: 0 = Normal, 1 = G-spline (not yet implemented)

    double* _bM;        // [_lbMarray] array with random effects stored as 
                        // (b[0,0], ..., b[0, _nRandom-1], ..., b[_nCluster-1,0],..., b[_nCluster-1,_nRandom-1])
    
    int* _diagI;           // [_nRandom]
    double* _covpar;       // [_larray] proposal covariance matrix (working space)
    double* _ichicovpar;   // [_larray] working space
    int* _indUpd;          // [_nRandom] array filled with 0, 1, ..., _nRandom-1; used in the call to 'rmvtnorm2'

    double* _Digamma;        // [_nRandom] working space for update of the random effects
    double* _propMean;       // [_nRandom] working space for update of either the random effects or used by the method 'sumSquare'
    double* _propMeanTemp;   // [_nRandom] working space for update of the random effects

  public:
    RandomEff();
    RandomEff(const int* parmI,  const double* parmD);
    RandomEff(const RandomEff& re);
    RandomEff&
    operator=(const RandomEff& re);
    ~RandomEff();

    void print() const;

    void RandomEff2initArray(int* parmI,  double* parmD) const;

    inline int
    nRandom() const { return _nRandom;}

    inline int
    nCluster() const { return _nCluster;}

    inline int
    nwithinCl(const int& i) const 
    { 
      if (i < 0 || i >= _nCluster) throw returnR("Error: i out of range in RandomEff::nwithinCl(i)", 1);
      return _nwithinCl[i];
    }

    inline int
    type_prior() const { return _type_prior;}
   
    inline int
    diagI(const int& i) const
    {
      if (i < 0 || i >= _nRandom) throw returnR("Error: i out of range in RandomEff::diagI(i)", 1);
      return _diagI[i];
    }

    inline int
    lbMarray() const { return _lbMarray; }

    inline int
    larray() const { return _larray; }

    inline double*
    bMP() const { return _bM; }

    inline int*
    nwithinClP() const { return _nwithinCl;}

    void
    GIBBSupdate(double* regresResM,    const int* nP,
                const double* XA,      const double* ZZtb,
                const Gspline* gg,     double** const mu,    const int* rM,
                const BetaGamma* bbg,  const CovMatrix* Dcm);

    void
    sumSquare(double* sumSq,  const BetaGamma* bbg) const;

    void
    predictNormalRE(const BetaGamma* bg,  const CovMatrix* DD);

    void
    Gspl_intcpt_update(double* regresResM,    const int* nP,
                       const Gspline* gg_b,   double** const mu_b,       const int* rM_b,
                       const Gspline* gg,     double** const mu,         const int* rM);

    void
    adjust_intcpt(const double* adj);

    void
    predictGspl_intcpt(const int* k_effect,      double* cum_w,  const double* prop_mu, 
                       const double* sig_scale,  int* rM_b);

    friend void
    Gspl_rho_intcpt_update(RandomEff *d,            RandomEff *b,              double *rho_zb,
                           double *regResOnset,     double *regResTime,        int *rho_accept,
                           const int *nP,           const int *rho_algor,      const double *rho_scaleL,
                           const Gspline *gg_d,     double** const mu_d,       const int *rM_d,
                           const Gspline *gg_b,     double** const mu_b,       const int *rM_b,
                           const Gspline *gg_zeta,  double** const mu_zeta,    const int *rM_zeta,
                           const Gspline *gg_eps,   double** const mu_eps,     const int *rM_eps);
};


#endif


