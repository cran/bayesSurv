// Derived class from BetaGamma
//  * added method to update means of random effects
//  * it has to be done via derived class since otherwise we had chicken-egg problem
//    - method RandomEff::GIBBSupdate has one of its arguments an object of class BetaGamma
//    - method BetaGammaExtend::GIBBSmeanRandom has one of its arguments an object of class RandomEff
//    - so without using derived classes approach, at the moment of declaration of method BetaGamma::GIBBSmeanRandom
//      inside the class BetaGamma, the compiler needs to know what class RandomEff is
//      and vice versa at the moment of declaration of method RandomEff::GIBBSupdate inside the class RandomEff,
//      the compiler needs to know what class BetaGamma is
//
#ifndef _CLASS_BETA_GAMMA_EXTEND_H_
#define _CLASS_BETA_GAMMA_EXTEND_H_

#include "classBetaGamma.h"
#include "classRandomEff.h"

class BetaGammaExtend : public BetaGamma
{
  private:

  public:
    BetaGammaExtend();
    BetaGammaExtend(const int* parmI,  const double* parmD);
    BetaGammaExtend(const BetaGammaExtend& bg);
    BetaGammaExtend(const BetaGamma& bg);
    BetaGammaExtend& operator=(const BetaGammaExtend& bg);

    void
    GIBBSmeanRandom(const RandomEff* b_obj,  const CovMatrix* Dcm);
}; 

#endif
