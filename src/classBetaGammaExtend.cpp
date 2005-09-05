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
// 29/01/2005: 'BetaGammaExtend::GIBBSmeanRandom'
//
#include "classBetaGammaExtend.h"

/***** Nonparametric constructor *****/
BetaGammaExtend::BetaGammaExtend() 
  : BetaGamma()
{}


/***** Parametric constructor  *****/
BetaGammaExtend::BetaGammaExtend(const int* parmI,  const double* parmD) 
  : BetaGamma(parmI, parmD)
{}


/***** Copy constructor 1 *****/
BetaGammaExtend::BetaGammaExtend(const BetaGammaExtend& bg) 
  : BetaGamma(bg)
{}


/***** Copy constructor 2 *****/
BetaGammaExtend::BetaGammaExtend(const BetaGamma& bg) 
  : BetaGamma(bg)
{}


/***** Assignment operator *****/
BetaGammaExtend& 
BetaGammaExtend::operator=(const BetaGammaExtend& bg)
{
  int i;

  if (!_nbeta && _randomIntcpt){
    free(_indbinXA);
  }
  if (_nbeta){
    free(_indbA);
    free(_beta);     free(_priorMean);   free(_priorSD);   free(_priorInvVar);
    if (_nFixed){
      free(_indFixed);
      free(_meanFixed);   free(_meanFixedTemp);   free(_covFixed);   free(_ichicovFixed);   free(_diagIFixed);    
    }
    if (_ngamma){
      free(_indgamma);
      free(_meangamma);   free(_meangammaTemp);   free(_covgamma);   free(_ichicovgamma);   free(_diagIgamma);
      free(_sumbM);       free(_indRandomUpdate);
      if (_nRandom > _ngamma){
        free(_sumgammab);  free(_indRandomKeep);
      }
    }
    if (_nRandom){
      free(_indbinXA);
    }
  }

  if (bg.nbeta() == 0){
    _nbeta =_nFixed = _ngamma = 0;  
    _indbA = _indFixed = NULL;
    _beta = _priorMean = _priorSD = _priorInvVar = NULL; 
    _lcovFixed = 0;  _meanFixed = _meanFixedTemp = _covFixed = _ichicovFixed = NULL;  _diagIFixed = NULL;
    _lcovgamma = 0;  _meangamma = _meangammaTemp = _covgamma = _ichicovgamma = NULL;  _diagIgamma = NULL;
    _sumbM = _sumgammab = NULL;
    _indRandomUpdate = _indRandomKeep = NULL;

    _randomIntcpt = bg.randomIntcpt();
    _nRandom      = bg.nRandom();
    if (_randomIntcpt){
      _indbinXA = (int*) malloc(sizeof(int));
      if (!_indbinXA) throw returnR("Not enough memory available in BetaGamma assignment operator (_indbinXA)", 1);
      *_indbinXA = -1;
    }
    else{
      _indbinXA = NULL;
    }
  }
  else{
    _nbeta        = bg.nbeta();
    _nFixed       = bg.nFixed();
    _ngamma       = bg.ngamma();
    _randomIntcpt = bg.randomIntcpt();
    _nRandom      = bg.nRandom();

    _indbA = (int*) calloc(_nbeta, sizeof(int));
    if (!_indbA) throw returnR("Not enough memory available in BetaGamma assignment operator (_indbA)", 1);
    for (i = 0; i < _nbeta; i++){
      _indbA[i] = bg._indbA[i];
    }    
    if (_nFixed > 0){
      _indFixed = (int*) calloc(_nFixed, sizeof(int));    
      if (!_indFixed) throw returnR("Not enough memory available in BetaGamma assignment operator (_indFixed)", 1);
      for (i = 0; i < _nFixed; i++) _indFixed[i] = bg._indFixed[i];
    }
    if (_ngamma){
      _indgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_indgamma) throw returnR("Not enough memory available in BetaGamma assignment operator (_indgamma)", 1);
      for (i = 0; i < _ngamma; i++) _indgamma[i] = bg._indgamma[i];
    }
    if (_nRandom > 0){
      _indbinXA = (int*) calloc(_nRandom, sizeof(int));    
      if (!_indbinXA) throw returnR("Not enough memory available in BetaGamma copy constructor (_indbinXA)", 1);
      for (i = 0; i < _nRandom; i++) _indbinXA[i] = bg._indbinXA[i];
    }

    _beta        = (double*) calloc(_nbeta, sizeof(double));
    _priorMean   = (double*) calloc(_nbeta, sizeof(double));
    _priorSD     = (double*) calloc(_nbeta, sizeof(double));
    _priorInvVar = (double*) calloc(_nbeta, sizeof(double));
    if (!_beta) throw returnR("Not enough memory available in BetaGamma assignment operator (_beta)", 1);
    if (!_priorMean || !_priorSD || !_priorInvVar) throw returnR("Not enough memory available in BetaGamma assign. oper. (_prior*)", 1);
    for (i = 0; i < _nbeta; i++){
      _beta[i]        = bg.beta(i);
      _priorMean[i]   = bg.priorMean(i); 
      _priorSD[i]     = bg.priorSD(i);
      _priorInvVar[i] = bg.priorInvVar(i);
    }

    _lcovFixed = bg.lcovFixed();
    if (_nFixed){
      _meanFixed     = (double*) calloc(_nFixed, sizeof(double));
      _meanFixedTemp = (double*) calloc(_nFixed, sizeof(double));
      if (!_meanFixed || !_meanFixedTemp) throw returnR("Not enough memory available in BetaGamma assignment operator (_meanFixed*)", 1);
      for (i = 0; i < _nFixed; i++){
        _meanFixed[i] = bg._meanFixed[i];
        _meanFixedTemp[i] = bg._meanFixedTemp[i];
      }

      _covFixed     = (double*) calloc(_lcovFixed, sizeof(double));
      _ichicovFixed = (double*) calloc(_lcovFixed, sizeof(double));
      if (!_covFixed || !_ichicovFixed) throw returnR("Not enough memory available in BetaGamma assignment operator (_*covFixed)", 1);
      for (i = 0; i < _lcovFixed; i++){
        _covFixed[i] = bg._covFixed[i];
        _ichicovFixed[i] = bg._ichicovFixed[i];
      }

      _diagIFixed = (int*) calloc(_nFixed, sizeof(int));
      if (!_diagIFixed) throw returnR("Not enough memory available in BetaGamma assignment operator (_diagIFixed)", 1);
      for (i = 0; i < _nFixed; i++) _diagIFixed[i] = bg._diagIFixed[i];
    }
    else{
      _meanFixed = _meanFixedTemp = _covFixed = _ichicovFixed = NULL;
      _diagIFixed = NULL;
    }        

    _lcovgamma = bg.lcovgamma();
    if (_ngamma){
      _meangamma     = (double*) calloc(_ngamma, sizeof(double));
      _meangammaTemp = (double*) calloc(_ngamma, sizeof(double));
      if (!_meangamma || !_meangammaTemp) throw returnR("Not enough memory available in BetaGamma assignment operator (_meangamma*)", 1);
      for (i = 0; i < _ngamma; i++){
        _meangamma[i] = bg._meangamma[i];
        _meangammaTemp[i] = bg._meangammaTemp[i];
      }

      _covgamma     = (double*) calloc(_lcovgamma, sizeof(double));
      _ichicovgamma = (double*) calloc(_lcovgamma, sizeof(double));
      if (!_covgamma || !_ichicovgamma) throw returnR("Not enough memory available in BetaGamma assignment operator (_*covgamma)", 1);
      for (i = 0; i < _lcovgamma; i++){
        _covgamma[i] = bg._covgamma[i];
        _ichicovgamma[i] = bg._ichicovgamma[i];
      }

      _diagIgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_diagIgamma) throw returnR("Not enough memory available in BetaGamma assignment operator (_diagIgamma)", 1);
      for (i = 0; i < _ngamma; i++) _diagIgamma[i] = bg._diagIgamma[i];

      _sumbM           = (double*) calloc(_ngamma, sizeof(double));
      _indRandomUpdate = (int*) calloc(_ngamma, sizeof(int));
      if (!_sumbM || !_indRandomUpdate) throw returnR("Not enough memory available in BetaGamma copy const. (_sumbM/_indRandomUpdate)", 1);
      for (i = 0; i < _ngamma; i++){
        _sumbM[i] = bg._sumbM[i];
        _indRandomUpdate[i] = bg._indRandomUpdate[i];
      }
      if (_nRandom > _ngamma){
        _sumgammab = (double*) calloc(_nRandom - _ngamma, sizeof(double));
        _indRandomKeep = (int*) calloc(_nRandom - _ngamma, sizeof(int));
        if (!_sumgammab || !_indRandomKeep) throw returnR("Not enough memory avail. in BetaGamma copy con. (_sumgammab/_indRandomKeep)", 1);
        _sumgammab[0]     = bg._sumgammab[0];        // _nRandom - _ngamma must be equal to 1
        _indRandomKeep[0] = bg._indRandomKeep[0];    // this must be an index of the random intercept
      }
      else{
        _sumgammab = NULL;
        _indRandomKeep = NULL;
      }
    }
    else{
      _meangamma = _meangammaTemp = _covgamma = _ichicovgamma = NULL;
      _diagIgamma = NULL;
      _sumbM = _sumgammab = NULL;
      _indRandomUpdate = _indRandomKeep = NULL;
    }
  }
  
  return *this;  
}    /** end of the assignment operator **/


//
// ***** GIBBSmeanRandom *****
//
// Update all means of random effects using a Gibbs move
//
void
BetaGammaExtend::GIBBSmeanRandom(const RandomEff* b_obj,  const CovMatrix* Dcm)
{
  if (!_ngamma) return;

  static int i, j, ii, jj, cl, rank;

    /** Inverse variance of full conditional distribution (Psi^{-1} + N*D^{-1}) (store it in _covgamma) AND           **/
    /** mean of the full conditional distribution, part 1 (Psi^{-1}*nu)                                               **/
  for (j = 0; j < _ngamma; j++){

      /* Diagonal */
    jj = _indbA[_indgamma[j]];
    if (jj < 0 || jj >= b_obj->nRandom()) throw returnR("BetaGammaExtend::GIBBSmeanRandom: Programming error, contact the author", 99);
    _covgamma[_diagIgamma[j]] = _priorInvVar[_indgamma[j]] + b_obj->nCluster() * (Dcm->icovm(Dcm->diagI(jj)));

      /* Off-diagonal in the jth column*/
    for (i = j + 1; i < _ngamma; i++){
      ii = _indbA[_indgamma[i]];
      if (ii > jj) _covgamma[_diagIgamma[j] + i - j] = b_obj->nCluster() * (Dcm->icovm(Dcm->diagI(jj) + ii - jj));
      else         _covgamma[_diagIgamma[j] + i - j] = b_obj->nCluster() * (Dcm->icovm(Dcm->diagI(ii) + jj - ii));
    }

      /* Part 1 of the mean */
    _meangammaTemp[j] = _priorInvVar[_indgamma[j]] * _priorMean[_indgamma[j]];
  }


    /** Cholesky decomposition of the inverse variance of full conditional distrib. **/
  cholesky(_covgamma, &rank, &_ngamma, _diagIgamma, &_toler_chol_BetaGamma);
 
    /** Variance of the full conditional distribution                               **/
    /**  and the inverse of the Cholesky decomposition of the inverse variance      **/
  chinv2(_covgamma, _ichicovgamma, &_ngamma, _diagIgamma);

    /** Mean of the full conditional distribution, part 2 (+ V_M*\sum b_M - W*\sum(gamma_{-M} - b_{-M}))   **/
  const double* bb;

    /*  a) \sum b_M (store it in _sumbM)                                                                     */
  for (j = 0; j < _ngamma; j++) _sumbM[j] = 0.0;
  bb = b_obj->bMP();
  for (cl = 0; cl < b_obj->nCluster(); cl++){
    for (j = 0; j < _ngamma; j++) _sumbM[j] += bb[_indbA[_indgamma[j]]];
    bb += b_obj->nRandom();
  }
  
    /*  b) += V_M * \sum b_M (store it first in _meangamma)                                                 */
  Mxa2(_meangamma, _sumbM, Dcm->icovmP(), _indRandomUpdate, &_ngamma, &_nRandom, Dcm->diagIP());    
  for (j = 0; j < _ngamma; j++) _meangammaTemp[j] += _meangamma[j];

    /*  c) \sum (gamma_{-M} - b_{-M}) (store it in _sumgammab)                                              */
    /*  d) -= W * \sum(gamma_{-M} - b_{-M}) (store it first in _meangamma)                                  */
  jj = _nRandom - _ngamma;
  if (jj > 0){
    if (jj != 1) throw returnR("Programming error in BetaGammaExtend::GIBBSmeanRandom, contact the author", 1);
    _sumgammab[0] = 0.0;
    bb = b_obj->bMP();
    for (cl = 0; cl < b_obj->nCluster(); cl++){
      _sumgammab[0] += (_Eb0_ - bb[0]);
      bb += b_obj->nRandom();    
    }
    Wxa(_meangamma, _sumgammab, Dcm->icovmP(), _indRandomUpdate, _indRandomKeep, &jj, &_nRandom, &_ngamma, Dcm->diagIP());
    for (j = 0; j < _ngamma; j++) _meangammaTemp[j] -= _meangamma[j];    
  }

    /** Mean of full conditional distribution, part 3 (* var(gamma(M)|...))    **/
  Mxa(_meangamma, _meangammaTemp, _covgamma, &ZERO_INT, &_ngamma, &_ngamma, _diagIgamma);

    /** Sample  **/
  rmvtnorm2(_beta, _meangamma, _ichicovgamma, &ZERO_INT, _indgamma, &_nbeta, &_ngamma, &_ngamma, &ONE_INT, _diagIgamma, &ZERO_INT);

  return;
}    /*** end of the function BetaGammaExtend::GIBBSmeanRandom  ***/

