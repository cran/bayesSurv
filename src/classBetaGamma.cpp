// Class to store regression parameters and means of random effects
//  *** to be used with G-spline error distribution
//  *** this should replace MHblocks in all future functions
//  *** it is intended to have somewhat simpler structure than MHblocks
//  *** only Gibbs move updating all parameters at once will be implemented
// 
// 20/12/2004: basics + GIBBSfixed
// 11/01/2005: basics finished
// 
//
#include "classBetaGamma.h"

/***** Nonparametric constructor *****/
BetaGamma::BetaGamma()
  : _nbeta(0), _nFixed(0), _ngamma(0), _randomIntcpt(0), _nRandom(0), 
    _beta(NULL),
    _indbA(NULL), _indFixed(NULL), _indgamma(NULL), _indbinXA(NULL),
    _priorMean(NULL), _priorSD(NULL), _priorInvVar(NULL),
    _lcovFixed(0), _meanFixed(NULL), _meanFixedTemp(NULL), _covFixed(NULL), _ichicovFixed(NULL), _diagIFixed(NULL),
    _lcovgamma(0), _meangamma(NULL), _meangammaTemp(NULL), _covgamma(NULL), _ichicovgamma(NULL), _diagIgamma(NULL),
    _sumbM(NULL), _sumgammab(NULL), 
    _indRandomUpdate(NULL), _indRandomKeep(NULL)
{  
}

/***** Parametric constructor  *****/
BetaGamma::BetaGamma(const int* parmI,  const double* parmD)
{
  int i;

  /** Set up integer pointers **/
  int pnbeta        = 0;
  int pnFixed       = pnbeta + 1;
  int pngamma       = pnFixed + 1;
  int prandomIntcpt = pngamma + 1;
  int pindbA        = prandomIntcpt + 1;
  // int pnext        = pindbA + parmI[pnbeta];

  /** Set up double pointers **/
  int pbeta         = 0;
  int ppriorMean    = pbeta + parmI[pnbeta]; 
  int ppriorVar     = ppriorMean + parmI[pnbeta]; 
  // int pnext        = ppriorVar + parmI[pnbeta];

  if (parmI[pnbeta] <= 0){
    _nbeta =_nFixed = _ngamma = 0;  
    _indbA = _indFixed = _indgamma = NULL;
    _beta = _priorMean = _priorSD = _priorInvVar = NULL; 
    _lcovFixed = 0;  _meanFixed = _meanFixedTemp = _covFixed = _ichicovFixed = NULL;  _diagIFixed = NULL;
    _lcovgamma = 0;  _meangamma = _meangammaTemp = _covgamma = _ichicovgamma = NULL;  _diagIgamma = NULL;
    _sumbM = _sumgammab = NULL;
    _indRandomUpdate = _indRandomKeep = NULL;

    _randomIntcpt = parmI[prandomIntcpt];
    _nRandom      = _randomIntcpt;
    if (_randomIntcpt){
      _indbinXA = (int*) malloc(sizeof(int));
      if (!_indbinXA) throw returnR("Not enough memory available in BetaGamma constructor (_indbinXA)", 1);
      *_indbinXA = -1;
    }
    else{
      _indbinXA = NULL;
    }
  }
  else{

    /** Dimension parameters **/
    _nbeta = parmI[pnbeta];
    if (parmI[pnFixed] < 0 || parmI[pnFixed] > _nbeta) throw returnR("C++ Error 1 in BetaGamma constructor ", 1);
    _nFixed = parmI[pnFixed];
    if (parmI[pngamma] < 0 || parmI[pngamma] > _nbeta) throw returnR("C++ Error 2 in BetaGamma constructor ", 1);
    _ngamma = parmI[pngamma];
    if (parmI[prandomIntcpt] < 0 || parmI[prandomIntcpt] > 1) throw returnR("C++ Error 3 in BetaGamma constructor ", 1);
    _randomIntcpt = parmI[prandomIntcpt];
    if (_nFixed + _ngamma != _nbeta) throw returnR("C++ Error 4 in BetaGamma constructor ", 1);
    _nRandom = _ngamma + _randomIntcpt;

    /** Indexing arrays         **/
    _indbA = (int*) calloc(_nbeta, sizeof(int));
    if (!_indbA) throw returnR("Not enough memory available in BetaGamma constructor (_indbA)", 1);
    if (_nFixed > 0){
      _indFixed = (int*) calloc(_nFixed, sizeof(int));    
      if (!_indFixed) throw returnR("Not enough memory available in BetaGamma constructor (_indFixed)", 1);
    }
    if (_ngamma){
      _indgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_indgamma) throw returnR("Not enough memory available in BetaGamma constructor (_indgamma)", 1);
    }
    if (_nRandom > 0){
      _indbinXA = (int*) calloc(_nRandom, sizeof(int));    
      if (!_indbinXA) throw returnR("Not enough memory available in BetaGamma constructor (_indbinXA)", 1);
    }

    int ifixed = 0;
    int irandom = _randomIntcpt;
    if (_randomIntcpt) _indbinXA[0] = -1;
    for (i = 0; i < _nbeta; i++){
      if (parmI[pindbA + i] < -1) throw returnR("C++ Error 5 in BetaGamma constructor ", 1);
      if (parmI[pindbA + i] == -1){
        if (ifixed >= _nFixed) throw returnR("C++ Error 6 in BetaGamma constructor ", 1);
        _indbA[i] = -1;
        _indFixed[ifixed] = i;
        ifixed++;
      }
      else{
        if (parmI[pindbA + i] == 0 && _randomIntcpt) throw returnR("C++ Error 7 in BetaGamma constructor ", 1);
        if (parmI[pindbA + i] >= _nRandom) throw returnR("C++ Error 8 in BetaGamma constructor ", 1);
        if (irandom >= _nRandom) throw returnR("C++ Error 9 in BetaGamma constructor ", 1);
        _indbA[i] = parmI[pindbA + i];
        _indbinXA[irandom] = i;
        _indgamma[irandom - _randomIntcpt] = i;
        irandom++;
      }
    }
    if (ifixed != _nFixed || irandom != _nRandom) throw returnR("C++ Error 10 in BetaGamma constructor ", 1);

    /** Priors and initials **/
    _beta        = (double*) calloc(_nbeta, sizeof(double));
    _priorMean   = (double*) calloc(_nbeta, sizeof(double));
    _priorSD     = (double*) calloc(_nbeta, sizeof(double));
    _priorInvVar = (double*) calloc(_nbeta, sizeof(double));
    if (!_beta) throw returnR("Not enough memory available in BetaGamma constructor (_beta)", 1);
    if (!_priorMean || !_priorSD || !_priorInvVar) throw returnR("Not enough memory available in BetaGamma constructor (_prior*)", 1);
    for (i = 0; i < _nbeta; i++){
      _beta[i]      = parmD[pbeta + i];
      _priorMean[i] = parmD[ppriorMean + i]; 
      if (parmD[ppriorVar + i] <= 0) throw returnR("C++ Error 11 in BetaGamma constructor ", 1);
      _priorSD[i]     = sqrt(parmD[ppriorVar + i]);
      _priorInvVar[i] = 1 / parmD[ppriorVar + i];
    }

    /** Stuff for update of fixed effects **/
    _lcovFixed = (_nFixed * (_nFixed + 1))/2;
    if (_nFixed){
      _meanFixed     = (double*) calloc(_nFixed, sizeof(double));
      _meanFixedTemp = (double*) calloc(_nFixed, sizeof(double));
      if (!_meanFixed || !_meanFixedTemp) throw returnR("Not enough memory available in BetaGamma constructor (_meanFixed*)", 1);
      for (i = 0; i < _nFixed; i++) _meanFixed[i] = _meanFixedTemp[i] = 0.0;

      _covFixed     = (double*) calloc(_lcovFixed, sizeof(double));
      _ichicovFixed = (double*) calloc(_lcovFixed, sizeof(double));
      if (!_covFixed || !_ichicovFixed) throw returnR("Not enough memory available in BetaGamma constructor (_*covFixed)", 1);
      for (i = 0; i < _lcovFixed; i++) _covFixed[i] = _ichicovFixed[i] = 0.0;

      _diagIFixed = (int*) calloc(_nFixed, sizeof(int));
      if (!_diagIFixed) throw returnR("Not enough memory available in BetaGamma constructor (_diagIFixed)", 1);
      _diagIFixed[0] = 0;
      for (i = 1; i < _nFixed; i++) _diagIFixed[i] = _diagIFixed[i-1] + _nFixed - i + 1;
    }
    else{
      _meanFixed = _meanFixedTemp = _covFixed = _ichicovFixed = NULL;
      _diagIFixed = NULL;
    }

    /** Stuff for update of means of random effects **/
    _lcovgamma = (_ngamma * (_ngamma + 1))/2;
    if (_ngamma){
      _meangamma     = (double*) calloc(_ngamma, sizeof(double));
      _meangammaTemp = (double*) calloc(_ngamma, sizeof(double));
      if (!_meangamma || !_meangammaTemp) throw returnR("Not enough memory available in BetaGamma constructor (_meangamma*)", 1);
      for (i = 0; i < _ngamma; i++) _meangamma[i] = _meangammaTemp[i] = 0.0;

      _covgamma     = (double*) calloc(_lcovgamma, sizeof(double));
      _ichicovgamma = (double*) calloc(_lcovgamma, sizeof(double));
      if (!_covgamma || !_ichicovgamma) throw returnR("Not enough memory available in BetaGamma constructor (_*covgamma)", 1);
      for (i = 0; i < _lcovgamma; i++) _covgamma[i] = _ichicovgamma[i] = 0.0;

      _diagIgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_diagIgamma) throw returnR("Not enough memory available in BetaGamma constructor (_diagIgamma)", 1);
      _diagIgamma[0] = 0;
      for (i = 1; i < _ngamma; i++) _diagIgamma[i] = _diagIgamma[i-1] + _ngamma - i + 1;

      _sumbM           = (double*) calloc(_ngamma, sizeof(double));
      _indRandomUpdate = (int*) calloc(_ngamma, sizeof(int));
      if (!_sumbM || !_indRandomUpdate) throw returnR("Not enough memory available in BetaGamma constructor (_sumbM/_indRandomUpdate)", 1);
      for (i = 0; i < _ngamma; i++){
        _sumbM[i] = 0.0;
        _indRandomUpdate[i] = _indbA[_indgamma[i]];
      }
      if (_nRandom > _ngamma){
        _sumgammab = (double*) calloc(_nRandom - _ngamma, sizeof(double));
        _indRandomKeep = (int*) calloc(_nRandom - _ngamma, sizeof(int));
        if (!_sumgammab || !_indRandomKeep) throw returnR("Not enough memory available in BetaGamma const. (_sumgammab/_indRandomKeep)", 1);
        _sumgammab[0]     = 0.0;      // _nRandom - _ngamma must be equal to 1
        _indRandomKeep[0] = 0;        // this must be an index of the random intercept
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
}    /** end of the parametric constructor **/


/***** Copy constructor *****/
BetaGamma::BetaGamma(const BetaGamma& bg)
{
  int i;

  if (bg.nbeta() == 0){
    _nbeta =_nFixed = _ngamma = 0;  
    _indbA = _indFixed = _indgamma = NULL;
    _beta = _priorMean = _priorSD = _priorInvVar = NULL; 
    _lcovFixed = 0;  _meanFixed = _meanFixedTemp = _covFixed = _ichicovFixed = NULL;  _diagIFixed = NULL;
    _lcovgamma = 0;  _meangamma = _meangammaTemp = _covgamma = _ichicovgamma = NULL;  _diagIgamma = NULL;
    _sumbM = _sumgammab = NULL;
    _indRandomUpdate = _indRandomKeep = NULL;

    _randomIntcpt = bg.randomIntcpt();
    _nRandom      = bg.nRandom();
    if (_randomIntcpt){
      _indbinXA = (int*) malloc(sizeof(int));
      if (!_indbinXA) throw returnR("Not enough memory available in BetaGamma copy constructor (_indbinXA)", 1);
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
    if (!_indbA) throw returnR("Not enough memory available in BetaGamma copy constructor (_indbA)", 1);
    for (i = 0; i < _nbeta; i++){
      _indbA[i] = bg._indbA[i];
    }    
    if (_nFixed > 0){
      _indFixed = (int*) calloc(_nFixed, sizeof(int));    
      if (!_indFixed) throw returnR("Not enough memory available in BetaGamma copy constructor (_indFixed)", 1);
      for (i = 0; i < _nFixed; i++) _indFixed[i] = bg._indFixed[i];
    }
    if (_ngamma){
      _indgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_indgamma) throw returnR("Not enough memory available in BetaGamma copy constructor (_indgamma)", 1);
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
    if (!_beta) throw returnR("Not enough memory available in BetaGamma copy constructor (_beta)", 1);
    if (!_priorMean || !_priorSD || !_priorInvVar) throw returnR("Not enough memory available in BetaGamma copy constructor (_prior*)", 1);
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
      if (!_meanFixed || !_meanFixedTemp) throw returnR("Not enough memory available in BetaGamma copy constructor (_meanFixed*)", 1);
      for (i = 0; i < _nFixed; i++){
        _meanFixed[i] = bg._meanFixed[i];
        _meanFixedTemp[i] = bg._meanFixedTemp[i];
      }

      _covFixed     = (double*) calloc(_lcovFixed, sizeof(double));
      _ichicovFixed = (double*) calloc(_lcovFixed, sizeof(double));
      if (!_covFixed || !_ichicovFixed) throw returnR("Not enough memory available in BetaGamma copy constructor (_*covFixed)", 1);
      for (i = 0; i < _lcovFixed; i++){
        _covFixed[i] = bg._covFixed[i];
        _ichicovFixed[i] = bg._ichicovFixed[i];
      }

      _diagIFixed = (int*) calloc(_nFixed, sizeof(int));
      if (!_diagIFixed) throw returnR("Not enough memory available in BetaGamma copy constructor (_diagIFixed)", 1);
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
      if (!_meangamma || !_meangammaTemp) throw returnR("Not enough memory available in BetaGamma copy constructor (_meangamma*)", 1);
      for (i = 0; i < _ngamma; i++){
        _meangamma[i] = bg._meangamma[i];
        _meangammaTemp[i] = bg._meangammaTemp[i];
      }

      _covgamma     = (double*) calloc(_lcovgamma, sizeof(double));
      _ichicovgamma = (double*) calloc(_lcovgamma, sizeof(double));
      if (!_covgamma || !_ichicovgamma) throw returnR("Not enough memory available in BetaGamma copy constructor (_*covgamma)", 1);
      for (i = 0; i < _lcovgamma; i++){
        _covgamma[i] = bg._covgamma[i];
        _ichicovgamma[i] = bg._ichicovgamma[i];
      }

      _diagIgamma = (int*) calloc(_ngamma, sizeof(int));
      if (!_diagIgamma) throw returnR("Not enough memory available in BetaGamma copy constructor (_diagIgamma)", 1);
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
}    /** end of the copy constructor **/


/***** Assignment operator *****/
BetaGamma&
BetaGamma::operator=(const BetaGamma& bg)
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


/***** Destructor *****/
BetaGamma::~BetaGamma()
{
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
}    /** end of the destructor **/


/***** BetaGamma2initArray *****/
void
BetaGamma::BetaGamma2initArray(int* parmI,  double* parmD) const
{
  int i;

  /** Set up integer pointers **/
  int pnbeta        = 0;
  int pnFixed       = pnbeta + 1;
  int pngamma       = pnFixed + 1;
  int prandomIntcpt = pngamma + 1;
  int pindbA        = prandomIntcpt + 1;
  // int pnext        = pindbA + parmI[pnbeta];

  /** Set up double pointers **/
  int pbeta         = 0;
  int ppriorMean    = pbeta + parmI[pnbeta]; 
  int ppriorVar     = ppriorMean + parmI[pnbeta]; 
  // int pnext        = ppriorVar + parmI[pnbeta];

  parmI[pnbeta]        = _nbeta;
  parmI[pnFixed]       = _nFixed;
  parmI[pngamma]       = _ngamma;
  parmI[prandomIntcpt] = _randomIntcpt;
  
  for (i = 0; i < _nbeta; i++){
    parmI[pindbA + i]     = _indbA[i];
    parmD[pbeta + i]      = _beta[i];
    parmD[ppriorMean + i] = _priorMean[i];
    parmD[ppriorVar + i]  = _priorSD[i] * _priorSD[i];
  }

  return;
}


//
// ***** GIBBSfixed *****
//
// Update all corresponding fixed effects using a Gibbs move
// Update also regression residuals
//
// regresResM ... [nP*gg->dim()] regression residuals y - beta'x - b'z 
//                  !!! possible intercept is not subtracted from y !!!
// nP ........... a) for univariate AFT model with random effects: number of all observations
//                b) for bivariate AFT model: number of pairs (bivariate observations)
// XA ........... [nobs*_nbeta] full design matrix for all observations
//                CHANGE COMPARED TO bayessurvreg1:
//                * stored in row-major order, i.e. XA[0,     ...,_nbeta-1]   = covariates for the first observation
//                                                  XA[_nbeta,...,2*_nbeta-1] = covariates for the second observation etc.
// XXtb ......... [nobs*_lcovFixed] x*x' (only components corresponding to fixed effects, only triangles) for each observation
//                CHANGE COMPARED TO bayessurvreg1:
//                * not stored in double-indexed array
// gg ........... G-spline giving the distribution of observations
// mu ........... [gg->dim(), gg->length(j)] already computed knots of the G-spline
// rM ........... [nP] allocation labels taking values from {0,..., g_y->total_length() - 1}
//
void
BetaGamma::GIBBSfixed(double* regresResM,  const int* nP,
                      const double* XA,    const double* XXtb,   
                      const Gspline* gg,   double** const mu,       const int* rM)
{
  if (!_nFixed) return;

  static int i, m, j, obs, helpi, rank;
  static double helpd;
  static double invsigscale2[_max_dim];
  

  /** Compute invsigscale2 **/
  for (j = 0; j < gg->dim(); j++) invsigscale2[j] = gg->invscale2(j) * gg->invsigma2(j);

  /** Mean of full conditional distribution, part 1 (Psi^{-1}*nu) AND                            **/
  /** Inverse variance of full conditional distribution, part 1 (Psi^{-1}) (store it in covpar)  **/
  for (m = 0; m < _nFixed; m++){
    _meanFixedTemp[m] = _priorInvVar[_indFixed[m]] * _priorMean[_indFixed[m]];
    _covFixed[_diagIFixed[m]] = _priorInvVar[_indFixed[m]];
    for (i = m + 1; i < _nFixed; i++){
      _covFixed[_diagIFixed[m] + i - m] = 0.0;
    }
  }

  /** Loop over observations **/
  double* regRes = regresResM;
  const double* Xm = XA;
  const double* XXm = XXtb;
  const int* rp = rM;
  for (obs = 0; obs < *nP; obs++){   

      /* Compute y - alpha - z'b - scale*mu, store it in regresResM */
    for (j = 0; j < gg->dim(); j++){
      helpi = j*_nbeta;
      for (m = 0; m < _nFixed; m++){
        *(regRes+j) += (*(Xm + _indFixed[m] + helpi)) * _beta[_indFixed[m]];
      }
    }
    switch (gg->dim()){
    case 1:
      *regRes -= (gg->intcpt(0) + gg->scale(0)*mu[0][*rp]);
      break;
    case 2:
      *regRes     -= (gg->intcpt(0) + gg->scale(0)*mu[0][*rp % gg->length(0)]);
      *(regRes+1) -= (gg->intcpt(1) + gg->scale(1)*mu[1][*rp / gg->length(0)]);
      break;
    default:
      throw returnR("C++ Error: BetaGamma::GIBBSfixed not implemented for G-spline dimension > 2", 1);
    }

      /* Inverse variance of full conditional distribution, part 2 (+ 1/sigma2 * x(M) * x(M)') */
    for (j = 0 ; j < gg->dim(); j++){
      helpi = j*_lcovFixed;
      for (m = 0; m < _lcovFixed; m++){
	_covFixed[m] += invsigscale2[j] * (*(XXm + m + helpi));
      }
    }

      /* Mean of full conditional distribution, part 2 (+ 1/sigma2 * x(M) * v) */
    for (j = 0; j < gg->dim(); j++){
      helpi = j*_nbeta;
      helpd = invsigscale2[j] * (*(regRes+j));
      for (m = 0; m < _nFixed; m++){
        _meanFixedTemp[m] += helpd * (*(Xm + _indFixed[m] + helpi));
      }
    }
  
    regRes += gg->dim();
    Xm     += gg->dim()*_nbeta;
    XXm    += gg->dim()*_lcovFixed;
    rp++;
  }  /** end of the loop over observations **/

  /** Cholesky decomposition of the inverse variance of full conditional distrib.  **/
  cholesky(_covFixed, &rank, &_nFixed, _diagIFixed, &_toler_chol_BetaGamma);
 
  /** Variance of the full conditional distribution                                **/
  /**  and the inverse of the Cholesky decomposition of the inverse variance       **/
  chinv2(_covFixed, _ichicovFixed, &_nFixed, _diagIFixed);

  /** Mean of full conditional distribution, part 3 (* var(beta(M)|...))           **/
  Mxa(_meanFixed, _meanFixedTemp, _covFixed, &ZERO_INT, &_nFixed, &_nFixed, _diagIFixed);

  /** Sample                                                                       **/
  rmvtnorm2(_beta, _meanFixed, _ichicovFixed, &ZERO_INT, _indFixed, &_nbeta, &_nFixed, &_nFixed, &ONE_INT, _diagIFixed, &ZERO_INT);

  /** Update regresResM                                                            **/
  regRes = regresResM;
  Xm     = XA;
  rp     = rM;
  for (obs = 0; obs < *nP; obs++){
    for (j = 0; j < gg->dim(); j++){
      helpi = j*_nbeta;
      for (m = 0; m < _nFixed; m++){
        *(regRes+j) -= (*(Xm + _indFixed[m] + helpi)) * _beta[_indFixed[m]];
      }
    }
    switch (gg->dim()){
    case 1:
      *regRes += (gg->intcpt(0) + gg->scale(0)*mu[0][*rp]);
      break;
    case 2:
      *regRes     += (gg->intcpt(0) + gg->scale(0)*mu[0][*rp % gg->length(0)]);
      *(regRes+1) += (gg->intcpt(1) + gg->scale(1)*mu[1][*rp / gg->length(0)]);
      break;
    default:
      throw returnR("C++ Error: BetaGamma::GIBBSfixed not implemented for dimension > 2", 1);
    }
    regRes += gg->dim();
    Xm     += gg->dim()*_nbeta;
    rp++;
  }

  return;
}    /*** end of function BetaGamma::GIBBSfixed ***/



