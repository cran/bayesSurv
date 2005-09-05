// Class to store random effects together with the prior information on them
//   and methods to update them
// * it should replace 'bblocks' class in all future version sof bayessurvreg
//
// 28/01/2005: basics
//             'GIBBSupdate'
// 29/01/2005: 'sumSquare'
// 31/01/2005: 'predictNormalRE'
//             'Gspl_intcpt_update'
//             'predictGspl_intcpt'
// 06/02/2005: 'predictGspl_intcpt' rewritten
// 20/02/2005: 'adjust_intcpt' added
//
#include "classRandomEff.h"

/***** Nonparametric constructor *****/
RandomEff::RandomEff()
  : _nRandom(0), _nCluster(0), _lbMarray(0), _larray(0), _nwithinCl(NULL), _type_prior(0), _bM(NULL),
    _diagI(NULL), _covpar(NULL), _ichicovpar(NULL), _indUpd(NULL),
    _Digamma(NULL), _propMean(NULL), _propMeanTemp(NULL)
{
}


/***** Parametric constructor  *****/
RandomEff::RandomEff(const int* parmI,  const double* parmD)
{
  int i;

  /** Set up integer pointers **/
  int ptype_prior = 0;
  int pnRandom    = ptype_prior + 1;
  int pnCluster   = pnRandom + 1;
  int pnwithinCl  = pnCluster + 1;
  //  int pnext = pnwithinCl + parmI[pnCluster];

  /** Set up double pointers **/
  int pbM = 0;
  // int pnext = pbM + parmI[pnRandom]*parmI[pnCluster];

  if (parmI[pnRandom] <= 0){
    _nRandom = _nCluster = _lbMarray = _larray = _type_prior = 0;
    _nwithinCl = _diagI = _indUpd = NULL;
    _bM = _covpar = _ichicovpar = _Digamma = _propMean = _propMeanTemp = NULL;
  }
  else{
    switch (parmI[ptype_prior]){
    case Normal_:
      _type_prior = Normal_;
      break;
    case Gspline_:
      _type_prior = Gspline_;
      break;
    default:
      throw returnR("Error: Unimplemented type of prior in RandomEff constructor", 1);
    }

    /** Dimensionality parameters **/
    _nRandom   = parmI[pnRandom];
    _nCluster  = parmI[pnCluster];
    _lbMarray  = _nRandom * _nCluster;
    _nwithinCl = (int*) calloc(_nCluster, sizeof(int));
    if (!_nwithinCl) throw returnR("Not enough memory available in RandomEff constructor (_nwithinCl)", 1);
    for (i = 0; i < _nCluster; i++) _nwithinCl[i] = parmI[pnwithinCl + i];

    /** Values of random effects **/
    _bM = (double*) calloc(_lbMarray, sizeof(double));
    if (!_bM) throw returnR("Not enough memory available in RandomEff constructor (_bM)", 1);
    for (i = 0; i < _lbMarray; i++) _bM[i] = parmD[pbM + i];

    /** Stuff for update of random effects **/
    _larray     = (_nRandom * (_nRandom + 1))/2;

    _diagI      = (int*) calloc(_nRandom, sizeof(int));
    _indUpd     = (int*) calloc(_nRandom, sizeof(int));
    if (!_diagI || !_indUpd) throw returnR("Not enough memory available in RandomEff constructor (_diagI/_indUpd)", 1);
    _diagI[0]  = 0;
    _indUpd[0] = 0;
    for (i = 1; i < _nRandom; i++){
      _diagI[i]  = _diagI[i-1] + _nRandom - i + 1;
      _indUpd[i] = i;
    }

    _covpar     = (double*) calloc(_larray, sizeof(double));
    _ichicovpar = (double*) calloc(_larray, sizeof(double));
    if (!_covpar || !_ichicovpar) throw returnR("Not enough memory available in RandomEff constructor (_covpar/_ichicovpar)", 1);
    for (i = 0; i < _larray; i++) _covpar[i] = _ichicovpar[i] = 0.0;

    _Digamma      = (double*) calloc(_nRandom, sizeof(double));
    _propMean     = (double*) calloc(_nRandom, sizeof(double));
    _propMeanTemp = (double*) calloc(_nRandom, sizeof(double));
    if (!_Digamma || !_propMean || !_propMeanTemp) 
      throw returnR("Not enough memory available in RandomEff constructor (_Digamma/_propMean/_propMeanTemp)", 1);
    for (i = 0; i < _nRandom; i++) _Digamma[i] = _propMean[i] = _propMeanTemp[i] = 0.0;
  }
}  /** end of the parametric constructor **/


/***** Copy constructor *****/
RandomEff::RandomEff(const RandomEff& re)
{
  int i;

  if (re._nRandom == 0){
    _nRandom = _nCluster = _lbMarray = _larray = _type_prior = 0;
    _nwithinCl = _diagI = _indUpd = NULL;
    _bM = _covpar = _ichicovpar = _Digamma = _propMean = _propMeanTemp = NULL;
  }
  else{
    _type_prior = re._type_prior;
    _nRandom    = re._nRandom;
    _nCluster   = re._nCluster;
    _lbMarray   = re._lbMarray;
    _larray     = re._larray;

    _nwithinCl = (int*) calloc(_nCluster, sizeof(int));
    if (!_nwithinCl) throw returnR("Not enough memory available in RandomEff copy constructor (_nwithinCl)", 1);
    for (i = 0; i < _nCluster; i++) _nwithinCl[i] = re._nwithinCl[i];

    _bM = (double*) calloc(_lbMarray, sizeof(double));
    if (!_bM) throw returnR("Not enough memory available in RandomEff copy constructor (_bM)", 1);
    for (i = 0; i < _lbMarray; i++) _bM[i] = re._bM[i];

    _diagI      = (int*) calloc(_nRandom, sizeof(int));
    _indUpd     = (int*) calloc(_nRandom, sizeof(int));
    if (!_diagI || !_indUpd) throw returnR("Not enough memory available in RandomEff copy constructor (_diagI/_indUpd)", 1);
    for (i = 0; i < _nRandom; i++){
      _diagI[i]  = re._diagI[i];
      _indUpd[i] = re._indUpd[i];
    }

    _covpar     = (double*) calloc(_larray, sizeof(double));
    _ichicovpar = (double*) calloc(_larray, sizeof(double));
    if (!_covpar || !_ichicovpar) throw returnR("Not enough memory available in RandomEff copy constructor (_covpar/_ichicovpar)", 1);
    for (i = 0; i < _larray; i++){
      _covpar[i] = re._covpar[i];
      _ichicovpar[i] = re._ichicovpar[i];
    }

    _Digamma      = (double*) calloc(_nRandom, sizeof(double));
    _propMean     = (double*) calloc(_nRandom, sizeof(double));
    _propMeanTemp = (double*) calloc(_nRandom, sizeof(double));
    if (!_Digamma || !_propMean || !_propMeanTemp) 
      throw returnR("Not enough memory available in RandomEff copy constructor (_Digamma/_propMean/_propMeanTemp)", 1);
    for (i = 0; i < _nRandom; i++){
      _Digamma[i] = re._Digamma[i];
      _propMean[i] = re._propMean[i];
      _propMeanTemp[i] = re._propMeanTemp[i];
    }
  }
}  /** end of the assignment operator **/


/***** Assignment operator *****/
RandomEff&
RandomEff::operator=(const RandomEff& re)
{
  int i;

  if (_nRandom){
    free(_propMeanTemp);  free(_propMean);  free(_Digamma);
    free(_ichicovpar);    free(_covpar);
    free(_indUpd);        free(_diagI);
    free(_bM);            free(_nwithinCl);
  }

  if (re._nRandom == 0){
    _nRandom = _nCluster = _lbMarray = _larray = _type_prior = 0;
    _nwithinCl = _diagI = _indUpd = NULL;
    _bM = _covpar = _ichicovpar = _Digamma = _propMean = _propMeanTemp = NULL;
  }
  else{
    _type_prior = re._type_prior;
    _nRandom    = re._nRandom;
    _nCluster   = re._nCluster;
    _lbMarray   = re._lbMarray;
    _larray     = re._larray;

    _nwithinCl = (int*) calloc(_nCluster, sizeof(int));
    if (!_nwithinCl) throw returnR("Not enough memory available in RandomEff assignment operator (_nwithinCl)", 1);
    for (i = 0; i < _nCluster; i++) _nwithinCl[i] = re._nwithinCl[i];

    _bM = (double*) calloc(_lbMarray, sizeof(double));
    if (!_bM) throw returnR("Not enough memory available in RandomEff assignment operator (_bM)", 1);
    for (i = 0; i < _lbMarray; i++) _bM[i] = re._bM[i];

    _diagI      = (int*) calloc(_nRandom, sizeof(int));
    _indUpd     = (int*) calloc(_nRandom, sizeof(int));
    if (!_diagI || !_indUpd) throw returnR("Not enough memory available in RandomEff assignment operator (_diagI/_indUpd)", 1);
    for (i = 0; i < _nRandom; i++){
      _diagI[i]  = re._diagI[i];
      _indUpd[i] = re._indUpd[i];
    }

    _covpar     = (double*) calloc(_larray, sizeof(double));
    _ichicovpar = (double*) calloc(_larray, sizeof(double));
    if (!_covpar || !_ichicovpar) throw returnR("Not enough memory available in RandomEff assignment operator (_covpar/_ichicovpar)", 1);
    for (i = 0; i < _larray; i++){
      _covpar[i] = re._covpar[i];
      _ichicovpar[i] = re._ichicovpar[i];
    }

    _Digamma      = (double*) calloc(_nRandom, sizeof(double));
    _propMean     = (double*) calloc(_nRandom, sizeof(double));
    _propMeanTemp = (double*) calloc(_nRandom, sizeof(double));
    if (!_Digamma || !_propMean || !_propMeanTemp) 
      throw returnR("Not enough memory available in RandomEff assignment operator (_Digamma/_propMean/_propMeanTemp)", 1);
    for (i = 0; i < _nRandom; i++){
      _Digamma[i] = re._Digamma[i];
      _propMean[i] = re._propMean[i];
      _propMeanTemp[i] = re._propMeanTemp[i];
    }
  }
  return *this;
}    /** end of the assignment operator **/


/***** Destructor *****/
RandomEff::~RandomEff()
{
  if (_nRandom){
    free(_propMeanTemp);  free(_propMean);  free(_Digamma);
    free(_ichicovpar);    free(_covpar);
    free(_indUpd);        free(_diagI);
    free(_bM);            free(_nwithinCl);
  }
}  /** end of the destructor **/


/***** print   *****/
void
RandomEff::print() const
{
  Rprintf("\nRandom effects object: \n");

  if (!_nRandom){
    Rprintf("   Empty RandomEff.\n");
    return;
  }

  int j, k;
  Rprintf("   nRandom = %d,  nCluster = %d\n   nwithinCl = ", _nRandom, _nCluster);
  for (j = 0; j < _nCluster; j++) Rprintf("%d,  ", _nwithinCl[j]); Rprintf("\n");
  char zz[10];
  switch (_type_prior){
  case Normal_ : strcpy(zz, "Normal"); break;
  case Gspline_: strcpy(zz, "G-spline"); break;
  default      : strcpy(zz, "unimplemented"); break;  
  }
  Rprintf("   Distribution = %s\n", zz);
  Rprintf("   lbMarray = %d\n   bM:", _lbMarray);
  const double* bb;
  for (j = 0; j < _nRandom; j++){
    bb = _bM + j;
    Rprintf("\n     b%d = ", j);
    for (k = 0; k < _nCluster; k++){
      Rprintf("%g,  ", *bb);
      bb += _nRandom;
    }
  }  

  Rprintf("\n   larray = %d", _larray);
  Rprintf("\n   diagI = ");        for (j = 0; j < _nRandom; j++) Rprintf("%d,  ", _diagI[j]);
  Rprintf("\n   indUpd = ");       for (j = 0; j < _nRandom; j++) Rprintf("%d,  ", _indUpd[j]);
  Rprintf("\n   covpar = ");       for (j = 0; j < _larray; j++) Rprintf("%g,  ", _covpar[j]);
  Rprintf("\n   ichicovpar = ");   for (j = 0; j < _larray; j++) Rprintf("%g,  ", _ichicovpar[j]);
  Rprintf("\n   Digamma = ");      for (j = 0; j < _nRandom; j++) Rprintf("%g,  ", _Digamma[j]);
  Rprintf("\n   propMean = ");     for (j = 0; j < _nRandom; j++) Rprintf("%g,  ", _propMean[j]);
  Rprintf("\n   propMeanTemp = "); for (j = 0; j < _nRandom; j++) Rprintf("%g,  ", _propMeanTemp[j]);
  Rprintf("\n");
  return;
}


/***** RandomEff2initArray *****/
void
RandomEff::RandomEff2initArray(int* parmI,  double* parmD) const
{
  int i;

  /** Set up integer pointers **/
  int ptype_prior = 0;
  int pnRandom    = ptype_prior + 1;
  int pnCluster   = pnRandom + 1;
  int pnwithinCl  = pnCluster + 1;
  //  int pnext = pnwithinCl + parmI[pnCluster];

  /** Set up double pointers **/
  int pbM = 0;
  // int pnext = pbM + parmI[pnRandom]*parmI[pnCluster];

  parmI[ptype_prior] = _type_prior;
  parmI[pnRandom]    = _nRandom;
  parmI[pnCluster]   = _nCluster;
  for (i = 0; i < _nCluster; i++) parmI[pnwithinCl + i] = _nwithinCl[i];

  for (i = 0; i < _lbMarray; i++) parmD[pbM + i] = _bM[i];
  
  return;
}


//
// ***** GIBBSupdate *****
//
// Update all random effects using a Gibbs move
// Update also regression residuals
//
// nP[1] ...... number of observations
// XA ........... [nobs*_nbeta] full design matrix for all observations
//                CHANGE COMPARED TO bayessurvreg1:
//                * stored in row-major order, i.e. XA[0,     ...,_nbeta-1]   = covariates for the first observation
//                                                  XA[_nbeta,...,2*_nbeta-1] = covariates for the second observation etc.
// ZZtb ......... [nobs*_larray] z*z' (only triangles) for each observation
//                CHANGE COMPARED TO bayessurvreg1:
//                * not stored in double-indexed array
// gg ........... G-spline giving the distribution of observations
// mu ........... [gg->dim(), gg->length(j)] already computed knots of the G-spline
// rM ........... [nP] allocation labels taking values from {0,..., g_y->total_length() - 1}
// bg ........... object with regression coefficients and means of random effects
// Dcm .......... object with covariance matrix of random effects
//
void
RandomEff::GIBBSupdate(double* regresResM,    const int* nP,
                       const double* XA,      const double* ZZtb,
                       const Gspline* gg,     double** const mu,         const int* rM,
                       const BetaGamma* bbg,  const CovMatrix* Dcm)
{
  if (!_nRandom) return;

  static double invsigscale2[_max_dim];
  static int i, j, cl, nobs, rank;
  int nX = bbg->nbeta();

  /*** Compute invsigscale2 ***/
  for (j = 0; j < gg->dim(); j++) invsigscale2[j] = gg->invscale2(j) * gg->invsigma2(j);

  /*** Mean of full conditional distribution, part 1 (D^{-1}gamma) ***/
  if (bbg->randomIntcpt() && _nRandom == 1){
    _Digamma[0] = (Dcm->icovm(0)) * _Eb0_;
  }
  else{
    if (bbg->randomIntcpt() && _nRandom > 1){         // gamma[0] = Eb0, but it is not contained in beta vector
        // first column (= row) of D^{-1} times gamma (the first element of the sum is sometimes zero)
      _Digamma[0] = Dcm->icovm(0) * _Eb0_;
      for (i = 1; i < _nRandom; i++) _Digamma[0] += Dcm->icovm(i) * bbg->beta(bbg->indbinXA(i));

        // all other columns of D^{-1} times gamma 
      for (j = 1; j < _nRandom; j++){
        _Digamma[j] = Dcm->icovm(Dcm->diagI(j)) * bbg->beta(bbg->indbinXA(j));
        for (i = j + 1; i < _nRandom; i++) _Digamma[j] += Dcm->icovm(Dcm->diagI(j) + i - j) * bbg->beta(bbg->indbinXA(i));
        for (i = 1; i < j; i++)            _Digamma[j] += Dcm->icovm(Dcm->diagI(i) + j - i) * bbg->beta(bbg->indbinXA(i));
        _Digamma[j] += Dcm->icovm(Dcm->diagI(0) + j) * _Eb0_;
      }
    }
    else{                                  // no random intercept, gamma = subvector of betaM
      Mxa(_Digamma, bbg->betaP(), Dcm->icovmP(), bbg->indbinXP(), &nX, &_nRandom, Dcm->diagIP());
    }
  }

  /*** Loop over clusters ***/
  double* regRes = regresResM;
  double* bb = _bM;
  const double* Xm = XA;
  const double* ZZm = ZZtb;
  const int* rp = rM;

  for (cl = 0; cl < _nCluster; cl++){
    nobs = _nwithinCl[cl];

      /* Initialize mean of the full conditional distribution (stored in _propMeanTemp) by zeros */
    for (j = 0; j < _nRandom; j++) _propMeanTemp[j] = 0.0;

      /* Initialize inverse variance of the full conditional distribution (stored in _covpar) by zeros */
    for (j = 0; j < _larray; j++) _covpar[j] = 0.0;    

      /* Loop over observations in a given cluster */
    for (i = 0; i < nobs; i++){

        /* Compute y - alpha - x'beta - scale*mu, store it in regresResM */
      if (bbg->randomIntcpt()) *regRes += (*bb);
      for (j = bbg->randomIntcpt(); j < _nRandom; j++){
        *regRes += Xm[bbg->indbinXA(j)] * bb[j];
      }
      *regRes -= (gg->intcpt(0) + gg->scale(0)*mu[0][*rp]);

        /* Inverse variance of the full conditional distribution, part 1 (+ z * z') */
      for (j = 0; j < _larray; j++){
        _covpar[j] += ZZm[j];
      }

        /* Mean of the full conditional distribution, part 1 (+ z * v) */
      if (bbg->randomIntcpt()) _propMeanTemp[0] += (*regRes);
      for (j = bbg->randomIntcpt(); j < _nRandom; j++){
        _propMeanTemp[j] += Xm[bbg->indbinXA(j)] * (*regRes);
      }

      regRes++;
      Xm  += bbg->nbeta();
      ZZm += _larray;
      rp++;
    }  /** end of the loop over observations in a given cluster **/

    //if (cl < 4){
    //      Rprintf("\nCluster %d, nobs = %d:", cl, nobs);
    //      Rprintf("\n  sum(z*v) = %g, %g;", _propMeanTemp[0], _propMeanTemp[1]);
    //    }

      /* Mean of the full conditional distribution, part 2: times 1/(sigma2*scale2) + (D^{-1}gamma) */
    for (j = 0; j < _nRandom; j++){
      _propMeanTemp[j] *= invsigscale2[0];
      _propMeanTemp[j] += _Digamma[j];
    }

    //    if (cl < 4){
    //      Rprintf("\n  Digamma + 1/sigscale2*sum(z*v) = %g, %g;", _propMeanTemp[0], _propMeanTemp[1]);      
    //    }

      /* Inverse variance of the full conditional distribution, part 2: times  1/(sigma2*scale2) + (D^{-1}) */
    for (j = 0; j < _larray; j++){
      _covpar[j] *= invsigscale2[0];
      _covpar[j] += Dcm->icovm(j);
    }

      /* Cholesky decomposition of the inverse variance of full conditional distrib. */
    cholesky(_covpar, &rank, &_nRandom, Dcm->diagIP(), &_toler_chol_CovMatrix);

      /* Variance of the full conditional distribution                               */
      /*  and the inverse of the Cholesky decomp. of the inverse variance            */
    chinv2(_covpar, _ichicovpar, &_nRandom, Dcm->diagIP());

      /* Mean of full conditional distribution, part 3 (* var(b|...))                */
    Mxa(_propMean, _propMeanTemp, _covpar, &ZERO_INT, &_nRandom, &_nRandom, Dcm->diagIP());

      /* Sample                                                                      */
    rmvtnorm2(bb, _propMean, _ichicovpar, &ZERO_INT, _indUpd, &_nRandom, &_nRandom, &_nRandom, &ONE_INT, Dcm->diagIP(), &ZERO_INT);

    //    if (cl < 4){
    //      Rprintf("\nCluster %d:", cl);
    //      Rprintf("\n  propMean = %g, %g;", _propMean[0], _propMean[1]);
    //      Rprintf("\n  propVar = %g, %g, %g;", _covpar[0], _covpar[1], _covpar[2]);
    //      Rprintf("\n  ichicovpar = %g, %g, %g;", _ichicovpar[0], _ichicovpar[1], _ichicovpar[2]);
    //      Rprintf("\n  sampled b = %g, %g", bb[0], bb[1]);
    //    }

      /* Update regresResM                                                           */
    regRes -= nobs;
    Xm     -= nobs*bbg->nbeta();
    rp     -= nobs;
    for (i = 0; i < nobs; i++){
      if (bbg->randomIntcpt()) *regRes -= (*bb);
      for (j = bbg->randomIntcpt(); j < _nRandom; j++){
        *regRes -= Xm[bbg->indbinXA(j)] * bb[j];
      }
      *regRes += (gg->intcpt(0) + gg->scale(0)*mu[0][*rp]);

      regRes++;
      Xm  += bbg->nbeta();
      rp++;
    }

    bb += _nRandom;
  }   /*** end of the loop over clusters ***/

  return;
}    /** end of the function RandomEff::GIBBSupdate **/


//
// ***** sumSquare
//
// Compute sum_{i=1}^{nCluster}(b_i - Eb) * (b_i - Eb)'
// 
// sumSq ...... [_larray] array to store the result (only lower triangle)
// bbg ........ an object holding means of random effects
//
void
RandomEff::sumSquare(double* sumSq,  const BetaGamma* bbg) const
{
  static int i, j, cl;

  for (i = 0; i < _larray; i++) sumSq[i] = 0.0;

  const double* bb = _bM;
  for (cl = 0; cl < _nCluster; cl++){

    /** Compute b_i - Eb, store it in _propMean  **/
    for (j = 0; j < _nRandom; j++){
      _propMean[j] = bb[j] - (bbg->indbinXA(j) < 0 ? _Eb0_ : bbg->beta(bbg->indbinXA(j)));
    }

    /** Compute the sum of squares **/
    for (j = 0; j < _nRandom; j++){
      for (i = j; i < _nRandom; i++){
        sumSq[_diagI[j] + i - j] += _propMean[i] * _propMean[j];
      }
    }
    bb += _nRandom;
  }

  return;
}    /***  end of function RandomEff::sumSquare   ***/


/***** predictNormalRE *****/
// 
// Sample values of random effects under the assumption of normality
// given their mean and covariance matrix
//
// * this function replaces 'predictRandom' from 'predictMisc.cpp'
//
void
RandomEff::predictNormalRE(const BetaGamma* bg,  const CovMatrix* DD)
{
  if (!_nRandom) return;
 
  static int j, kk, cl;

  /*** Mean of the normal distribution to sample from ***/
  for (j = 0; j < _nRandom; j++){
    kk = bg->indbinXA(j);
    _propMean[j] = (kk < 0 ? _Eb0_ : bg->beta(kk));
  }

  /*** Loop over clusters, sample ***/
  double* bA = _bM;
  for (cl = 0; cl < _nCluster; cl++){
    rmvtnorm2(bA, _propMean, DD->ichicovmP(), &ZERO_INT, _indUpd, &_nRandom, &_nRandom, &_nRandom, &ONE_INT, DD->diagIP(), &ZERO_INT);
    bA += _nRandom;
  }
  return;
}


//
// ***** Gspl_intcpt_update *****
//
// This function assumes that _nRandom == 1
//
// Update random intercept under the assumption that it is distributed as a G-spline
// Update also regression residuals
//
// nP[1] ........ number of observations

// gg_b ......... G-spline giving the distribution of the random intercept
// mu_b ......... [gg_b->dim(), gg_b->length(j)] already computed knots of the random intercept G-spline
// rM_b ......... [nCluster] allocation labels for the random intercept

// gg ........... G-spline giving the distribution of observations
// mu ........... [gg->dim(), gg->length(j)] already computed knots of the observational G-spline
// rM ........... [nP] allocation labels taking values from {0,..., g_y->total_length() - 1} for observations
//
void
RandomEff::Gspl_intcpt_update(double* regresResM,    const int* nP,
                              const Gspline* gg_b,   double** const mu_b,       const int* rM_b,
                              const Gspline* gg,     double** const mu,         const int* rM)

{
  if (!_nRandom) return;

  static double invsigscale2[_max_dim];
  static double invsigscale2_b;
  static int i, j, cl;

  /*** Compute invsigscale2 and invsigscale2_b ***/
  for (j = 0; j < gg->dim(); j++) invsigscale2[j] = gg->invscale2(j) * gg->invsigma2(j);
  invsigscale2_b = gg_b->invscale2(0) * gg_b->invsigma2(0);

  /*** Loop over clusters ***/
  double* regRes = regresResM;
  double* bb = _bM;
  const int* rp_b = rM_b;
  const int* rp = rM;

  for (cl = 0; cl < _nCluster; cl++){

    /** Variance of the full conditional distribution **/
    _covpar[0] = invsigscale2_b + _nwithinCl[cl]*invsigscale2[0];
    if (_covpar[0] <= 0) throw returnR("Trap: Non-positive proposal variance for update of the random intercept", 1);
    _covpar[0] = 1/_covpar[0];

    /** Mean of the full conditional distribution, part 1 **/
    _propMean[0] = invsigscale2_b * (gg_b->intcpt(0) + gg_b->scale(0)*mu_b[0][*rp_b]);

    /** Loop over observations in a given cluster **/
    _propMeanTemp[0] = 0.0;
    for (i = 0; i < _nwithinCl[cl]; i++){
      
      /* Add old value of the random intercept to regRes                       */
      /* Compute sum(y - alpha - x'beta - scale*mu), store it in _propMeanTemp */
      *regRes += (*bb);
      _propMeanTemp[0] += (*regRes) - (gg->intcpt(0) + gg->scale(0)*mu[0][*rp]);

      regRes++;
      rp++;
    }       /** end of the loop over observations in a given cluster **/
    
    /** Mean of the full conditional distribution, finish **/
    _propMean[0] += invsigscale2[0] * _propMeanTemp[0];
    _propMean[0] *= _covpar[0];

    /** Standard deviation of the full conditional distribution **/
    _ichicovpar[0] = sqrt(_covpar[0]);
  
    /** Sample new value of the random intercept **/
    *bb = rnorm(_propMean[0], _ichicovpar[0]);

    /** Update regresResM  **/
    regRes -= _nwithinCl[cl];
    for (i = 0; i < _nwithinCl[cl]; i++){
      *regRes -= (*bb);
      regRes++;
    }
    bb++;
    rp_b++;
  }    /*** end of the loop over clusters ***/

  return;
}


//
// ***** adjust_intcpt *****
//
// From all values of _bM subtract the value of *adj
//
// * used in the model with G-spline random intercept
// * called by the function 'adjust_intcpt' defined in 'bayessurvreg2.cpp'
//
void
RandomEff::adjust_intcpt(const double* adj)
{
  if (_nRandom != 1) throw returnR("C++ Error: RandomEff::adjust_intcpt can be called only for univariate random effects", 1);

  static int cl;
  double* bb = _bM;
  for (cl = 0; cl < _nCluster; cl++){
    (*bb) -= (*adj);
    bb++;
  }

  return;
}


//
// ***** predictGspl_intcpt ******
//
// Sample random intercept from the G-spline distribution
// * used by predictive functions
//
// k_effect ........ effective number of mixture components
// cum_w ........... cumulative weights
// prop_mu ......... proposal means for each mixture component (intcpt + scale*knot)
// sig_scale ....... proposal standard deviation (sigma*scale)
// rM_b ............ [gg_b->nCluster()] working space to sample allocations
//
void
RandomEff::predictGspl_intcpt(const int* k_effect,      double* cum_w,  const double* prop_mu, 
                              const double* sig_scale,  int* rM_b)
{
  if (!_nRandom) return;  

  static int cl;

  /** Sample allocations   **/
  discreteSampler2(rM_b, cum_w, k_effect, &_nCluster, &ONE_INT, &ZERO_INT);

  /** Sample values of the random intercept  **/
  double* bb = _bM;
  int* rp_b = rM_b;
  for (cl = 0; cl < _nCluster; cl++){
    *bb = rnorm(prop_mu[*rp_b], sig_scale[0]);
    bb++;
    rp_b++;
  }

  return;
}
