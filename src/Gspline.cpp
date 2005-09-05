// Class to store information on a G-spline
//  and to update it
//
// 07/10/2004: start working
// 20/02/2005: 'adjust_intcpt' added (immediately decided not to use it)
//
#include "Gspline.h"

using namespace std;

/***** Nonparametric constructor *****/
Gspline::Gspline()
  : _dim(0), _neighbor_system(0), _equal_lambda(true), _total_length(0), _length(NULL), _K(NULL), _total_izero(0), _izero(NULL), _order(0),
    _log_null_w(0.0),
    _lambda(NULL), _a(NULL), 
    _a_max(0.0),   _expa(NULL), _sumexpa(0.0), _sumexpa_margin(NULL), _penalty(NULL), _k_effect(0), _ind_w_effect(NULL),
    _abscis(NULL), _iwv(NULL), _rwv(NULL), _hx(NULL), _hpx(NULL), _type_update_a(0), _k_overrelax_a(1),
    _gamma(NULL), _invsigma2(NULL), _sigma(NULL), _c4delta(NULL), _delta(NULL),
    _intcpt(NULL), _invscale2(NULL), _scale(NULL)
{
  int j;
  for (j = 0; j < _max_dim; j++){
     _prior_for_lambda[j] = 0;
     _prior_for_gamma[j] = 0;
     _prior_for_sigma[j] = 0;
     _prior_for_intcpt[j] = 0;
     _prior_for_scale[j] = 0;
     _k_overrelax_sigma[j] = 1;
     _k_overrelax_scale[j] = 1;
  }
  for (j = 0; j < 2* _max_dim; j++){
    _prior_lambda[j] = 0.0;
    _prior_gamma[j] = 0.0;
    _prior_sigma[j] = 0.0;
    _prior_intcpt[j] = 0.0;
    _prior_scale[j] = 0.0;
  }
}

/***** Parametric constructor *****/
//
// pprior_lambda ....... pairs (shape, rate) for each gamma distribution 
//                       or pairs (S, whatever) if uniform (0, S) is used for lambda^{-1/2}
// pprior_gamma ........ pairs (mean, variance) for each normal distribution
//
Gspline::Gspline(const int* parmI,   const double* parmD)
{
  int i, j, k;

  /** Set up integer pointers **/
  int ddim               = 0;
  int nneighbor          = ddim               + 1;
  int eequal_lambda      = nneighbor          + 1;
  int KK                 = eequal_lambda      + 1;
  int iizero             = KK                 + parmI[ddim];        /* on INPUT: scale -K, ..., K; in the G-spline: scale 0, ..., 2*K */
  int oorder             = iizero             + parmI[ddim];
  int pprior_for_lambda  = oorder             + 1;
  int pprior_for_gamma   = pprior_for_lambda  + parmI[ddim];
  int pprior_for_sigma   = pprior_for_gamma   + parmI[ddim];
  int pprior_for_intcpt  = pprior_for_sigma   + parmI[ddim];
  int pprior_for_scale   = pprior_for_intcpt  + parmI[ddim];
  int ttype_update_a     = pprior_for_scale   + parmI[ddim];
  int kk_overrelax_a     = ttype_update_a     + 1;
  int kk_overrelax_sigma = kk_overrelax_a     + 1;
  int kk_overrelax_scale = kk_overrelax_sigma + parmI[ddim];
  // int inext = kk_overrelax_scale + parmI[ddim];

  if (parmI[ddim] < 0 || parmI[ddim] > 2) throw returnR("C++ Error: G-spline of incorrect/unimplemented dimension asked", 1);
  _dim = parmI[ddim];

  if (_dim == 0){
    _neighbor_system = 0;  _equal_lambda = true;
    _total_length = 0;     _total_izero = 0;
    _K = NULL;             _length = NULL;      _izero = NULL;
    _order = 0;            _log_null_w = 0.0;          
    _lambda = NULL;        _a = NULL;           _a_max = 0.0;
    _expa = NULL;          _sumexpa = 0.0;      _sumexpa_margin = NULL, 
    _penalty = NULL,       _k_effect = 0;       _ind_w_effect = NULL;
    _abscis = NULL;        _iwv = NULL;         _rwv = NULL;            _hx = NULL;         _hpx = NULL;
    _type_update_a = 0;    _k_overrelax_a = 1;
    _gamma = NULL;         _invsigma2 = NULL;   _sigma = NULL;          _c4delta = NULL;    _delta = NULL;
    _intcpt = NULL;        _invscale2 = NULL;   _scale = NULL;
    for (j = 0; j < _max_dim; j++){
      _prior_for_lambda[j] = 0;
      _prior_for_gamma[j] = 0;
      _prior_for_sigma[j] = 0;
      _prior_for_intcpt[j] = 0;
      _prior_for_scale[j] = 0;
      _k_overrelax_sigma[j] = 1;
      _k_overrelax_scale[j] = 1;
    }
    for (j = 0; j < 2* _max_dim; j++){
      _prior_lambda[j] = 0.0;
      _prior_gamma[j] = 0.0;
      _prior_sigma[j] = 0.0;
      _prior_intcpt[j] = 0.0;
      _prior_scale[j] = 0.0;
    }
  }
  else{
    if (_dim == 1){
      _neighbor_system = uniCAR;
      _equal_lambda = true;
      _order = parmI[oorder];
    }
    else{
      switch (parmI[nneighbor]){
      case uniCAR:           
        _neighbor_system = uniCAR;
        _equal_lambda = parmI[eequal_lambda];
        _order = parmI[oorder];
        break;
      case eight_neighbors:  
        _neighbor_system = eight_neighbors;
	_equal_lambda = true;
        _order = 2;       /* two here to update optimally 'a' (with as small as possible autocorrelation) */
        break;
      case twelve_neighbors: 
        _neighbor_system = twelve_neighbors;
	_equal_lambda = true;
        _order = 3;       /* three here to update optimally 'a' (with as small as possible autocorrelation) */
        break;
      default:               
        throw returnR("C++ Error: Unimplemented neighboring system of the G-spline supplied", 1);
      }
    }  /* end of _dim > 1 */
    if (_order < 0 || _order > 3) throw returnR("C++ Error: Unimplemented order of autoregression of the G-spline", 1);   // order == 0 -> fixed a's

    _length = new int[_dim];
    if (_length == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _K = new int[_dim];
    if (_K == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _izero = new int[_dim];
    if (_izero == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      if (parmI[KK+i] < 0) throw returnR("C++ Error: G-spline must have a positive length (i.e. K >= 0)", 1);
      if (parmI[KK+i] < _order) throw returnR("C++ Error: All _K's in the G-spline must be at least equal to_order", 1);
      _K[i] = parmI[KK+i];
      _length[i] = 2*_K[i] + 1;
      if (parmI[iizero+i] < -_K[i] || parmI[iizero+i] > _K[i]) throw returnR("C++ Error: Supplied izero for the G-spline out of range", 1);     
      _izero[i] = parmI[iizero+i] + _K[i];     /* change R index (-K,...,K) to G-spline index (0,...,2K) */
    }
    
    switch (_dim){
    case 1:
      _total_length = _length[0];
      _total_izero = _izero[0];
      break;
    case 2:
      if (_length[0] <= 0 || _length[1] <= 0) throw returnR("C++ Error: G-spline must have positive lengths", 1);
      _total_length = _length[0] * _length[1];
      _total_izero = _izero[1]*_length[0] + _izero[0];
      break;
    default:
      throw returnR("C++ Error: Unimplemented dimension of the G-spline", 1);
    }

    _log_null_w = log(_null_mass) - log(1 - _null_mass) - log(double(_total_length));

    /** Set up double pointers **/
    int aa             = 0;
    int llambda        = aa            + _total_length;
    int ggamma         = llambda       + _dim;
    int ssigma         = ggamma        + _dim;
    int iintcpt        = ssigma        + _dim;
    int sscale         = iintcpt       + _dim;
    int cc4delta       = sscale        + _dim;
    int pprior_lambda  = cc4delta      + _dim;
    int pprior_gamma   = pprior_lambda + 2*_dim;
    int pprior_sigma   = pprior_gamma  + 2*_dim;
    int pprior_intcpt  = pprior_sigma  + 2*_dim;
    int pprior_scale   = pprior_intcpt + 2*_dim;
    // int dnext = pprior_scale + 2*_dim;

  /*** Storage of G-spline variables ***/
    _a = new double[_total_length];
    if (_a == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _expa = new double[_total_length];
    if (_expa == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _ind_w_effect = new int[_total_length];
    if (_ind_w_effect == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++) _a[i] = parmD[aa + i];
    if (_dim == 1) _sumexpa_margin = NULL;
    else{
      _sumexpa_margin = new double*[_dim];
      if (_sumexpa_margin == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      for (j = 0; j < _dim; j++){
        _sumexpa_margin[j] = new double[_length[j]];
        if (_sumexpa_margin[j] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      }
    }
    a2expa();
    update_a_max_center_and_k_effect();

    _penalty = new double[_dim];
    if (_penalty == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    if (_equal_lambda && _dim > 1) for (j = 1; j < _dim; j++) _penalty[j] = 0.0;
    penalty();

    _lambda = new double[_dim];
    if (_lambda == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _gamma = new double[_dim];
    if (_gamma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _sigma = new double[_dim];
    if (_sigma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invsigma2 = new double[_dim];
    if (_invsigma2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _c4delta = new double[_dim];
    if (_c4delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _delta = new double[_dim];
    if (_delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _intcpt = new double[_dim];
    if (_intcpt == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _scale = new double[_dim];
    if (_scale == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invscale2 = new double[_dim];
    if (_invscale2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      if (!_equal_lambda || i == 0){
        if (parmD[llambda + i] <= 0) throw returnR("C++ Error: Non-positive lambda for the G-spline supplied", 1);
        _lambda[i] = parmD[llambda + i];
      }
      else
        _lambda[i] = _lambda[0];

      _gamma[i] = parmD[ggamma + i];
      if (parmD[ssigma + i] <= 0) throw returnR("C++ Error: Non-positive standard deviation of the G-spline supplied", 1);
      _sigma[i] = parmD[ssigma + i];
      _invsigma2[i] = 1 / (_sigma[i]*_sigma[i]);

      _intcpt[i] = parmD[iintcpt + i];
      if (parmD[sscale + i] <= 0) throw returnR("C++ Error: Non-positive scale of the G-spline supplied", 1);
      _scale[i] = parmD[sscale + i];
      _invscale2[i] = 1 / (_scale[i]*_scale[i]);

      if (parmD[cc4delta + i] <= 0) throw returnR("C++ Error: Non-positive 'c' coefficient for the G-spline supplied", 1);
      _c4delta[i] = parmD[cc4delta + i];
      _delta[i] = _c4delta[i] * _sigma[i];
    }

  /*** Specifications of priors ***/
    /** lambda **/    
    for (k = 0; k <= (_dim-1)*(!_equal_lambda); k++){
      switch (parmI[pprior_for_lambda + k]){
      case Fixed_:
        _prior_for_lambda[k] = Fixed_;
        for (i = 0; i < 2; i++) _prior_lambda[2*k+i] = 0.0;
        break;
      case Gamma:  
        _prior_for_lambda[k] = Gamma; 
        if (parmD[pprior_lambda + 2*k] < 0)   throw returnR("C++ Error: Shape parameter for G-spline lambda prior must be positive", 1);
        if (parmD[pprior_lambda + 2*k+1] < 0) throw returnR("C++ Error: Rate parameter for G-spline lambda prior must be positive", 1);
        for (i = 0; i < 2; i++) _prior_lambda[2*k+i] = parmD[pprior_lambda + 2*k+i];
        break;
      case SDUnif: 
        _prior_for_lambda[k] = SDUnif; 
        if (parmD[pprior_lambda + 2*k] < 0) throw returnR("C++ Error: Upper limit for G-spline lambda^{-1/2} prior must be positive", 1);
        _prior_lambda[2*k] = 0.0;
        _prior_lambda[2*k+1] = 1/(parmD[pprior_lambda + 2*k+1]*parmD[pprior_lambda + 2*k+1]);    /* store 1/S^2 where S is the upper limit for lambda^{-1/2} */
        break;
      default:     
        throw returnR("C++ Error: Unimplemented prior for lambda of the G-spline asked", 1);
      }
    }
    for (k = 1*(_equal_lambda) + _dim*(!_equal_lambda); k < _max_dim; k++){
      _prior_for_lambda[k] = Fixed_;
      for (i = 0; i < 2; i++) _prior_lambda[2*k+i] = 0.0;
    }

    /** gamma **/
    for (k = 0; k < _dim; k++){
      switch (parmI[pprior_for_gamma + k]){
      case _Fixed_:
        _prior_for_gamma[k] = _Fixed_;
        for (i = 0; i < 2; i++) _prior_gamma[2*k+i] = 0.0;
        break;
      case Normal:
        _prior_for_gamma[k] = Normal;
        _prior_gamma[2*k] = parmD[pprior_gamma + 2*k];
        if (parmD[pprior_gamma + 2*k+1] <= 0) throw returnR("C++ Error: Non-positive prior variance of gamma G-spline parameter supplied", 1);
        _prior_gamma[2*k+1] = 1/parmD[pprior_gamma + 2*k+1];       /* store inverted variance */
        break;
      default:
        throw returnR("C++ Error: Unimplemented prior for sigma of the G-spline asked", 1);
      }
    }
    for (k = _dim; k < _max_dim; k++){
      _prior_for_gamma[k] = _Fixed_;
      for (i = 0; i < 2; i++) _prior_gamma[2*k+i] = 0.0;
    }

    /** sigma **/
    for (k = 0; k < _dim; k++){
      switch (parmI[pprior_for_sigma + k]){
      case Fixed_:
        _prior_for_sigma[k] = Fixed_;
        for (i = 0; i < 2; i++) _prior_sigma[2*k+i] = 0.0;
        _k_overrelax_sigma[k] = 1;
        break;
      case Gamma:  
        _prior_for_sigma[k] = Gamma; 
        if (parmD[pprior_sigma + 2*k] < 0)   throw returnR("C++ Error: Shape parameter for G-spline sigma{-2} prior must be positive", 1);
        if (parmD[pprior_sigma + 2*k+1] < 0) throw returnR("C++ Error: Rate parameter for G-spline sigma^{-2} prior must be positive", 1);
        for (i = 0; i < 2; i++) _prior_sigma[2*k+i] = parmD[pprior_sigma + 2*k+i];       
        if (parmI[kk_overrelax_sigma+k] <= 0) throw returnR("C++ Error: _k_overrelax_sigma must be positive in a constructor of Gspline", 1);
        _k_overrelax_sigma[k] = parmI[kk_overrelax_sigma+k];
        break;
      case SDUnif: 
        _prior_for_sigma[k] = SDUnif; 
        if (parmD[pprior_sigma + 2*k+1] < 0) throw returnR("C++ Error: Upper limit for G-spline sigma prior must be positive", 1);
        _prior_sigma[2*k] = 0.0;
        _prior_sigma[2*k+1] = 1/(parmD[pprior_sigma + 2*k+1]*parmD[pprior_sigma + 2*k+1]);    /* store 1/S^2 where S is the upper limit for sigma */
        if (parmI[kk_overrelax_sigma+k] <= 0) throw returnR("C++ Error: _k_overrelax_sigma must be positive in a constructor of Gspline", 1);
        _k_overrelax_sigma[k] = parmI[kk_overrelax_sigma+k];
        break;
      default:     
        throw returnR("C++ Error: Unimplemented prior for sigma of the G-spline asked", 1);
      }
    }        
    for (k = _dim; k < _max_dim; k++){
      _prior_for_sigma[k] = Fixed_;
      _k_overrelax_sigma[k] = 1;
      for (i = 0; i < 2; i++) _prior_sigma[2*k+i] = 0.0;
    }

    /** intcpt **/
    for (k = 0; k < _dim; k++){
      switch (parmI[pprior_for_intcpt + k]){
      case _Fixed_:
        _prior_for_intcpt[k] = _Fixed_;
        for (i = 0; i < 2; i++) _prior_intcpt[2*k+i] = 0.0;
        break;
      case Normal:
        _prior_for_intcpt[k] = Normal;
        _prior_intcpt[2*k] = parmD[pprior_intcpt + 2*k];
        if (parmD[pprior_intcpt + 2*k+1] <= 0) throw returnR("C++ Error: Non-positive prior variance of intcpt G-spline parameter supplied", 1);
        _prior_intcpt[2*k+1] = 1/parmD[pprior_intcpt + 2*k+1];       /* store inverted variance */
        break;
      default:
        throw returnR("C++ Error: Unimplemented prior for sigma of the G-spline asked", 1);
      }
    }
    for (k = _dim; k < _max_dim; k++){
      _prior_for_intcpt[k] = _Fixed_;
      for (i = 0; i < 2; i++) _prior_intcpt[2*k+i] = 0.0;
    }

    /** scale **/
    for (k = 0; k < _dim; k++){
      switch (parmI[pprior_for_scale + k]){
      case Fixed_:
        _prior_for_scale[k] = Fixed_;
        for (i = 0; i < 2; i++) _prior_scale[2*k+i] = 0.0;
        _k_overrelax_scale[k] = 1;
        break;
      case Gamma:  
        _prior_for_scale[k] = Gamma; 
        if (parmD[pprior_scale + 2*k] < 0)   throw returnR("C++ Error: Shape parameter for G-spline scale{-2} prior must be positive", 1);
        if (parmD[pprior_scale + 2*k+1] < 0) throw returnR("C++ Error: Rate parameter for G-spline scale^{-2} prior must be positive", 1);
        for (i = 0; i < 2; i++) _prior_scale[2*k+i] = parmD[pprior_scale + 2*k+i];       
        if (parmI[kk_overrelax_scale+k] <= 0) throw returnR("C++ Error: _k_overrelax_scale must be positive in a constructor of Gspline", 1);
        _k_overrelax_scale[k] = parmI[kk_overrelax_scale+k];
        break;
      case SDUnif: 
        _prior_for_scale[k] = SDUnif; 
        if (parmD[pprior_scale + 2*k+1] < 0) throw returnR("C++ Error: Upper limit for G-spline scale prior must be positive", 1);
        _prior_scale[2*k] = 0.0;
        _prior_scale[2*k+1] = 1/(parmD[pprior_scale + 2*k+1]*parmD[pprior_scale + 2*k+1]);    /* store 1/S^2 where S is the upper limit for scale */
        if (parmI[kk_overrelax_scale+k] <= 0) throw returnR("C++ Error: _k_overrelax_scale must be positive in a constructor of Gspline", 1);
        _k_overrelax_scale[k] = parmI[kk_overrelax_scale+k];
        break;
      default:     
        throw returnR("C++ Error: Unimplemented prior for scale of the G-spline asked", 1);
      }
    }        
    for (k = _dim; k < _max_dim; k++){
      _prior_for_scale[k] = Fixed_;
      _k_overrelax_scale[k] = 1;
      for (i = 0; i < 2; i++) _prior_scale[2*k+i] = 0.0;
    }


  /*** Stuff for adaptive rejection/slice  sampling of 'a' coefficients ***/
    _abscis = new double*[_total_length];
    if (_abscis == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++){
      _abscis[i] = new double[_nabscis];
      if (_abscis[i] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      find_start_abscis(&i);
    }

    _iwv = new int[7 + _ns];
    if (_iwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _rwv = new double[9 + 6*(1+_ns)];
    if (_rwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _hx = new double[_nabscis];
    if (_hx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _hpx = new double[_nabscis];
    if (_hpx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);

    switch (parmI[ttype_update_a]){
    case Slice:
      _type_update_a = parmI[ttype_update_a];
      if (parmI[kk_overrelax_a] <= 0) throw returnR("C++ Error: _k_overrelax_a must be positive in a constructor of Gspline", 1);
      _k_overrelax_a = parmI[kk_overrelax_a];
      break;
    case ARS_quantile:
    case ARS_mode:
      _type_update_a = parmI[ttype_update_a];
      _k_overrelax_a = 1;
      break;
    default:
      throw returnR("C++ Error: Unimplemented _type_update_a appeared in a constructor of Gspline", 1);
    }
  }  /** end of else (_dim == 0) **/

}  /** end of the parametric constructor **/


/***** Copy constructor *****/
Gspline::Gspline(const Gspline& gg)
{
  int i, j, k;

  if (gg.dim() == 0){
    _neighbor_system = 0;  _equal_lambda = true;
    _total_length = 0;     _total_izero = 0;
    _K = NULL;             _length = NULL;      _izero = NULL;
    _order = 0;            _log_null_w = 0.0;
    _lambda = NULL;        _a = NULL;           _a_max = 0.0;
    _expa = NULL;          _sumexpa = 0.0;      _sumexpa_margin = NULL, _penalty = NULL,    _k_effect = 0;      _ind_w_effect = NULL;
    _abscis = NULL;        _iwv = NULL;         _rwv = NULL;            _hx = NULL;         _hpx = NULL;
    _type_update_a = 0;    _k_overrelax_a = 1;
    _gamma = NULL;         _invsigma2 = NULL;   _sigma = NULL;          _c4delta = NULL;    _delta = NULL;
    _intcpt = NULL;        _invscale2 = NULL;   _scale = NULL;
    for (j = 0; j < _max_dim; j++){
      _prior_for_lambda[j] = 0;
      _prior_for_gamma[j] = 0;
      _prior_for_sigma[j] = 0;
      _prior_for_intcpt[j] = 0;
      _prior_for_scale[j] = 0;
      _k_overrelax_sigma[j] = 1;
      _k_overrelax_scale[j] = 1;
    }
    for (j = 0; j < 2* _max_dim; j++){
      _prior_lambda[j] = 0.0;
      _prior_gamma[j] = 0.0;
      _prior_sigma[j] = 0.0;
      _prior_intcpt[j] = 0.0;
      _prior_scale[j] = 0.0;
    }
  }
  else{
    _dim = gg.dim();
    _neighbor_system = gg._neighbor_system;
    _equal_lambda = gg._equal_lambda;

    _total_length = gg.total_length();
    _total_izero = gg.total_izero();
    _length = new int[_dim];
    if (_length == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _K = new int[_dim];
    if (_K == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _izero = new int[_dim];
    if (_izero == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      _K[i] = gg.K(i);
      _length[i] = gg.length(i);
      _izero[i] = gg.izero(i);
    }

    _order = gg.order();
    _log_null_w = gg.log_null_w();

    _a = new double[_total_length];
    if (_a == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _expa = new double[_total_length];
    if (_expa == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _ind_w_effect = new int[_total_length];
    if (_ind_w_effect == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++){
      _a[i] = gg.a(i);
      _expa[i] = gg.expa(i);    
    }
    _a_max = gg.a_max();
    _sumexpa = gg.sumexpa();
    _k_effect = gg.k_effect();
    for (i = 0; i < _k_effect; i++) _ind_w_effect[i] = gg.ind_w_effect(i);
    if (_dim == 1) _sumexpa_margin = NULL;
    else{
      _sumexpa_margin = new double*[_dim];
      if (_sumexpa_margin == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      for (j = 0; j < _dim; j++){
        _sumexpa_margin[j] = new double[_length[j]];
        if (_sumexpa_margin[j] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
        for (k = 0; k < _length[j]; k++) _sumexpa_margin[j][k] = gg.sumexpa_margin(j, k);
      }
    }

    _abscis = new double*[_total_length];
    if (_abscis == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++){
      _abscis[i] = new double[_nabscis];
      if (_abscis[i] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      for (k = 0; k < _nabscis; k++) _abscis[i][k] = gg.abscis(i, k);
    }

    _iwv = new int[7 + _ns];
    if (_iwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < 7 + _ns; i++) _iwv[i] = gg._iwv[i];
    _rwv = new double[9 + 6*(1+_ns)];
    if (_rwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < 9 + 6*(1+_ns); i++) _rwv[i] = gg._rwv[i];
    _hx = new double[_nabscis];
    if (_hx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _hpx = new double[_nabscis];
    if (_hpx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _nabscis; i++){
      _hx[i] = gg._hx[i];
      _hpx[i] = gg._hpx[i];
    }
    _type_update_a = gg._type_update_a;
    _k_overrelax_a = gg._k_overrelax_a;

    _penalty = new double[_dim];
    if (_penalty == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _lambda = new double[_dim];
    if (_lambda == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _gamma = new double[_dim];
    if (_gamma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _sigma = new double[_dim];
    if (_sigma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invsigma2 = new double[_dim];
    if (_invsigma2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _c4delta = new double[_dim];
    if (_c4delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _delta = new double[_dim];
    if (_delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _intcpt = new double[_dim];
    if (_intcpt == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _scale = new double[_dim];
    if (_scale == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invscale2 = new double[_dim];
    if (_invscale2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      _penalty[i] = gg.penalty(i);
      _lambda[i] = gg.lambda(i);
      _gamma[i] = gg.gamma(i);
      _sigma[i] = gg.sigma(i);
      _invsigma2[i] = gg.invsigma2(i);
      _c4delta[i] = gg.c4delta(i);
      _delta[i] = gg.delta(i);
      _intcpt[i] = gg.intcpt(i);
      _scale[i] = gg.scale(i);
      _invscale2[i] = gg.invscale2(i);
    }
    for (j = 0; j < _max_dim; j++){
      _prior_for_lambda[j] = gg._prior_for_lambda[j];
      _prior_for_gamma[j] = gg._prior_for_gamma[j];
      _prior_for_sigma[j] = gg._prior_for_sigma[j];
      _prior_for_intcpt[j] = gg._prior_for_intcpt[j];
      _prior_for_scale[j] = gg._prior_for_scale[j];
      _k_overrelax_sigma[j] = gg._k_overrelax_sigma[j];
      _k_overrelax_scale[j] = gg._k_overrelax_scale[j];
    }
    for (j = 0; j < 2* _max_dim; j++){
      _prior_lambda[j] = gg._prior_lambda[j];
      _prior_gamma[j] = gg._prior_gamma[j];
      _prior_sigma[j] = gg._prior_sigma[j];
      _prior_intcpt[j] = gg._prior_intcpt[j];
      _prior_scale[j] = gg._prior_scale[j];
    }
  }
}  /** end of the copy constructor **/


/***** Assignment operator *****/
Gspline&
Gspline::operator=(const Gspline& gg)
{
  int i, j, k;

  if (_dim > 0){
    delete [] _K,           delete [] _length;     delete [] _izero;
    delete [] _lambda;      delete [] _a;          
    delete [] _expa;        delete [] _penalty;    delete [] _ind_w_effect;
    for (i = 0; i < _total_length; i++) delete [] _abscis[i];
    delete [] _abscis;
    delete [] _iwv;         delete [] _rwv;        delete [] _hx;      delete [] _hpx;
    delete [] _gamma;       delete [] _invsigma2;  delete [] _sigma;   delete [] _c4delta;   delete [] _delta;
    delete [] _intcpt;      delete [] _invscale2;  delete [] _scale;
    if (_dim > 1){
      for (j = 0; j < _dim; j++) delete [] _sumexpa_margin[j];
      delete [] _sumexpa_margin;
    }
  }

  if (gg.dim() == 0){
    _neighbor_system = 0;  _equal_lambda = true;
    _total_length = 0;     _total_izero = 0;
    _K = NULL;             _length = NULL;      _izero = NULL;
    _order = 0;            _log_null_w = 0.0;
    _lambda = NULL;        _a = NULL;           _a_max = 0.0;
    _expa = NULL;          _sumexpa = 0.0;      _sumexpa_margin = NULL, _penalty = NULL;    _k_effect = 0;      _ind_w_effect = NULL;
    _abscis = NULL;        _iwv = NULL;         _rwv = NULL;            _hx = NULL;         _hpx = NULL;
    _type_update_a = 0;    _k_overrelax_a = 1;
    _gamma = NULL;         _invsigma2 = NULL;   _sigma = NULL;          _c4delta = NULL;    _delta = NULL;
    _intcpt = NULL;        _invscale2 = NULL;   _scale = NULL;
    for (j = 0; j < _max_dim; j++){
      _prior_for_lambda[j] = 0;
      _prior_for_gamma[j] = 0;
      _prior_for_sigma[j] = 0;
      _prior_for_intcpt[j] = 0;
      _prior_for_scale[j] = 0;
      _k_overrelax_sigma[j] = 1;
      _k_overrelax_scale[j] = 1;
    }
    for (j = 0; j < 2* _max_dim; j++){
      _prior_lambda[j] = 0.0;
      _prior_gamma[j] = 0.0;
      _prior_sigma[j] = 0.0;
      _prior_intcpt[j] = 0.0;
      _prior_scale[j] = 0.0;
    }
  }
  else{
    _dim = gg.dim();
    _neighbor_system = gg._neighbor_system;
    _equal_lambda = gg._equal_lambda;

    _total_length = gg.total_length();
    _total_izero = gg.total_izero();
    _length = new int[_dim];
    if (_length == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _K = new int[_dim];
    if (_K == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _izero = new int[_dim];
    if (_izero == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      _K[i] = gg.K(i);
      _length[i] = gg.length(i);
      _izero[i] = gg.izero(i);
    }

    _order = gg.order();
    _log_null_w = gg.log_null_w();

    _a = new double[_total_length];
    if (_a == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _expa = new double[_total_length];
    if (_expa == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _ind_w_effect = new int[_total_length];
    if (_ind_w_effect == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++){
      _a[i] = gg.a(i);
      _expa[i] = gg.expa(i);    
    }
    _a_max = gg.a_max();
    _sumexpa = gg.sumexpa();
    _k_effect = gg.k_effect();
    for (i = 0; i < _k_effect; i++) _ind_w_effect[i] = gg.ind_w_effect(i);
    if (_dim == 1) _sumexpa_margin = NULL;
    else{
      _sumexpa_margin = new double*[_dim];
      if (_sumexpa_margin == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      for (j = 0; j < _dim; j++){
        _sumexpa_margin[j] = new double[_length[j]];
        if (_sumexpa_margin[j] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
        for (k = 0; k < _length[j]; k++) _sumexpa_margin[j][k] = gg.sumexpa_margin(j, k);
      }
    }

    _abscis = new double*[_total_length];
    if (_abscis == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _total_length; i++){
      _abscis[i] = new double[_nabscis];
      if (_abscis[i] == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
      for (k = 0; k < _nabscis; k++) _abscis[i][k] = gg.abscis(i, k);
    }

    _iwv = new int[7 + _ns];
    if (_iwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < 7 + _ns; i++) _iwv[i] = gg._iwv[i];
    _rwv = new double[9 + 6*(1+_ns)];
    if (_rwv == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < 9 + 6*(1+_ns); i++) _rwv[i] = gg._rwv[i];
    _hx = new double[_nabscis];
    if (_hx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _hpx = new double[_nabscis];
    if (_hpx == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _nabscis; i++){
      _hx[i] = gg._hx[i];
      _hpx[i] = gg._hpx[i];
    }
    _type_update_a = gg._type_update_a;
    _k_overrelax_a = gg._k_overrelax_a;

    _penalty = new double[_dim];
    if (_penalty == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _lambda = new double[_dim];
    if (_lambda == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _gamma = new double[_dim];
    if (_gamma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _sigma = new double[_dim];
    if (_sigma == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invsigma2 = new double[_dim];
    if (_invsigma2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _c4delta = new double[_dim];
    if (_c4delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _delta = new double[_dim];
    if (_delta == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _intcpt = new double[_dim];
    if (_intcpt == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _scale = new double[_dim];
    if (_scale == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    _invscale2 = new double[_dim];
    if (_invscale2 == NULL) throw returnR("C++ Error: Could not allocate needed memory", 1);
    for (i = 0; i < _dim; i++){
      _penalty[i] = gg.penalty(i);
      _lambda[i] = gg.lambda(i);
      _gamma[i] = gg.gamma(i);
      _sigma[i] = gg.sigma(i);
      _invsigma2[i] = gg.invsigma2(i);
      _c4delta[i] = gg.c4delta(i);
      _delta[i] = gg.delta(i);
      _intcpt[i] = gg.intcpt(i);
      _scale[i] = gg.scale(i);
      _invscale2[i] = gg.invscale2(i);
    }
    for (j = 0; j < _max_dim; j++){
      _prior_for_lambda[j] = gg._prior_for_lambda[j];
      _prior_for_gamma[j] = gg._prior_for_gamma[j];
      _prior_for_sigma[j] = gg._prior_for_sigma[j];
      _prior_for_intcpt[j] = gg._prior_for_intcpt[j];
      _prior_for_scale[j] = gg._prior_for_scale[j];
      _k_overrelax_sigma[j] = gg._k_overrelax_sigma[j];
      _k_overrelax_scale[j] = gg._k_overrelax_scale[j];
    }
    for (j = 0; j < 2* _max_dim; j++){
      _prior_lambda[j] = gg._prior_lambda[j];
      _prior_gamma[j] = gg._prior_gamma[j];
      _prior_sigma[j] = gg._prior_sigma[j];
      _prior_intcpt[j] = gg._prior_intcpt[j];
      _prior_scale[j] = gg._prior_scale[j];
    }
  }
  return *this;
}  /** end of the assignment operator **/


/***** Destructor *****/
Gspline::~Gspline()
{
  if (_dim > 0){
    int i, j;
    delete [] _K;           delete [] _length;     delete [] _izero;
    delete [] _lambda;      delete [] _a;          
    delete [] _expa;        delete [] _penalty;    delete [] _ind_w_effect;
    for (i = 0; i < _total_length; i++) delete [] _abscis[i];
    delete [] _abscis;
    delete [] _iwv;         delete [] _rwv;        delete [] _hx;      delete [] _hpx;
    delete [] _gamma;       delete [] _invsigma2;  delete [] _sigma;   delete [] _c4delta;   delete [] _delta;
    delete [] _intcpt;      delete [] _invscale2;  delete [] _scale;
    if (_dim > 1){
      for (j = 0; j < _dim; j++) delete [] _sumexpa_margin[j];
      delete [] _sumexpa_margin;
    }
  }
}


/***** print *****/
void
Gspline::print() const
{
  Rprintf("G-spline object:\n");

  if (_dim == 0){
    Rprintf("   G-spline of dimension 0.\n");
    return;
  }

  int j;
  Rprintf("   Dimension = %d, ", _dim);
  Rprintf("   Total length = %d\n      Lengths in each dimension: ", _total_length);
  for (j = 0; j < _dim-1; j++) Rprintf("%d,  ", _length[j]);
  Rprintf("%d\n      K in each dimension: ", _length[_dim-1]);
  for (j = 0; j < _dim-1; j++) Rprintf("%d,  ", _K[j]);
  Rprintf("%d\n", _K[_dim-1]);

  char z[20];
  char z2[10];
  switch (_neighbor_system){
  case 0:  strcpy(z, "univariate CAR"); break;
  case 1:  strcpy(z, "eight neighbors"); break;
  case 2:  strcpy(z, "twelve neighbors"); break;
  default: strcpy(z, "unimplemented"); break;
  }
  strcpy(z2, (_equal_lambda ? "TRUE" : "FALSE"));
  Rprintf("   Neighboring system = %s,  order = %d,  equal lambda = %s\n", z, _order, z2);
  Rprintf("   Total index of the reference a = %d\n      Indeces per dimension: ", _total_izero);
  for (j = 0; j < _dim-1; j++) Rprintf("%d,  ", _izero[j]);
  Rprintf("%d\n", _izero[_dim-1]);
  Rprintf("   Reference log(null weight) = %g\n", _log_null_w);
  Rprintf("   c for delta = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _c4delta[j]); Rprintf("%g\n\n", _c4delta[_dim-1]);

  Rprintf("   intercept   = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _intcpt[j]); Rprintf("%g\n", _intcpt[_dim-1]);
  Rprintf("   inv. scale2 = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _invscale2[j]); Rprintf("%g\n", _invscale2[_dim-1]);
  Rprintf("   scale       = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _scale[j]); Rprintf("%g\n", _scale[_dim-1]);
  Rprintf("   gamma       = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _gamma[j]); Rprintf("%g\n", _gamma[_dim-1]);
  Rprintf("   inv. sigma2 = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _invsigma2[j]); Rprintf("%g\n", _invsigma2[_dim-1]);
  Rprintf("   sigma       = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _sigma[j]); Rprintf("%g\n", _sigma[_dim-1]);
  Rprintf("   delta       = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _delta[j]); Rprintf("%g\n", _delta[_dim-1]);
  Rprintf("   lambda      = "); for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _lambda[j]); Rprintf("%g\n", _lambda[_dim-1]);
  Rprintf("   penalty     = "); 
  if (_equal_lambda) Rprintf("%g\n", _penalty[0]);
  else{
    for (j = 0; j < _dim-1; j++) Rprintf("%g,  ", _penalty[j]); Rprintf("%g\n", _penalty[_dim-1]);
  }

  Rprintf("   a           = "); for (j = 0; j < _total_length-1; j++) Rprintf("%g,  ", _a[j]); Rprintf("%g\n", _a[_total_length-1]);
  Rprintf("   exp(a)      = "); for (j = 0; j < _total_length-1; j++) Rprintf("%g,  ", _expa[j]); Rprintf("%g\n", _expa[_total_length-1]);
  Rprintf("   sum(exp(a)) = %g\n", _sumexpa);
  Rprintf("   max(a)      = %g\n", _a_max);
  Rprintf("   Effective number of mixture components = %d\n", _k_effect);
  switch (_type_update_a){
  case Slice:         strcpy(z, "Slice sampler"); break;
  case ARS_quantile:  strcpy(z, "ARS with quantiles as starting abscissae"); break;
  case ARS_mode:      strcpy(z, "ARS with mode+- approx. sd as starting abscissae"); break;
  default:            strcpy(z, "unimplemented"); break;
  }
  Rprintf("   Type of update for 'a' = %s,  k for overrelaxation = %d\n", z, _k_overrelax_a);
  Rprintf("   Update of scale: k for overrelaxation = "); 
  for (j = 0; j < _dim-1; j++) Rprintf("%d,  ", _k_overrelax_scale[j]);  Rprintf("%d\n", _k_overrelax_scale[_dim-1]);
  Rprintf("   Update of sigma: k for overrelaxation = "); 
  for (j = 0; j < _dim-1; j++) Rprintf("%d,  ", _k_overrelax_sigma[j]);  Rprintf("%d\n", _k_overrelax_sigma[_dim-1]);

  return;
}


/***** Gspline2initArray *****/
// Write appropriate components of the G-spline to arrays parmI, parmD
// Something like an inverse of the parametric constructor
//
void
Gspline::Gspline2initArray(int* parmI,  double* parmD) const
{
  int j, i, k;
  if (_dim == 0) return;

  /** Set up integer pointers **/
  int ddim               = 0;
  int nneighbor          = ddim               + 1;
  int eequal_lambda      = nneighbor          + 1;
  int KK                 = eequal_lambda      + 1;
  int iizero             = KK                 + _dim;      /* in the G-spline: scale 0, ..., 2*K; on OUTPUT: scale -K, ..., K;  */
  int oorder             = iizero             + _dim;
  int pprior_for_lambda  = oorder             + 1;
  int pprior_for_gamma   = pprior_for_lambda  + _dim;
  int pprior_for_sigma   = pprior_for_gamma   + _dim;
  int pprior_for_intcpt  = pprior_for_sigma   + _dim;
  int pprior_for_scale   = pprior_for_intcpt  + _dim;
  int ttype_update_a     = pprior_for_scale   + _dim;
  int kk_overrelax_a     = ttype_update_a     + 1;
  int kk_overrelax_sigma = kk_overrelax_a     + 1;
  int kk_overrelax_scale = kk_overrelax_sigma + _dim;
  // int inext = kk_overrelax_scale + _dim;

  /** Set up double pointers **/
  int aa            = 0;
  int llambda       = aa            + _total_length;
  int ggamma        = llambda       + _dim;
  int ssigma        = ggamma        + _dim;
  int iintcpt       = ssigma        + _dim;
  int sscale        = iintcpt       + _dim;
  int cc4delta      = sscale        + _dim;
  int pprior_lambda = cc4delta      + _dim;
  int pprior_gamma  = pprior_lambda + 2*_dim;
  int pprior_sigma  = pprior_gamma  + 2*_dim;
  int pprior_intcpt = pprior_sigma  + 2*_dim;
  int pprior_scale  = pprior_intcpt + 2*_dim;
  // int dnext = pprior_scale + 2*_dim;

  parmI[ddim] = _dim;
  parmI[nneighbor] = _neighbor_system;
  parmI[eequal_lambda] = (_equal_lambda ? 1 : 0);
  for (j = 0; j < _dim; j++){
    parmI[KK + j] = _K[j];
    parmI[iizero + j] = _izero[j] - _K[j];                     /* give R index on the scale -K,...,K */
    parmI[pprior_for_lambda + j]  = _prior_for_lambda[j];
    parmI[pprior_for_gamma + j]   = _prior_for_gamma[j];
    parmI[pprior_for_sigma + j]   = _prior_for_sigma[j];
    parmI[pprior_for_intcpt + j]  = _prior_for_intcpt[j];
    parmI[pprior_for_scale + j]   = _prior_for_scale[j];
    parmI[kk_overrelax_sigma + j] = _k_overrelax_sigma[j];
    parmI[kk_overrelax_scale + j] = _k_overrelax_scale[j];

    parmD[llambda + j]  = _lambda[j];
    parmD[ggamma + j]   = _gamma[j];
    parmD[ssigma + j]   = _sigma[j];
    parmD[iintcpt + j]  = _intcpt[j];
    parmD[sscale + j]   = _scale[j];
    parmD[cc4delta + j] = _c4delta[j];

    switch (_prior_for_gamma[j]){
    case _Fixed_:
      for (i = 0; i < 2; i++) parmD[pprior_gamma + 2*j] = _prior_gamma[2*j+i];  /* = 0.0   */
      break;
    case Normal:
      parmD[pprior_gamma + 2*j] = _prior_gamma[2*j];
      parmD[pprior_gamma + 2*j + 1] = 1/_prior_gamma[2*j + 1];    /* give prior variance, whereas prior inverse variance was stored in _prior_gamma  */
      break;
    default:
      throw returnR("C++ Error: Unimplemented prior for gamma of the G-spline appeared in Gspline::Gspline2initArray()", 1);
    }

    switch (_prior_for_intcpt[j]){
    case _Fixed_:
      for (i = 0; i < 2; i++) parmD[pprior_intcpt + 2*j] = _prior_intcpt[2*j+i];  /* = 0.0   */
      break;
    case Normal:
      parmD[pprior_intcpt + 2*j] = _prior_intcpt[2*j];
      parmD[pprior_intcpt + 2*j + 1] = 1/_prior_intcpt[2*j + 1];    /* give prior variance, whereas prior inverse variance was stored in _prior_intcpt  */
      break;
    default:
      throw returnR("C++ Error: Unimplemented prior for intcpt of the G-spline appeared in Gspline::Gspline2initArray()", 1);
    }

    switch (_prior_for_sigma[j]){
    case Fixed_:
    case Gamma:  
      for (i = 0; i < 2; i++) parmD[pprior_sigma + 2*j+i] = _prior_sigma[2*j+i];
      break;
    case SDUnif: 
      parmD[pprior_sigma + 2*j] = _prior_sigma[2*j];                   /* = 0.0                                             */
      parmD[pprior_sigma + 2*j + 1] = sqrt(1/_prior_sigma[2*j+1]);     /* give S, whereas 1/S^2 was stored in _prior_sigma  */
      break;
    default:     
      throw returnR("C++ Error: Unimplemented prior for sigma of the G-spline appeared in Gspline::Gspline2initArray()", 1);
    }

    switch (_prior_for_scale[j]){
    case Fixed_:
    case Gamma:  
      for (i = 0; i < 2; i++) parmD[pprior_scale + 2*j+i] = _prior_scale[2*j+i];
      break;
    case SDUnif: 
      parmD[pprior_scale + 2*j] = _prior_scale[2*j];                   /* = 0.0                                             */
      parmD[pprior_scale + 2*j + 1] = sqrt(1/_prior_scale[2*j+1]);     /* give S, whereas 1/S^2 was stored in _prior_scale  */
      break;
    default:     
      throw returnR("C++ Error: Unimplemented prior for scale of the G-spline appeared in Gspline::Gspline2initArray()", 1);
    }

  }
  parmI[oorder]         = _order;
  parmI[ttype_update_a] = _type_update_a;
  parmI[kk_overrelax_a] = _k_overrelax_a;

  for (i = 0; i < _total_length; i++){
    parmD[aa + i] = _a[i]; 
  }

  for (k = 0; k <= (_dim-1)*(!_equal_lambda); k++){
    switch (_prior_for_lambda[k]){
    case Fixed_:
    case Gamma:  
      for (i = 0; i < 2; i++) parmD[pprior_lambda + 2*k+i] = _prior_lambda[2*k+i];
      break;
    case SDUnif: 
      parmD[pprior_lambda + 2*k] = _prior_lambda[2*k];                   /* = 0.0                                             */
      parmD[pprior_lambda + 2*k + 1] = sqrt(1/_prior_lambda[2*k+1]);     /* give S, whereas 1/S^2 was stored in _prior_lambda */
      break;
    default:     
      throw returnR("C++ Error: Unimplemented prior for lambda of the G-spline appeared in Gspline::Gspline2initArray()", 1);
    }
  }
  if (_equal_lambda){
    for (k = 1; k < _dim; k++){
      for (i = 0; i < 2; i++) parmD[pprior_lambda + 2*k+i] = _prior_lambda[2*k+i];     /* = 0.0 */
    }
  }

  return;
}

/***** a2expa *****/
// Compute an array _expa and variables _sumexpa, _sumexpa_margin,_k_effect, _ind_w_effect from _a
//
void
Gspline::a2expa()
{
  int ia, j, k0, k1;

  _sumexpa = 0.0;
  _k_effect = 0;

  bool a_inf = false;
  switch (_dim){
  case 1:
    for (ia = 0; ia < _total_length; ia++){
      if (_a[ia] - _a_max > _log_null_w){
        _ind_w_effect[_k_effect] = ia;
        _k_effect++;
      }

      if (_a[ia] >= _log_inf){
        _expa[ia] = FLT_MAX;
        a_inf = true;
      }
      else{
        _expa[ia] = exp(_a[ia]);
        _sumexpa += _expa[ia];
      }
    }
    break;

  case 2:
    for (j = 0; j < _dim; j++){
      for (k1 = 0; k1 < _length[j]; k1++) _sumexpa_margin[j][k1] = 0.0;
    }
    for (ia = 0; ia < _total_length; ia++){
      k0 = ia % _length[0];
      k1 = ia / _length[0];
      if (_a[ia] - _a_max > _log_null_w){
        _ind_w_effect[_k_effect] = ia;
        _k_effect++;
      }
       
      if (_a[ia] >= _log_inf){
        _expa[ia] = FLT_MAX;
        a_inf = true;
        _sumexpa_margin[0][k0] = FLT_MAX;
        _sumexpa_margin[1][k1] = FLT_MAX;
      }
      else{
        _expa[ia] = exp(_a[ia]);
        _sumexpa += _expa[ia];
        _sumexpa_margin[0][k0] += _expa[ia];
        _sumexpa_margin[1][k1] += _expa[ia];
      }
    }
    break;

  default:
    throw returnR("C++ Error: Function Gspline::a2expa() not yet implemented for _dim > 2", 1);
  }

  if (a_inf){
    _sumexpa = FLT_MAX;
  }
  return;
}  /** end of function a2expa **/


/***** change_a *****/
// Update an array _expa and a variable _sumexpa if one 'a' changes
// Do not update _k_effect and _ind_w_effect 
//
// newa ... new value of 'a' coefficient
// ia ..... index of a new 'a'
//
void
Gspline::change_a(const double* newa, const int& ia)
{
  int k0, k1;
  if (ia < 0 || ia >= _total_length) throw returnR("C++ Error: Incorrect ia in Gspline:change_a", 1);

  _a[ia] = *newa;

  _sumexpa -= _expa[ia];
  switch (_dim){
  case 1:
    if (_a[ia] >= _log_inf){
       _expa[ia] = FLT_MAX;
       _sumexpa = FLT_MAX;
    }
    else{
      _expa[ia] = exp(_a[ia]);
      _sumexpa += _expa[ia];
    }
    break;

  case 2:
    k0 = ia % _length[0];
    k1 = ia / _length[0];
    _sumexpa_margin[0][k0] -= _expa[ia];
    _sumexpa_margin[1][k1] -= _expa[ia];

    if (_a[ia] >= _log_inf){
       _expa[ia] = FLT_MAX;
       _sumexpa = FLT_MAX;
       _sumexpa_margin[0][k0] = FLT_MAX;
       _sumexpa_margin[1][k1] = FLT_MAX;
    }
    else{
      _expa[ia] = exp(_a[ia]);
      _sumexpa += _expa[ia];
      _sumexpa_margin[0][k0] += _expa[ia];
      _sumexpa_margin[1][k1] += _expa[ia];
    }
    break;

  default:
    throw returnR("C++ Error: Function Gspline::change_a() not yet implemented for _dim > 2", 1);
  }

  return;
}  /** end of function change_a **/


/***** update_k_effect *****/
// Update _k_effect and an array _ind_w_effect 
//
void
Gspline::update_k_effect()
{
  int ia;
  _k_effect = 0;

  for (ia = 0; ia < _total_length; ia++){
    if (_a[ia] - _a_max > _log_null_w){
      _ind_w_effect[_k_effect] = ia;
      _k_effect++;
    }
  }
  return;
}


/***** update_a_max *****/
// * update _a_max and adjust a's if _a_max is out of its limits (do it only if (_center_a == false))
// * update also _k_effect
// * if (_center_a == true) perform also centering
// * if needed, update also _expa, _sumexpa etc.
//
void
Gspline::update_a_max_center_and_k_effect()
{
  int i;
  double tmp;

  _a_max = _a[0];
  for (i = 1; i < _total_length; i++) if (_a[i] > _a_max) _a_max = _a[i];

  if (_center_a){
    tmp = 0.0;
    for (i = 0; i < _total_length; i++) tmp += _a[i];
    tmp /= _total_length;
    for (i = 0; i < _total_length; i++) _a[i] -= tmp;
    _a_max -= tmp;
    a2expa();    /** this updates also _k_effect **/
  }
  else{
    if (_a_max > _a_ceil){
      tmp = _a_max - _a_ceil;
      for (i = 0; i < _total_length; i++) _a[i] -= tmp;
      _a_max = _a_ceil;
      a2expa();    /** this updates also _k_effect **/
    }
    else{
      if (_a_max < _amax_floor){
        tmp = _amax_floor -_a_max;
        for (i = 0; i < _total_length; i++) _a[i] += tmp;
        _a_max = _amax_floor;
        a2expa();    /** this updates also _k_effect **/
      }
      else{
        update_k_effect();
      }
    }
  }

  return;
}


/***** moments *****/
// Compute means and a covariance (matrix) of a univariate or bivariate G-spline
//
// EY ........ mean for each dimension                            [_dim]
// varY ...... covariance matrix (varY1, cov(Y1, Y2), var(Y2))    [0.5*_dim*(_dim+1)]
//
void
Gspline::moments(double* EY, double* varY) const
{
  int i, j, k0, k1;
  double mu_jk, mu_jk1;

  switch (_dim){
  case 1:
    /* expectation of epsilon */
    EY[0] = 0.0;
    for (i = 0; i < _k_effect; i++) EY[0] += _expa[_ind_w_effect[i]] * mu(0, _ind_w_effect[i]);
    EY[0] /=_sumexpa;
    
    /* variance */
    varY[0] = 0.0;
    for (i = 0; i < _k_effect; i++){
      mu_jk = mu(0, _ind_w_effect[i]) - EY[0];
      varY[0] += _expa[_ind_w_effect[i]] * mu_jk * mu_jk;
    }
    varY[0] /= _sumexpa;
    varY[0] += _sigma[0] * _sigma[0];
    varY[0] *= _scale[0] * _scale[0];

    /* expectation of Y */
    EY[0] *= _scale[0];
    EY[0] += _intcpt[0];
    break;

  case 2:
    for (j = 0; j < _dim; j++){
      /* expectation of epsilon */
      EY[j] = 0.0;
      mu_jk = _gamma[j] - _K[j]*_delta[j];
      for (k0 = 0; k0 < _length[j]; k0++){
        EY[j] += _sumexpa_margin[j][k0] * mu_jk;
        mu_jk += _delta[j];
      }
      EY[j] /= _sumexpa;

      /* variance */
      varY[2*j] = 0.0;
      mu_jk = _gamma[j] - _K[j]*_delta[j];
      for (k0 = 0; k0 < _length[j]; k0++){
        varY[2*j] += _sumexpa_margin[j][k0] * (mu_jk - EY[j]) * (mu_jk - EY[j]);
        mu_jk += _delta[j];
      }
      varY[2*j] /= _sumexpa;
      varY[2*j] += _sigma[j] * _sigma[j];
      varY[2*j] *= _scale[j] * _scale[j];
    }

    /* covariance */
    varY[1] = 0.0;
    mu_jk = _gamma[0] - _K[0]*_delta[0];    
    for (k0 = 0; k0 < _length[0]; k0++){
      mu_jk1 = _gamma[1] - _K[1]*_delta[1];
      for (k1 = 0; k1 < _length[1]; k1++){
        varY[1] += _expa[k1*_length[0] + k0] * (mu_jk - EY[0]) * (mu_jk1 - EY[1]);
        mu_jk1 += _delta[1];
      }
      mu_jk += _delta[0];
    }
    varY[1] /= _sumexpa;
    varY[1] *= _scale[0] * _scale[1];

    /* expectation of Y */
    for (j = 0; j < _dim; j++){
      EY[j] *= _scale[j];
      EY[j] += _intcpt[j];
    }

    break;

  default:
    throw returnR("C++ Error: Function Gspline::moments not yet implemented for _dim > 2", 1);
  }

  return;
}


/***** mean_univariate *****/
// Compute a mean of a univariate G-spline
// * this function is primarily used by ''adjust_intcpt' function inside 'bayessurvreg2'
//   in the case of G-spline random intercept
//
// EY ........ computed mean
//
void
Gspline::mean_univariate(double* EY) const
{
  if (_dim != 1) throw returnR("C++ Error: Function Gspline::mean_univariate is only implemented for _dim = 1", 1);

  int i;

  *EY = 0.0;
  for (i = 0; i < _k_effect; i++) (*EY) += _expa[_ind_w_effect[i]] * mu(0, _ind_w_effect[i]);
  (*EY) /=_sumexpa;
  (*EY) *= _scale[0];
  (*EY) += _intcpt[0];

  return;
}


/***** adjust_intcpt *****/
//
// Add value(s) of adj to _intcpt
//
// * used in the model with G-spline random intercept
// * called by the function 'adjust_intcpt' defined in 'bayessurvreg2.cpp'
//
void
Gspline::adjust_intcpt(double* adj)
{
  int j;
  for (j = 0; j < _dim; j++) _intcpt[j] += adj[j];
  return;
}


