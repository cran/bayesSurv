/*** rhoNorm.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  04/12/2006
//
//    PURPOSE:  Routines to sample correlation coefficient (or its Fisher Z transform)
//              with the bivariate normal likelihood
//
//    DOCUMENTATION: rhoNorm.pdf in ~/Rlib/akGALMM/doc
//
//    z2rho:         04/12/2006
//    rho2z:         04/12/2006
//    rho2zError:    06/12/2006
//    update_pUnif:  04/12/2006
//    ML_est:        04/12/2006
//    lposter0:      04/12/2006
//    lposter1:      05/12/2006
//    lposter2:      04/12/2006
//
//    mcmc_rhoNorm:  05/12/2006
//
#include "rhoNorm.h"

namespace rhoNorm {

/*** z2rho                              ***/
/* rho = 2/(1 + exp(-2*z)) - 1 = tanh(z)  */
/*                                        */
void
z2rho(double *rho,  const double *z)
{
  if (*z < -rhoNorm::zeps) *rho = -1;
  else{
    if (*z > rhoNorm::zeps) *rho = 1;
    else                    *rho = 2/(1 + exp(-2*(*z))) - 1;
  }
  
  return;
}


/*** rho2z                                     ***/
/* z = -0.5*log((1-rho)/(1+rho)) = arctanh(rho)  */
/*                                               */
void
rho2z(double *z,  const double *rho)
{
  if (*rho > rhoNorm::rhoONE) *z = rhoNorm::zeps;
  else{
    if (*rho < rhoNorm::rhoMinONE) *z = -rhoNorm::zeps;
    else                           *z = -0.5*log((1 - (*rho))/(1 + (*rho)));      
  }

  return;
}


/*** rho2zError                                ***/
/* z = -0.5*log((1-rho)/(1+rho)) = arctanh(rho)  */
/*                                               */
void
rho2zError(double *z,  const double *rho)
{
  if (*rho > rhoNorm::rhoONE) throw returnR("Error in rhoNorm.cpp: rho2zError. rho is too close to 1", 1);
  else{
    if (*rho < rhoNorm::rhoMinONE) throw returnR("Error in rhoNorm.cpp: rho2zError. rho is too close to -1", 1);
    else                           *z = -0.5*log((1 - (*rho))/(1 + (*rho)));      
  }

  return;
}


/***** update_pUnif:  Update the correlation coefficient (1 iteration of MCMC)                   *****/
/*****     * uniform prior on rho is assumed                                                     *****/
/*****     * sampling is performed on the scale of the Fisher transformation of rho              *****/
/*****     * sampling strategy:                                                                  *****/
/*****       a) normal approximation to the full conditional which is used as a proposal         *****/
/*****          in the Metropolis-Hastings algorithm                                             *****/
/*****       b) Langevin algorithm (Robert and Casella, 2004, Monte Carlo Statistical Methods,   *****/
/*****          Second Edition, p. 318-320)                                                      *****/
/*****       c) if above fails, ARMS (Gilks, Best, Tan, 1995, Applied Statistics 44, 455-472)    *****/
/*****          is tried                                                                         *****/
/*****          NOT IMPLEMENTED                                                                  *****/
//
// z[1]:     INPUT:  Current value of z
//           OUTPUT: New value of z
// rho[1]:   INPUT:  Arbitrary
//           OUTPUT: New value of rho
// work[3]:  INPUT:  Arbitrary
//           OUTPUT: Values of posterior log-density, its first and minus second derivatives at current z
//
// algorithm[1]:  Type of the algorithm used to update rho
//                See rhoNorm.h for possible values
// scaleL[1]:     Scaling factor sigma for the Langevin algorithm
//
//
void
update_pUnif(int *accept,           double *z,             double *rho,          double *work,
             const double *sumu2,   const double *sumv2,   const double *sumuv,  const int *nobs,
             const int *algorithm,  const double *scaleL)
{
  static double *ll, *dll, *ddll;
  static double scaleL2, u, log_scaleL2, prop_log_q, log_q, log_A;
  static double prop_log_poster, prop_dlog_poster, prop_ddlog_poster;
  static double temp_log_poster, temp_dlog_poster, temp_ddlog_poster;
  static double prop_mean[1], prop_rho[1], temp_rho[1], prop_z[1];
  static int err[1], niter[1];
  static bool do_arms;

  ll   = work;
  dll  = ll  + 1;
  ddll = dll + 1;

  switch (*algorithm){
  case _Normal_Around_Mode_:

    /** Construct normal approximation **/
    rhoNorm::lposter2(ll, dll, ddll, rho, z, sumu2, sumv2, sumuv, nobs);
    *prop_mean        = *z;  
    *prop_rho         = *rho;
    prop_log_poster   = *ll;
    prop_dlog_poster  = *dll;
    prop_ddlog_poster = *ddll;  
    rhoNorm::ML_est(&prop_log_poster, &prop_dlog_poster, &prop_ddlog_poster, 
                    prop_mean, prop_rho, niter, err, sumu2, sumv2, sumuv, nobs, &rhoNorm::_maxiter);
    //Rprintf("\nniter=%d,  z(%g)=%g,  mode=%g,  ll=%g,  dll=%g,  min ddll=%g", 
    //        *niter, *ll, *z, *prop_mean, prop_log_poster, prop_dlog_poster, prop_ddlog_poster);
    if (*err >= 2 || prop_ddlog_poster <= 0){
      throw returnR("Trap in rhoNorm.cpp: update_pUnif. Not possible to construct normal approximation. Consider usage of the Langevin algorithm.", 1);
    }
    else{
      /** Propose new value of z by sampling from the normal approximation **/
      *prop_z    = norm_rand();
      prop_log_q = 0.5*(log_AK(prop_ddlog_poster) - (*prop_z)*(*prop_z));
      *prop_z   /= sqrt(prop_ddlog_poster);
      *prop_z   += *prop_mean;

      /** Construct the reversal normal approximation **/
      rhoNorm::lposter2(&prop_log_poster, &prop_dlog_poster, &prop_ddlog_poster, prop_rho, prop_z, sumu2, sumv2, sumuv, nobs);
      *prop_mean        = *prop_z;
      *temp_rho         = *prop_rho;
      temp_log_poster   = prop_log_poster;
      temp_dlog_poster  = prop_dlog_poster;
      temp_ddlog_poster = prop_ddlog_poster;
      rhoNorm::ML_est(&temp_log_poster, &temp_dlog_poster, &temp_ddlog_poster,
                      prop_mean, temp_rho, niter, err, sumu2, sumv2, sumuv, nobs, &rhoNorm::_maxiter);
      //Rprintf("\n   REVERSAL: niter=%d,  z=%g,  mode=%g,  ll=%g,  dll=%g,  min ddll=%g", 
      //        *niter, *prop_z, *prop_mean, prop_log_poster, prop_dlog_poster, prop_ddlog_poster);
      if (*err >= 2 || prop_ddlog_poster <= 0){
        log_A   = _AK_EMIN - 1;
        *accept = 0;
        return;
      }
      else{
        u     = (*z - (*prop_mean))*sqrt(prop_ddlog_poster);
        log_q = 0.5*(log_AK(prop_ddlog_poster) - u*u);
      }
    }
    break;

  case _Langevin_:
    scaleL2     = (*scaleL)*(*scaleL);
    log_scaleL2 = log_AK(scaleL2);

    /** Construct Langevin normal approximation **/
    rhoNorm::lposter1(ll, dll, rho, z, sumu2, sumv2, sumuv, nobs);
    if (!R_finite(*ll)){
      throw returnR("Trap in rhoNorm.cpp: update_pUnif. Value of the correlation is too close to +-1.", 1);
    }
    else{
      *prop_mean = *z + scaleL2*(*dll); 

      /** Propose new value of z by sampling from the normal approximation **/
      *prop_z    = norm_rand();
      prop_log_q = 0.5*(-log_scaleL2 - (*prop_z)*(*prop_z));
      *prop_z   *= *scaleL;
      *prop_z   += *prop_mean;

      /** Construct the reversal Langevin normal approximation **/
      rhoNorm::lposter1(&prop_log_poster, &prop_dlog_poster, prop_rho, prop_z, sumu2, sumv2, sumuv, nobs);
      if (!R_finite(prop_log_poster)){
        log_A   = _AK_EMIN - 1;
        *accept = 0;
        return;
      }
      else{
        *prop_mean = *prop_z + scaleL2*prop_dlog_poster;         
        u          = (*z - (*prop_mean))/(*scaleL);
        log_q      = 0.5*(-log_scaleL2 - u*u);
      }
    }
    break;

  default:
    throw returnR("Error in rhoNorm.cpp: update_pUnif. Unknown algorithm required.", 1);
  }
  
  /** Logarithm of the acceptance ratio and acceptance test **/
  log_A = prop_log_poster + log_q - (*ll) - prop_log_q;
  //Rprintf("\n   log_A=%g,  prop_log_poster=%g,  log_q=%g,  ll=%g,  prop_log_q=%g",  log_A, prop_log_poster, log_q, *ll, prop_log_q);
  if (log_A < _AK_EMIN){
    *accept = 0;
    return;
  }
  if (log_A >= 0){
    *accept = 1;
  }
  else{  /* decide by sampling from exponential distribution */
    *ll = exp_rand();    
    *accept = (*ll > -log_A ? 1 : 0);
  }
    
  if (*accept){
    *z    = *prop_z;
    *rho  = *prop_rho;
    *ll   = prop_log_poster;
    *dll  = prop_dlog_poster;
    *ddll = prop_ddlog_poster;
  }
  return;
}


/***** ML_est:  Find the mode of the posterior log-density                                 *****/
//
// INPUT:
//    ll[1]:    Current value of the posterior log-density
//    dll[1]:   Current value of the first derivative of the posterior log-density w.r.t. z
//    ddll[1]:  Current value of MINUS the second derivative of the posterior log-density w.r.t. z
//    z[1]:     Current value of z
//    rho[1]:   Current value of rho
void
ML_est(double *ll,             double *dll,           double *ddll,
       double *z,              double *rho,           int *niter,           int *err,
       const double *sumu2,    const double *sumv2,   const double *sumuv,  const int *nobs,
       const int *maxiter)
{
  static double NR_step, old_z, old_ll, relat_diff;
  static int halfstep;

  *err = 0;

  /** Check initials for NaN **/
  if (!R_finite(*ll)){
    //REprintf("Trap in rhoNorm.cpp: ML_est. Initial values lead to the log-likelihood of -Inf.\n");
    *err = 4;
    return;
  }

  /*** Iterate (Newton-Raphson) ***/
  for (*niter = 0; *niter < *maxiter; (*niter)++){

    if (*ddll <= 0){
      //REprintf("Trap in rhoNorm.cpp: ML_est. Non-positive minus second derivative.\n");
      *err = 3;
      return;
    }

    /*** Newton-Raphson step ***/
    NR_step = (*dll)/(*ddll);

    /*** Compute a new value of z and rho ***/
    old_z = *z;
    *z += NR_step;

    /*** Update derivatives and rho ***/
    old_ll = *ll;
    rhoNorm::lposter2(ll, dll, ddll, rho, z, sumu2, sumv2, sumuv, nobs);

    /*** Check convergence ***/
    relat_diff = R_finite(*ll) ? fabs(1 - old_ll/(*ll)) : R_PosInf;
    if (relat_diff <= rhoNorm::_toler){
      break;
    }

    /*** If not yet convergence, check whether the objective function increases ***/
    /***   if no increase perform step-halving                                  ***/
    if (!R_finite(*ll) || *ll < old_ll){      
      for (halfstep = 0; halfstep < rhoNorm::_max_stephalf; halfstep++){
        NR_step *= 0.5;
        *z -= NR_step;
	rhoNorm::lposter0(ll, rho, z, sumu2, sumv2, sumuv, nobs);

        if (*ll >= old_ll){
          rhoNorm::lposter2(ll, dll, ddll, rho, z, sumu2, sumv2, sumuv, nobs);
          break;
        }
      }
      if (halfstep == rhoNorm::_max_stephalf){
        *z = old_z;
        rhoNorm::lposter2(ll, dll, ddll, rho, z, sumu2, sumv2, sumuv, nobs);
        *err = 2;
        break;
      }
    }
  }

  /*** Check niter ***/
  if (*maxiter && *niter == *maxiter) *err = 1;
  else                                (*niter)++;     /* to get the number of iterations really performed */

  return;
}


/***** lposter0, lposter1, lposter2                                                                 ***/
/*** Posterior log-density of z and its derivatives                                                 ***/
/*** when uniform prior on rho is combined with the bivariate normal likelihood                     ***/ 
/***      (y[1],y[2]) ~ N((mu[1], mu[2]), Sigma), Sigma[1,1] = sigma[1]^2                           ***/
/***                                              Sigma[2,2] = sigma[2]^2                           ***/
/***                                              Sigma[1,2] = Sigma[2,1] = rho*sigma[1]*sigma[2]   ***/
/***                                                                                                ***/
//
// ll[1]:     Computed value of log-density
// dll[1]:    Computed value of the first derivative of the log-density
// ddll[1]:   Computed value of MINUS the second derivative of the log-density
// rho[1]:    Computed value of rho corresponding to z
//
// z[1]:      Input z
// sumu2[1]:  sum(u[i]^2), where u[i] = (y[i,1] - mu[1])/sigma[1]
// sumv2[1]:  sum(v[i]^2), where v[i] = (y[i,2] - mu[2])/sigma[2]
// sumuv[1]:  sum(u[i]*v[i])
// nobs:      number of observations entering the log-likelihood
//
void
lposter0(double *ll,           double *rho,          const double *z,      
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs)
{
  static double one_rho2;

  if (*z < -rhoNorm::zeps){
    *rho = -1;
    *ll = R_NegInf;
    return;
  }
  if (*z > rhoNorm::zeps){
    *rho = 1;
    *ll = R_NegInf;
    return;
  }

  *rho     = 2/(1 + exp(-2*(*z))) - 1;
  one_rho2 = 1 - (*rho)*(*rho);
  *ll      = (1 - (*nobs)/2)*log(one_rho2) - (1/(2*one_rho2))*(*sumu2 + (*sumv2)) + (*rho/one_rho2)*(*sumuv);

  return;    
}


void
lposter1(double *ll,           double *dll,          double *rho,          const double *z,  
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs)
{
  static double rho2, one_rho2, dd_rho1, dd_rho2;

  if (*z < -rhoNorm::zeps){
    *rho = -1;
    *ll   = R_NegInf;
    *dll  = R_NegInf;
    return;
  }
  if (*z > rhoNorm::zeps){
    *rho = 1;
    *ll   = R_NegInf;
    *dll  = R_NegInf;
    return;
  }

  *rho     = 2/(1 + exp(-2*(*z))) - 1;
  rho2     = (*rho)*(*rho);
  one_rho2 = 1 - rho2;
  dd_rho1  = *rho/one_rho2;
  dd_rho2  = (1 + rho2)/one_rho2;

  *ll   = (1 - (*nobs)/2)*log(one_rho2) - (1/(2*one_rho2))*(*sumu2 + (*sumv2)) + dd_rho1*(*sumuv);
  *dll  = (*nobs - 2)*(*rho) - dd_rho1*(*sumu2 + (*sumv2)) + dd_rho2*(*sumuv);

  return;    
}


void
lposter2(double *ll,           double *dll,          double *ddll,         double *rho,          const double *z,  
         const double *sumu2,  const double *sumv2,  const double *sumuv,  const int *nobs)
{
  static double rho2, one_rho2, dd_rho1, dd_rho2;

  if (*z < -rhoNorm::zeps){
    *rho = -1;
    *ll   = R_NegInf;
    *dll  = R_NegInf;
    *ddll = R_NegInf;
    return;
  }
  if (*z > rhoNorm::zeps){
    *rho = 1;
    *ll   = R_NegInf;
    *dll  = R_NegInf;
    *ddll = R_NegInf;
    return;
  }

  *rho     = 2/(1 + exp(-2*(*z))) - 1;
  rho2     = (*rho)*(*rho);
  one_rho2 = 1 - rho2;
  dd_rho1  = *rho/one_rho2;
  dd_rho2  = (1 + rho2)/one_rho2;

  *ll   = (1 - (*nobs)/2)*log(one_rho2) - (1/(2*one_rho2))*(*sumu2 + (*sumv2)) + dd_rho1*(*sumuv);
  *dll  = (*nobs - 2)*(*rho) - dd_rho1*(*sumu2 + (*sumv2)) + dd_rho2*(*sumuv);
  *ddll = -(*nobs - 2)*one_rho2 + dd_rho2*(*sumu2 + (*sumv2)) - 4*dd_rho1*(*sumuv);

  return;    
}


//
// algorithm[1]:  0 = _Normal_Around_Mode_
//                1 = _Langevin_
// scaleL[1]:     Scale factor for the Langevin algorithm
//
extern "C"{
  void
  mcmc_rhoNorm(int *acceptSample,     double *zSample,      double *rhoSample,     int *iter,
               const double *sumu2,   const double *sumv2,  const double *sumuv,   const int *nobs,  const int *nsimul,
               const int *algorithm,  const double *scaleL)
  {
    try{
      int i;

      GetRNGstate(); 

      /*** Numbers of iterations etc. ***/
      const int niter = nsimul[0];
      const int nthin = nsimul[1];
      const int nwrite = nsimul[2];
 

      /*** Current values of parameters, related quantities and their initialization ***/
      double z[1]   = {*zSample};
      double rho[1]; 
      double ll_dll_ddll[3];
      lposter2(ll_dll_ddll, ll_dll_ddll+1, ll_dll_ddll+2, rho, z, sumu2, sumv2, sumuv, nobs);
      *rhoSample = *rho;


      /*** MCMC ***/
      int accept[1];

      int *acceptSampleP    = acceptSample;
      double *zSampleP      = zSample;
      double *rhoSampleP    = rhoSample;

      int nullthIter   = *iter;
      int lastIter     = nullthIter + niter;
      int iterTotal    = nullthIter*nthin;
      int iterTotalNow = 0;
      int backs        = 0;
      int writeAll     = 0;
      int witer;
  
      Rprintf("Iteration ");    
      for (*iter=nullthIter + 1; *iter <= lastIter; (*iter)++){
        for (witer = 1; witer <= nthin; witer++){               /**  thinning cycle  **/
          iterTotal++;                                          /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                                       /* = (iter-1)*nthin + witer - nullthIter*nthin */

	  rhoNorm::update_pUnif(accept, z, rho, ll_dll_ddll, sumu2, sumv2, sumuv, nobs, algorithm, scaleL);
          //goto CLEANING;
        }

        *acceptSampleP = *accept;
        acceptSampleP++;

        *zSampleP = *z;
        zSampleP++;

        *rhoSampleP = *rho;
        rhoSampleP++;

        if (!(*iter % nwrite) || *iter == lastIter){
          writeAll = 1;
          for (i = 0; i < backs; i++) Rprintf("\b");
            Rprintf("%d", *iter);
            backs = int(log10(double(*iter))) + 1;
        }
      }

      /*** Cleaning ***/
    CLEANING:
      Rprintf("\n");

      PutRNGstate();
      return;
    }
    catch (returnR rr){
      PutRNGstate();
      return;
    }
  }
}

} /*** end of the namespace rhoNorm ***/
