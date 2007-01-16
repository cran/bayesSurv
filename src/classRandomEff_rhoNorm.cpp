/*** classRandomEff_rhoNorm.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//    CREATED:  06/12/2006
//
//    PURPOSE:  Routines for version = 31 of bayessurvreg2
//              * doubly-interval-censored data
//              * univariate random effect d for both onset and univariate random effect b for time-to-event
//              * distribution of (d, b) is modelled as a Gspline with correlated basis
//            
#include "classRandomEff_rhoNorm.h"

//
// ***** Gspl_rho_intcpt_update *****
//
// * This function assumes that: d->_nRandom == 1 = b->_nRandom == 1
//                               gg_d->_dim = gg_b->_dim = gg_zeta->_dim = gg_eps->_dim = 1
//
// * Update random intercepts d and b in the model with doubly-interval-censored data
//      under the assumption that (d, b) is distributed as a G-spline with correlated basis
//
// * Update also the value of the correlation coefficient
//
// * Update also regression residuals
//
//
// d:                   random effect object for onset
// b:                   random effect object for time-to-event
// rho_zb[2]:           rho_zb[0] = rho
//                      rho_zb[1] = z = atanh(rho)
//
// regResOnset[nP]:     regression residuals for onset
// regResTime[nP]:      regression residuals for time-to-event
//
// rho_accept[1]:       acceptance indicator for the update of rho
//
// nP[1]:               total number of observations (over all clusters)
// rho_algor[1]:        type of the algorithm used to update rho (see rhoNorm.h for possible values)
// rho_scaleL[1]:       value of the scale parameter for the Langevin update of rho             
//
// gg_d:                                G-spline giving the distribution of the random intercept in the onset part
// mu_d[gg_d->dim(), gg_d->length(j)]:  already computed knots of the onset random intercept G-spline
// rM_d[d._nCluster]:                   allocation labels for the onset random intercept
//                                      \in {0, ..., gg_d->total_length() - 1}
//
// gg_b:                                G-spline giving the distribution of the random intercept in the time-to-event part
// mu_b[gg_b->dim(), gg_b->length(j)]:  already computed knots of the time-to-event random intercept G-spline
// rM_b[b._nCluster]:                   allocation labels for the time-to-event random intercept
//                                      \in {0, ..., gg_b->total_length() - 1}
//
// gg_zeta:                                      G-spline giving the distribution of the error in the onset part
// mu_zeta[gg_zeta->dim(), gg_zeta->length(j)]:  already computed knots of the onset error G-spline
// rM_zeta[nP]:                                  allocation labels for the onset error
//                                               \in {0, ..., gg_zeta->total_length() - 1}
//
// gg_eps:                                    G-spline giving the distribution of the error in the time-to-event part
// mu_eps[gg_eps->dim(), gg_eps->length(j)]:  already computed knots of the time-to-event error G-spline
// rM_eps[nP]:                                allocation labels for the time-to-event error
//                                            \in {0, ..., gg_eps->total_length() - 1}
//
void
Gspl_rho_intcpt_update(RandomEff *d,            RandomEff *b,              double *rho_zb,
                       double *regResOnset,     double *regResTime,        int *rho_accept,
                       const int *nP,           const int *rho_algor,      const double *rho_scaleL,
                       const Gspline *gg_d,     double** const mu_d,       const int *rM_d,
                       const Gspline *gg_b,     double** const mu_b,       const int *rM_b,
                       const Gspline *gg_zeta,  double** const mu_zeta,    const int *rM_zeta,
                       const Gspline *gg_eps,   double** const mu_eps,     const int *rM_eps)
{
  if (!d->_nRandom || !b->_nRandom) return;
  if (gg_d->dim() != 1 || gg_b->dim() != 1 || gg_zeta->dim() != 1 || gg_eps->dim() != 1)
    throw returnR("Error in classRandomEff_rhoNorm.cpp: Gspl_rho_intcpt_update. Not implemented for multivariate G-splines", 1);
  if (d->_nCluster != b->_nCluster)
    throw returnR("Error in classRandomEff_rhoNorm.cpp: Gspl_rho_intcpt_update. d->_nCluster != b->_nCluster", 1);

  static const int _TWO_[1] = {2};
  static int info[1];
  static double ll_dll_ddll[3];

  static int cl, i;
  static double one_rho2, inv_one_rho2;
  static double invsigscale_d, invsigscale_b, invsigscale2_d, invsigscale2_b, invsigscale2_zeta, invsigscale2_eps;
  static double mean_d, mean_b, u, v;
  static double sumu2, sumv2, sumuv;
  static double *tempP, *temp2P;
  static double *regResOnsetP, *regResTimeP;
  static double *dP, *bP;
  static double ivar_prop[3], ivar_propCommon[3], mean_prop[2], mean_propTemp[2];
        /* used to store lower triangle of the inverse variance of the full conditional distribution of (d, b) */
        /* and the canonical mean of the full conditional distribution of (d, b)                               */

  static const int *rdP, *rbP, *rzetaP, *repsP, *nwithinClP;  

  /***** UPDATE OF RANDOM EFFECTS *****/
  /***** ======================== *****/

  /*** Compute invsigscale2's ***/
  invsigscale_d     = 1/(gg_d->scale(0) * gg_d->sigma(0));
  invsigscale_b     = 1/(gg_b->scale(0) * gg_b->sigma(0));
  invsigscale2_d    = gg_d->invscale2(0)    * gg_d->invsigma2(0);
  invsigscale2_b    = gg_b->invscale2(0)    * gg_b->invsigma2(0);
  invsigscale2_zeta = gg_zeta->invscale2(0) * gg_zeta->invsigma2(0);
  invsigscale2_eps  = gg_eps->invscale2(0)  * gg_eps->invsigma2(0);

  /*** Compute one_rho2 and inv_one_rho2 ***/
  if (rho_zb[0] > rhoNorm::rhoONE || rho_zb[0] < rhoNorm::rhoMinONE)
    throw returnR("Trap in classRandomEff_rhoNorm.cpp: Gspl_rho_intcpt_update. Value of rho is too close to +-1", 1);
  one_rho2     = 1 - rho_zb[0]*rho_zb[0];
  inv_one_rho2 = 1/one_rho2;

  /*** Compute the common part of ivar_prop, see page 54 of red notes ***/
  tempP  = ivar_propCommon;
  *tempP = invsigscale2_d * inv_one_rho2;                               /* ivar_propCommon[0,0] */
  tempP++;
  *tempP = -(rho_zb[0]*inv_one_rho2) * invsigscale_d * invsigscale_b;   /* ivar_propCommon[1,0] */
  tempP++;
  *tempP = invsigscale2_b * inv_one_rho2;                               /* ivar_propCommon[1,1] */

  /*** Loop over clusters                                                                         ***/
  /*** Within the loop compute also sumu2, sumv2, sumuv needed for the update of rho afterwards   ***/
  regResOnsetP = regResOnset;
  regResTimeP  = regResTime;
  dP           = d->_bM;
  bP           = b->_bM;  
  rdP          = rM_d;
  rbP          = rM_b;
  rzetaP       = rM_zeta;
  repsP        = rM_eps;  
  nwithinClP   = d->_nwithinCl;

  sumu2 = 0;
  sumv2 = 0;
  sumuv = 0;

  for (cl = 0; cl < d->_nCluster; cl++){
  
    /*** Compute ivar_prop: inverse variance of the full conditional distribution, see p. 54 of the red notes  ***/
    tempP  = ivar_prop;
    temp2P = ivar_propCommon;
    *tempP = *temp2P + (*nwithinClP) * invsigscale2_zeta;              /* ivar_prop[0,0]  */
    tempP++;
    temp2P++;
    *tempP = *temp2P;                                                  /* ivar_prop[1,0]  */
    tempP++;
    temp2P++;
    *tempP = *temp2P + (*nwithinClP) * invsigscale2_eps;               /* ivar_prop[1,1]  */

    /*** Compute canonical mean of the full conditional distribution, see p. 55 of the red notes ***/
    /*** Part 1 comming from the prior distribution on (d, b)                                    ***/
    mean_d = gg_d->intcpt(0) + gg_d->scale(0) * mu_d[0][*rdP];
    mean_b = gg_b->intcpt(0) + gg_b->scale(0) * mu_b[0][*rbP];
    mean_prop[0] = ivar_propCommon[0]*mean_d + ivar_propCommon[1]*mean_b;
    mean_prop[1] = ivar_propCommon[1]*mean_d + ivar_propCommon[2]*mean_b;
    rdP++;
    rbP++;    

    /*** Part 2 comming from the likelihood                                                      ***/
    mean_propTemp[0] = 0.0;
    mean_propTemp[1] = 0.0;
    for (i = 0; i < *nwithinClP; i++){             /** loop over the observations in a given cluster **/
      
      /* Add old value of the random intercept to regRes                       */
      /* Compute sum(y - alpha - x'beta - scale*mu), store it in mean_PropTemp[0] */
      *regResOnsetP += (*dP);
      mean_propTemp[0] += (*regResOnsetP) - (gg_zeta->intcpt(0) + gg_zeta->scale(0)*mu_zeta[0][*rzetaP]);
      regResOnsetP++;
      rzetaP++;

      *regResTimeP += (*bP);
      mean_propTemp[1] += (*regResTimeP) - (gg_eps->intcpt(0) + gg_eps->scale(0)*mu_eps[0][*repsP]);
      regResTimeP++;
      repsP++;
    }                                              /** end of the loop over observations in a given cluster **/
    mean_prop[0] += invsigscale2_zeta * mean_propTemp[0];
    mean_prop[1] += invsigscale2_eps  * mean_propTemp[1];

    /** Sample new value of the random intercepts (d, b) **/
    AK_BLAS_LAPACK::chol_dpptrf(ivar_prop, _TWO_, info);
    if (*info) throw returnR("Trap in classRandomEff_rhoNorm.cpp: Gspl_rho_intcpt_update. Singular covariance matrix of the full conditional distribution of the random effects", 1);
    Mvtdist3::rmvnormC2006(mean_propTemp, mean_prop, ivar_prop, _TWO_);
    *dP = mean_propTemp[0];
    *bP = mean_propTemp[1];

    /** Update sumu2, sumv2, sumuv **/
    u = (*dP - mean_d)*invsigscale_d;
    v = (*bP - mean_b)*invsigscale_b;
    sumu2 += u*u;
    sumv2 += v*v;
    sumuv += u*v;

    /** Update regResOnset and regResTime **/
    regResOnsetP -= *nwithinClP;
    for (i = 0; i < *nwithinClP; i++){
      *regResOnsetP -= (*dP);
      regResOnsetP++;
    }
    dP++;

    regResTimeP -= *nwithinClP;
    for (i = 0; i < *nwithinClP; i++){
      *regResTimeP -= (*bP);
      regResTimeP++;
    }
    bP++;

    nwithinClP++;
  }    /*** end of the loop over clusters ***/


  /***** UPDATE OF THE CORRELATION COEFFICIENT RHO *****/
  /***** ========================================= *****/
  rhoNorm::update_pUnif(rho_accept, rho_zb+1, rho_zb, ll_dll_ddll, &sumu2, &sumv2, &sumuv, &(d->_nCluster), rho_algor, rho_scaleL);

  return;
}

