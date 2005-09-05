// Functions to update label components for the use with G-splines
//
// 19/10/2004: 'update_Alloc_GS_bak' (removed version)
// 02/11/2004: 'update_Alloc_GS'
// 05/11/2004: 'update_Alloc_GS' extended to return the log-likelihood
// 10/12/2004: 'update_Alloc_GS' extended for the second specification 
//
#include "update_Alloc_GS.h"

using namespace std;

//
// rM ....................... labels for each observation (integers from 0, ..., gg.total_length()-1)         [nP]
// mixtureNM ................ numbers of observations belonging to each mixture component                     [gg.total_length()]
// mu ....................... on OUTPUT: computed knots                                                       [gg.dim(), gg.length(j)]
// loglik ................... on OUTPUT: log-likelihood
// logpr .................... on OUTPUT: log(prod_{i=1}^N p(r|w))
// gg ....................... G-spline
// regresResM[nP*_dim] ...... * in the case of i.i.d data, these are observations y_1,...,y_n,
//                            * if dim(y_i) = d > 1 then it is assumed that they are stored as y_1[1],...,y_1[d], ..., y_n[1], ..., y_n[d],
//                              i.e. y_i[j] = regresResM[i*d + j]
//                            * in the case of regression, these are y_i - beta'x_i - b_i'z_i
//                              (where beta does not contain the intercept term)
//                              - these residuals are not devided by the scale!!!
// nP ....................... number of observational vectors                
// iwork .................... working array                                                                   [gg.total_length()]
// dwork .................... working array                                                                   [gg.total_length()]
//
void
update_Alloc_GS(int* rM,            int* mixtureNM,            double** mu,    double* loglik,  double* logpr,
                const Gspline* gg,  const double* regresResM,  const int* nP,   
                int* iwork,         double* dwork)
{
  int i, j, ia, row, col, obs;

  /***   Reset mixtureNM and loglik  ***/
  for (ia = 0; ia < gg->total_length(); ia++) mixtureNM[ia] = 0;
  *loglik = -(gg->dim() * LOG_SQRT_2PI);

  /*** Compute knots and add sigma part to the log-likelihood***/
  for (j = 0; j < gg->dim(); j++){
    *loglik -= log(gg->sigma(j));
    gg->muP(j, mu[j]);
  }
  (*loglik) *= (*nP);

  /*** Loop over observations ***/ 
  int k_obs, label;
  double max_logw;
  double y_mu;

  int* rp = rM;
  const double* regRes = regresResM;
  for (obs = 0; obs < *nP; obs++){
    max_logw = -FLT_MAX;

    /* Compute log-weights */
    switch (gg->dim()){
    case 1:
      for (ia = 0; ia < gg->total_length(); ia++){
        y_mu = *regRes - gg->intcpt(0) - gg->scale(0)*mu[0][ia];
        dwork[ia] = -0.5 * (gg->invsigma2(0)*gg->invscale2(0) * y_mu * y_mu) + gg->a(ia);
        if (dwork[ia] > max_logw) max_logw = dwork[ia];
      }
      break;
    case 2:
      for (col = 0; col < gg->length(1); col++){
	i = col*gg->length(0);
        for (row = 0; row < gg->length(0); row++){
          ia = i + row;
          y_mu = *regRes - gg->intcpt(0) - gg->scale(0)*mu[0][row];
          dwork[ia] = -(gg->invsigma2(0)*gg->invscale2(0) * y_mu * y_mu);
          y_mu = *(regRes+1) - gg->intcpt(1) - gg->scale(1)*mu[1][col];
          dwork[ia] -= (gg->invsigma2(1)*gg->invscale2(1) * y_mu * y_mu);
          dwork[ia] *= 0.5;
          dwork[ia] += gg->a(ia);
          if (dwork[ia] > max_logw) max_logw = dwork[ia];
        }
      }
      break;
    default:
      throw returnR("Unimplemented dimension appeared in update_Alloc_GS", 1);
    }

    /* Rescale log-weights such that the highest one will be equal to zero                                                    */
    /*  take now only weights that are higher than 'null' weights -> put cumsum(exp(log-weights)) at the beginning of dwork   */
    k_obs = 0;
    ia = 0;
    while (k_obs == 0){
      dwork[ia] -= max_logw;
      if (dwork[ia] > gg->log_null_w()){
        dwork[0] = exp(dwork[ia]);                                  /* Note that 0 is always <= ia */
        iwork[0] = ia;
        k_obs++;
      }
      ia++;
    }
    while (ia < gg->total_length()){
      dwork[ia] -= max_logw;
      if (dwork[ia] > gg->log_null_w()){
        dwork[k_obs] = dwork[k_obs-1] + exp(dwork[ia]);             /* Note that k_obs is always <= ia */
        iwork[k_obs] = ia;
        k_obs++;
      }
      ia++;
    }
    if (k_obs == 0){                               /* This is impossible when log_null_w() < 0 */
      throw returnR("C++ Error: Trap in update_Alloc_GS, all weights are zero for some observation", 1);
    }

    /* Sample new label and fill in new values into appropriate arrays */
    discreteSampler2(&label, dwork, &k_obs, &ONE_INT, &ONE_INT, &ZERO_INT);
    *rp = iwork[label];
    mixtureNM[*rp]++;

    /* Update log-likelihood */
    switch (gg->dim()){
    case 1:
      y_mu = *regRes - gg->intcpt(0) - gg->scale(0)*mu[0][*rp];
      (*loglik) -= 0.5 * (gg->invsigma2(0)*gg->invscale2(0) * y_mu * y_mu);
      regRes++;
      break;
    case 2:
      row = *rp % gg->length(0);
      col = *rp / gg->length(0);      
      y_mu = *regRes - gg->intcpt(0) - gg->scale(0)*mu[0][row];
      dwork[0] = -(gg->invsigma2(0)*gg->invscale2(0) * y_mu * y_mu);
      y_mu = *(regRes+1) - gg->intcpt(1) - gg->scale(1)*mu[1][col];
      dwork[0] -= (gg->invsigma2(1)*gg->invscale2(1) * y_mu * y_mu);
      dwork[0] *= 0.5;
      (*loglik) += dwork[0];
      regRes += 2;
      break;
    default:
      throw returnR("Unimplemented dimension appeared in update_Alloc_GS", 1);
    }

    rp++;
  }           /* end of the loop over observations */

  /* Update log(p(r|w))    */
  *logpr = -(*nP) * log(gg->sumexpa());
  for (ia = 0; ia < gg->total_length(); ia++) (*logpr) += (mixtureNM[ia] ? mixtureNM[ia]*gg->a(ia) : 0.0);

  return;
}

//
// OLDER VERSION (currently not used and most likely no more working):
//
// rM ............ labels for each observation (integers from 0, ..., gg.total_length()-1)         [nP]
// mixtureNM ..... numbers of observations belonging to each mixture component                     [gg.total_length()]
// gg ............ G-spline
// regresResM .... observations                                                                    [gg.dim() * nP]
// nP ............ number of observational vectors                
// iwork ......... working array                                                                   [gg.total_length()]
// dwork ......... working array                                                                   [gg.total_length()]
//
void
update_Alloc_GS_bak(int* rM,            int* mixtureNM,            
                    const Gspline* gg,  const double* regresResM,  const int* nP,   
                    int* iwork,         double* dwork)
{
  int i, j, ia, obs;

  /***   Reset mixtureNM   ***/
  for (ia = 0; ia < gg->total_length(); ia++) mixtureNM[ia] = 0;

  /*** Loop over observations ***/ 
  int k_obs, label;
  double y_mu;
  for (obs = 0; obs < *nP; obs++){

    /* Compute weights, take only these that are higher than 'null' weights */
    k_obs = 0;
    for (i = 0; i < gg->k_effect(); i++){
      dwork[k_obs] = 0.0;
      for (j = 0; j < gg->dim(); j++){
        y_mu = regresResM[obs*gg->dim()+j]-gg->mu_component(j,gg->ind_w_effect(i));
        dwork[k_obs] -= gg->invsigma2(j) * y_mu * y_mu;
      }
      dwork[k_obs] *= 0.5; 
      dwork[k_obs] += gg->a(gg->ind_w_effect(i));
      if (dwork[k_obs] > gg->log_null_w()){
        dwork[k_obs] = exp(dwork[k_obs]);
        iwork[k_obs] = gg->ind_w_effect(i);
        k_obs++;
      }
    }
    if (k_obs == 0){
      throw returnR("C++ Error: Trap in update_Alloc_GS, all weights are zero for some observation", 1);
    }

    /* Sample new label and fill in new values into appropriate arrays */
    discreteSampler(&label, dwork, &k_obs, &ONE_INT, &ZERO_INT, &ZERO_INT);
    rM[obs] = iwork[label];
    mixtureNM[rM[obs]]++;
  }
  return;
}
