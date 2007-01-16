/*** structRandomEff32.cpp ***/
//
//     AUTHOR:  Arnost Komarek (my name in TeX: Arno\v{s}t Kom\'arek)
//              akom@email.cz
//
//     CREATED:  11/12/2006
//
//     PURPOSE:  Structure to hold random intercepts for the bayessurvreg2, version 32 
//
//                 init:  12/12/2006
//               update:  12/12/2006
//   updateAfterChangeD:  13/12/2006
//        readDfromFile:  13/12/2006
//           predict_db:  13/12/2006
//
//

#include "structRandomEff32.h"

namespace RandomEff32 {

/*** =================================================== ***/
/*** init:  Initialize the structure                     ***/
/***                                                     ***/
/*** =================================================== ***/
//
//  data:         Initialized structure
//
//  dVal:         Initial values of the onset random intercept
//                See  'priorb1D' argument of bayessurvreg2 function  
//
//  bVal:         Initial values of the time-to-event random intercept
//                See  'priorb2D' argument of bayessurvreg2 function
//
//  parD[4]:      parD[0,1,2] = initial value of var(d,b) = D (lower triangle)
//                parD[3]     = prior degrees of freedom of the Wishart prior
//                parD[4,5,6] = prior scale matrix of the Wishart prior (lower triangle)
//
//  pardI:        Integer parameters for the onset random intercept
//                See  'priorb1I' argument of bayessurvreg2 function
//               parI[0]     = type of prior (0 = Normal, 1 = Gspline), it MUST BE 0
//               parI[1]     = number of random effects, it MUST BE 1
//               parI[2]     = number of clusters
//     parI[3,...2+nCluster] = numbers of observations within each cluster
//
//  parbI:        Integer parameters for the time-to-event random intercept
//                Structure the same as for pardI
//                
//
void
init(RandomEff32::RE *data,  double *dVal,   double *bVal,  double *parD,  const int *pardI,  const int *parbI)
{
  int i, info;
  const int *nCld, *nClb;
  const double *cdP;
  double *dP;

  /*** Type of the distribution of random effects ***/
  if (pardI[0] != 0 || parbI[0] != 0){
    throw returnR("Error in structRandomEff32.cpp: init. Type of prior of random effects must me 0 (normal).", 1);
  }

  /*** Dimension of the random effects ***/
  if (pardI[1] != 1 || parbI[1] != 1){
    throw returnR("Error in structRandomEff32.cpp: init. There must be exactly 1 random effect in each part of the model.", 1);
  }
  data->_nRandom = pardI[1] + parbI[1];
  data->_lD      = (data->_nRandom * (data->_nRandom + 1))/2;

  /*** Number of clusters ***/
  if (pardI[2] <= 0 || parbI[2] <= 0 || pardI[2] != parbI[2]){
    throw returnR("Error in structRandomEff32.cpp: init. Number of clusters must be positive and the same in both parts of the model.", 1);
  }
  data->_nCluster = pardI[2];  

  /*** Numbers of observations within each cluster ***/
  nCld = pardI + 3; 
  nClb = parbI + 3; 
  for (i = 0; i < data->_nCluster; i++){
    if (*nCld <= 0 || *nClb <= 0 || *nCld != *nClb){
      throw returnR("Error in structRandomEff32.cpp: init. Numbers of observations within each clusters must be positive and the same in both part sof the model.", 1);
    }
    nCld++;
    nClb++;
  }
  data->_nwithinCl = pardI + 3;

  /*** Values of random effects ***/
  data->_d = dVal;
  data->_b = bVal;

  /*** Value of the covariance matrix of random effects ***/
  data->_D = parD;

  /*** Covariance matrix -> its determinant and inverse ***/
  cdP = data->_D;
  dP  = data->_Di; 
  for (i = 0; i < data->_lD; i++){
    *dP = *cdP;
    dP++;
    cdP++;
  }
  AK_BLAS_LAPACK::chol_dpptrf(data->_Di, &data->_nRandom, &info);
  if (info){
    throw returnR("Error in structRandomEff32.cpp: init. Initial covariance matrix is not positive definite.", 1);
  }
  data->_detD = data->_Di[0] * data->_Di[0] * data->_Di[2] * data->_Di[2];
  AK_BLAS_LAPACK::chol_dpptri(data->_Di, &data->_nRandom, &info);

  /*** Parameters of the prior of the covariance matrix of random effects ***/
  /** Degrees of freedom **/
  if (parD[3] <= data->_nRandom - 1){
    throw returnR("Error in structRandomEff32.cpp: init. Prior Wishart degrees of freedom must be higher than 1.", 1);
  }
  data->_priorDF = parD[3];

  /** Scale matrix -> inverse scale matrix **/
  cdP = parD + 4;
  dP  = data->_priorSi; 
  for (i = 0; i < data->_lD; i++){
    *dP = *cdP;
    dP++;
    cdP++;
  }
  AK_BLAS_LAPACK::chol_dpptrf(data->_priorSi, &data->_nRandom, &info);
  if (info){
    throw returnR("Error in structRandomEff32.cpp: init. Prior Wishart scale matrix is not positive definite.", 1);
  }
  AK_BLAS_LAPACK::chol_dpptri(data->_priorSi, &data->_nRandom, &info);

  /** Degrees of freedom of the full conditional of the covariance matrix of random effects **/
  data->_propDF = data->_priorDF + data->_nCluster;

  return;
}


/*** =================================================== ***/
/*** update:  Update the members of the structure        ***/
/***                                                     ***/
/***  1) update (d, b)                                   ***/
/***  2) update var(d, b)                                ***/
/***                                                     ***/
/*** =================================================== ***/
//
// regResOnset[nP]:     regression residuals for onset
// regResTime[nP]:      regression residuals for time-to-event
//
// nP[1]:               total number of observations (over all clusters)
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
update(RandomEff32::RE *data,  
       double *regResOnset,     double *regResTime,
       const int *nP,
       const Gspline *gg_zeta,  double** const mu_zeta,    const int *rM_zeta,
       const Gspline *gg_eps,   double** const mu_eps,     const int *rM_eps)
{
  static int info[1];

  static int cl, i;
  static double invsigscale2_zeta, invsigscale2_eps;
  static double *sumd2, *sumb2, *sumdb;
  static double *tempP, *temp2P;
  static double *regResOnsetP, *regResTimeP;
  static double *propMean_d, *propMean_b;
  static double *dP, *bP;
  static const double *cdP;
  static const int *rzetaP, *repsP, *nwithinClP;  

  /***** UPDATE OF RANDOM EFFECTS *****/
  /***** ======================== *****/

  /*** Compute invsigscale2's ***/
  invsigscale2_zeta = gg_zeta->invscale2(0) * gg_zeta->invsigma2(0);
  invsigscale2_eps  = gg_eps->invscale2(0)  * gg_eps->invsigma2(0);

  /*** Loop over clusters                                                                         ***/
  /*** Within the loop compute also sumd2, sumb2, sumbd needed for the update of D afterwards     ***/
  regResOnsetP = regResOnset;
  regResTimeP  = regResTime;
  dP           = data->_d;
  bP           = data->_b;  
  rzetaP       = rM_zeta;
  repsP        = rM_eps;  
  nwithinClP   = data->_nwithinCl;

  sumd2 = data->_propSi;
  sumdb = sumd2 + 1;
  sumb2 = sumdb + 1;
  *sumd2 = 0;
  *sumdb = 0;
  *sumb2 = 0;

  propMean_d = data->_propMean;
  propMean_b = propMean_d + 1;

  for (cl = 0; cl < data->_nCluster; cl++){

    /*** Compute the inverse variance of the full conditional distribution, see page 57 of red notes ***/
    tempP  = data->_propVar;
    temp2P = data->_Di;
    *tempP = *temp2P + (*nwithinClP) * invsigscale2_zeta;              /* _propVar[0,0]  */
    tempP++;
    temp2P++;
    *tempP = *temp2P;                                                  /* _propVar[1,0]  */
    tempP++;
    temp2P++;
    *tempP = *temp2P + (*nwithinClP) * invsigscale2_eps;               /* _propVar[1,1]  */

    /*** Compute canonical mean of the full conditional distribution, see p. 57 of the red notes ***/
    /*** Part comming from the likelihood                                                        ***/    
    *propMean_d = 0.0;
    *propMean_b = 0.0;
    for (i = 0; i < *nwithinClP; i++){             /** loop over the observations in a given cluster **/
      
      /* Add old value of the random intercept to regRes                       */
      /* Compute sum(y - alpha - x'beta - scale*mu), store it in _propMean     */
      *regResOnsetP += (*dP);
      *propMean_d   += (*regResOnsetP) - (gg_zeta->intcpt(0) + gg_zeta->scale(0)*mu_zeta[0][*rzetaP]);
      regResOnsetP++;
      rzetaP++;

      *regResTimeP += (*bP);
      *propMean_b  += (*regResTimeP) - (gg_eps->intcpt(0) + gg_eps->scale(0)*mu_eps[0][*repsP]);
      regResTimeP++;
      repsP++;
    }                                              /** end of the loop over observations in a given cluster **/
    *propMean_d *= invsigscale2_zeta;
    *propMean_b *= invsigscale2_eps;

    /** Sample new value of the random intercepts (d, b) in given cluster **/
    AK_BLAS_LAPACK::chol_dpptrf(data->_propVar, &data->_nRandom, info);
    if (*info) throw returnR("Trap in structRandomEff32.cpp: update. Singular covariance matrix of the full conditional distribution of the random effects", 1);
    Mvtdist3::rmvnormC2006(data->_propValue, data->_propMean, data->_propVar, &data->_nRandom);
    *dP = data->_propValue[0];
    *bP = data->_propValue[1];

    /** Update sumd2, sumb2, sumdb **/
    *sumd2 += (*dP)*(*dP);
    *sumb2 += (*bP)*(*bP);
    *sumdb += (*dP)*(*bP);

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
  }


  /***** UPDATE OF THE COVARIANCE MATRIX OF RANDOM EFFECTS *****/
  /***** ================================================= *****/

  /*** Inverse scale matrix of the Wishart full conditional ***/
  tempP  = data->_propSi;
  temp2P = data->_priorSi;
  *tempP = *temp2P + (*sumd2);                                          /* _propSi[0,0]  */
  tempP++;
  temp2P++;
  *tempP = *temp2P + (*sumdb);                                          /* _propSi[1,0]  */
  tempP++;
  temp2P++;
  *tempP = *temp2P + (*sumb2);                                          /* _propSi[1,1]  */

  /*** Sample from the Wishart distribution ***/
  Mvtdist3::rwishart3(data->_Di, data->_workWishart, &data->_propDF, data->_propSi, &data->_nRandom, 1);

  /*** Inverse covariance matrix -> covariance matrix and its determinant ***/
  cdP = data->_Di;
  dP  = data->_D; 
  for (i = 0; i < data->_lD; i++){
    *dP = *cdP;
    dP++;
    cdP++;
  }
  AK_BLAS_LAPACK::chol_dpptrf(data->_D, &data->_nRandom, info);
  if (*info){
    throw returnR("Error in structRandomEff32.cpp: update. Sampled covariance matrix is not positive definite.", 1);
  }
  data->_detD = 1/(data->_D[0] * data->_D[0] * data->_D[2] * data->_D[2]);
  AK_BLAS_LAPACK::chol_dpptri(data->_D, &data->_nRandom, info);

  return;
}


/*** =========================================================================== ***/
/*** updateAfterChangeD:  Update _Di and _detD after _D has been changed         ***/
/***                                                                             ***/
/*** =========================================================================== ***/
void
updateAfterChangeD(RandomEff32::RE *data)
{
  static const double *cdP;
  static double *dP;
  static int i;
  static int info[1];

  /*** Covariance matrix -> inverse covariance matrix and its determinant ***/
  cdP = data->_D;
  dP  = data->_Di; 
  for (i = 0; i < data->_lD; i++){
    *dP = *cdP;
    dP++;
    cdP++;
  }
  AK_BLAS_LAPACK::chol_dpptrf(data->_Di, &data->_nRandom, info);
  if (*info){
    throw returnR("Error in structRandomEff32.cpp: updateAfterChangeD. Covariance matrix is not positive definite.", 1);
  }
  data->_detD = data->_Di[0] * data->_Di[0] * data->_Di[2] * data->_Di[2];
  AK_BLAS_LAPACK::chol_dpptri(data->_Di, &data->_nRandom, info);

  return;
}


/*** =========================================================================== ***/
/*** predict_db:  Sample from N(0, D) for each cluster                           ***/
/***                                                                             ***/
/*** =========================================================================== ***/
void
predict_db(RandomEff32::RE *data)
{
  static const double *cdP;
  static double *dP, *bP;
  static int i, cl;
  static int info[1];

  /*** Covariance matrix -> Cholesky decomposition ***/
  cdP = data->_D;
  dP  = data->_propVar; 
  for (i = 0; i < data->_lD; i++){
    *dP = *cdP;
    dP++;
    cdP++;
  }
  AK_BLAS_LAPACK::chol_dpptrf(data->_propVar, &data->_nRandom, info);
  if (*info){
    throw returnR("Error in structRandomEff32.cpp: predict_db. Covariance matrix is not positive definite.", 1);
  }

  /*** Mean ***/
  data->_propMean[0] = 0;
  data->_propMean[1] = 0;

  /*** Sample ***/
  dP = data->_d;
  bP = data->_b;
  for (cl = 0; cl < data->_nCluster; cl++){
    Mvtdist3::rmvnorm2006(data->_propValue, data->_propMean, data->_propVar, &data->_nRandom);
    *dP = data->_propValue[0];
    *bP = data->_propValue[1];
    dP++;
    bP++;
  }

  return;
}

}  /** end of the namespace RandomEff32 **/


