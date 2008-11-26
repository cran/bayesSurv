// Main function to run MCMC for an AFT model with random effects
//  and univariate G-spline error
// * allow for doubly censoring
// * distribution of random effects: either NORMAL
//                                   or G-spline if only random intercept in the model
//
// 27/01/2005:  start working on it
// 31/01/2005:  version for normal random effects finished
// 02/02/2005:  extended to handle also G-spline random intercept
// 20/02/2005:  'adjust_intcpt' function added and adjustment of the G-spline intercept terms
//              included in the case when there is G-spline random intercept in the model
//              * immediately again removed due to bad performance
// 25/04/2005:  option mainSimul added to write simulated values during burn-in only after nwrite iteration
// 28/11/2006:  implementing simultaneous update of all 'a' coefficients
// 06/12/2006:  implementing version 31 (see below)
// 11/12/2006:  implementing version 32 (see below)
// 26/11/2008:  * adjustments to allow for a model with (random) intercept only
//              * previous versions took take this info in beta1->randomIntcpt() and beta2->randomIntcpt()
//                which was always set to 0 if there are no covariates in the model  
//              * this caused SEGFAULT in RandomEff::GIBBSupdate
//
#include "bayessurvreg2.h"

extern "C"{

// dimsP[] ................ dimensionality parameters
//    dimsP[0]                   = nP (sample size, number of observations)
//    dimsP[1]                   = doubly, 0/1 indiating whether we have doubly censored observations

// priorb1I[] ............. integer parameters for constructor of the random effects object
//    priorb1I[0]                    = type of prior for b (0 = NORMAL)
//                                                         (1 = G-spline)
//    priorb1I[1]                    = nRandom 
//    priorb1I[2]                    = nCluster
//    priorb1I[3, ..., 3+nCluster-1] = nwithinCl
//
// For version 32: priorb1I[0] = 0
//                 priorb1I[1] = 1
//
//                 priorb2I[0] = 0
//                 priorb2I[1] = 1

// priorb1D[] ............. double parameters for constructor of the random effects object (RandomEff class or RandomEff32 structure)
//    priorb1D[0, ..., nRandom*nCluster-1] = initial values for random effects for onset
//
// For version 32: priorb1D[0,...,nCluster-1] = initial values of the random intercept for onset
//                 priorb2D[0,...,nCluster-1] = initial values of the random intercept for time-to-event

// priorCovMat
// ============
// 1) in version for normal random effects, these are arguments for the constructor of the class CovMatrix
// 2) in version for G-spline random intercept, these are arguments for the constructor of the class Gspline
//    related to the distribution of the random intercept
//
//    1) version with normal random effects:
//    priorCovMat1I[] ........ integer parameters for constructor of the covariance matrix of random effects
//       priorCovMat1I[0] = nRandom 
//       priorCovMat1I[1] = type of prior for D (0 = Inverse Wishart)
//                                              (1 = SDUniform)
//
//    For version 32:  priorCovMat1I[], priorCovMat2I[] = arbitrary
//
//    priorCovMat1D[] ........ double parameters for constructor of the covariance matrix of random effects
//       priorCovMat1D[0, ..., 0.5*nRandom*(nRandom+1)-1] = initial value of D (lower triangle)
//       priorCovMat1D[0.5*nRandom*(nRandom+1)]           = dfD
//       priorCovMat1D[]                                  = scaleD (lower triangle)
//
//    For version 32: priorCovMat1D[] = as above
//                    priorCovMat2D[] = arbitrary
//
//    2) version with G-spline random intercept
//    priorCovMatI1[] ......... integer parameters for constructor of the G-spline defining the distribution of random effects
//
//    priorCovMatD1[] ......... double parameters for constructor of the G-spline defining the distribution of random effects
//
// rhob 
// ====
//    rhob[1]:         INPUT:   initial value of rho in the version = 31
//                     OUTPUT:  last sampled value of rho in the version = 31
//    rho_accept[]:    OUTPUT:  acceptance indicators for rho (correlation coefficient)
//    rhobI[1]:        Type of the algorithm that is used to update rho
//                     see 'rhoNorm.h' for possible values
//    rhobD[1]:        Scale parameter of the proposal for the case rhobI[0] = rhoNorm::_Langevin_      
//
// ----------------------------------------------------------------------------------------------------------------------------

// storeP[8 or 12] ............ parameters defining what should be stored (a1, y1, r1, b1, a2, y2, r2, b2, a_b, r_b, a_b2, r_b2)
//                              (last four elements are only present when the version for G-spline random intercept is in use)

// version .................... 2 or 3, specifying which version of bayessurvreg called this C++ routine from R
//    2  = bayessurvreg2, normal random effects
//    3  = bayessurvreg3, G-spline random effects
//    31 = bayessurvreg3, doubly interval-censored data, d = univariate random effect for onset
//                                                       b = univariate random effect for time-to-event
//                        (b,d) ~ G-spline with basis having the covariance matrix D = (d[j,k]),
//                                d[1,1] = sigma[1]^1 ...... fixed basis variance in the first margin                             
//                                d[2,2] = sigma[2]^1 ...... fixed basis variance in the second margin
//                                d[1,2] = rho*sigma[1]*sigma[2], where rho is unknown correlation coefficient
//                                         apriori rho ~ Unif(-1, 1)
//                                         * sampling is done on z = -0.5*log((1-rho)/(1+rho)) = atanh(rho)
//   32 = bayessurvreg3, doubly interval-censored data, d = univariate random effect for onset
//                                                      b = univariate random effect for time-to-event
//                        (b,d) ~ Normal with zero mean and covariance matrix D (general)
//                                         apriori D ~ Wishart(nu, S)
//
// mainSimul .................. 0/1, if 0 then all parameters are written to files only every nwrite iteration

void
bayessurvreg2(char **dirP,             const int *dimsP,         const double *X1,     const double *X2,
              const double *y1_left,   const double *y1_right,   const int *status1,
              const double *t2_left,   const double *t2_right,   const int *status2,               
              double *Y1,              double *Y2,
              int *r1,                 int *r2,                  const int *specif,    
              int *r_b1,               int *r_b2,                const int *specif_b,
              int *GsplineI1,          double *GsplineD1,
              int *GsplineI2,          double *GsplineD2,
              int *priorBeta1I,        double *priorBeta1D,
              int *priorBeta2I,        double *priorBeta2D,
              int *priorb1I,           double *priorb1D,
              int *priorb2I,           double *priorb2D,
              int *priorCovMat1I,      double *priorCovMat1D,
              int *priorCovMat2I,      double *priorCovMat2D,              
              double *rhob,            int *rho_accept,          const int *rhobI,     const double *rhobD,
              int *iterM,              int *nsimulP,             int *storeP,
              const int *version,      const int *mainSimul,     int *errP)
{
  try{
    int out_format[2] = {6, 1};                /** precision and width for output **/
    double dtemp;
    int itemp;
    double *ddtemp = (double*)malloc(sizeof(double));

    const int check_k_effect = 0;    /** option for writeToFiles_bayesHistogram function **/

    int i, j, k;
    GetRNGstate();

    *errP = 0;
    std::string dir = *dirP;
    switch (*version){
    case 2:
    case 3:
    case 31:
    case 32:
      break;

    default:
      throw returnR("Error: unimplemented version appeared in bayessurvreg2", 1);
    }

    /** Numbers of iterations etc.  **/
    int niter = nsimulP[0]; 
    int nthin = nsimulP[1];
    int nwrite = nsimulP[2];


    /** Mandatory storing of simulated values **/
    int *storea1P = storeP;
    int *storey1P = storeP + 1;
    int *storer1P = storeP + 2;
    int *storeb1P = storeP + 3;
    int *storea2P = storeP + 4;
    int *storey2P = storeP + 5;
    int *storer2P = storeP + 6;
    int *storeb2P = storeP + 7;
    int storea_b1P, storer_b1P, storea_b2P, storer_b2P;
    switch (*version){
    case 3:
    case 31:
    case 32:
      storea_b1P = storeP[8];
      storer_b1P = storeP[9];
      storea_b2P = storeP[10];
      storer_b2P = storeP[11];
      break;
    }

    /** Dimensionality parameters  **/        
    const int *nP        = dimsP;                                      /* number of observational uni- or bivariate vectors  */
    const int *doubly    = dimsP + 1;                                  /* 0/1 .... is doubly censored?                       */


    /** Dimensionality parameters for compatibility with bayesBisurvreg **/
    const int *dimP   = GsplineI1 + 0;                                 /* dimension of the response                          */
    if (*dimP != 1) throw returnR("Error: Dimension of the G-spline must be 1 (bayessurvreg2; 1)", 1);
    if (*doubly){
      k = GsplineI2[0];
      if (k != 1) throw returnR("Error: Dimension of the G-spline must be 1 (bayessurvreg2; 2)", 1);
    }
    const int nobs    = *nP;                                           /* number of observations                             */
    if (!(*doubly)) *storea2P = *storey2P = *storer2P = *storeb2P = 0;
    if ((*version == 3) && !(*doubly)) storea_b2P = storer_b2P = 0;
    if ((*version == 31) && !(*doubly))
      throw returnR("Error: Version 31 of Cpp bayessurvreg2 expects doubly-interval-censored data", 1);
    if ((*version == 32) && !(*doubly))
      throw returnR("Error: Version 32 of Cpp bayessurvreg2 expects doubly-interval-censored data", 1);


    /** Find out how many censored observations we have **/
    int n1_interval = 0;
    int n1_censored = 0;
    for (i = 0; i < nobs; i++){
      switch (status1[i]){
      case 1:
        break;
      case 0:
      case 2:
        n1_censored++;
        break;
      case 3:
        n1_censored++;
        n1_interval++;
        break;
      default:
        throw returnR("Incorrect status indicator supplied", 1);
      }
    }
    if (!n1_censored) *storey1P = 0;

    int n2_interval = 0;
    int n2_censored = 0;
    if (*doubly){
      for (i = 0; i < nobs; i++){
        switch (status2[i]){
        case 1:
          break;
        case 0:
        case 2:
          n2_censored++;
          break;
        case 3:
          n2_censored++;
          n2_interval++;
          break;
        default:
          throw returnR("Incorrect status indicator supplied", 1);
        }
      }
    }


    /** Objects for beta1 and beta2          **/
    const int *RandomIntcpt1 = priorBeta1I + 3;       /** indication of random intercept in model   **/
    const int *RandomIntcpt2 = priorBeta2I + 3;       /** introduced on 26/11/2008                  **/

    BetaGammaExtend *beta1 = new BetaGammaExtend;
    BetaGammaExtend *beta2 = new BetaGammaExtend;
    if (!beta1) throw returnR("Not enough memory available in bayessurvreg2 (beta1)", 1);
    if (!beta2) throw returnR("Not enough memory available in bayessurvreg2 (beta2)", 1);
    //if (priorBeta1I[0])             *beta1 = BetaGammaExtend(priorBeta1I, priorBeta1D);   // commented on 26/11/2008
    //if (*doubly && priorBeta2I[0])  *beta2 = BetaGammaExtend(priorBeta2I, priorBeta2D);   // commented on 26/11/2008
    *beta1 = BetaGammaExtend(priorBeta1I, priorBeta1D);                                     // added on 26/11/2008
    if (*doubly) *beta2 = BetaGammaExtend(priorBeta2I, priorBeta2D);                        // added on 26/11/2008

    /** Object for b1 & D1/random intcpt G-spline  and b2 & D2/random intcpt G-spline   **/
    RandomEff *bb1       = new RandomEff;
    RandomEff *bb2       = new RandomEff;

    CovMatrixExtend *DD1 = new CovMatrixExtend;
    CovMatrixExtend *DD2 = new CovMatrixExtend;

    Gspline *g_bb1       = new Gspline;
    Gspline *g_bb2       = new Gspline;

    RandomEff32::RE *db  = new RandomEff32::RE;

    if (!bb1 || !DD1 || !g_bb1) throw returnR("Not enough memory available in bayessurvreg2 (bb1/DD1/g_bb1)", 1);
    if (!bb2 || !DD2 || !g_bb2) throw returnR("Not enough memory available in bayessurvreg2 (bb2/DD2/g_bb2)", 1);
    if (!db) throw returnR("Not enough memory available in bayessurvreg2 (db)", 1);
    double _null_weight1_b, _null_weight2_b; 
    switch (*version){
    case 2:
      if (priorb1I[1] && priorCovMat1I[0]){
        *bb1 = RandomEff(priorb1I, priorb1D);
        *DD1 = CovMatrixExtend(priorCovMat1I, priorCovMat1D);
      }    
      if (*doubly && priorb2I[1] && priorCovMat2I[0]){
        *bb2 = RandomEff(priorb2I, priorb2D);
        *DD2 = CovMatrixExtend(priorCovMat2I, priorCovMat2D);
      }
      break;

    case 3:
    case 31:
      if (priorb1I[1] && priorCovMat1I[0]){
        *bb1   = RandomEff(priorb1I, priorb1D);
        *g_bb1 = Gspline(priorCovMat1I, priorCovMat1D);
        _null_weight1_b = _null_mass/g_bb1->total_length();
      }
      if (*doubly && priorb2I[1] && priorCovMat2I[0]){
        *bb2   = RandomEff(priorb2I, priorb2D);
        *g_bb2 = Gspline(priorCovMat2I, priorCovMat2D);
        _null_weight2_b = _null_mass/g_bb2->total_length();
      }
      break;

    case 32:
      RandomEff32::init(db, priorb1D, priorb2D, priorCovMat1D, priorb1I, priorb2I);
      break;               
    }

    /** Objects for rho and z (if needed) **/
    double rho_zb[2];
    if (*version == 31){
      rho_zb[0] = *rhob;
      rhoNorm::rho2zError(rho_zb+1, rho_zb);      /* this will also check whether rho is not too close to +-1 */      
    }

    /** Regression residuals       **/
    double *regresRes1 = (double*) calloc(nobs, sizeof(double));
    if (!regresRes1) throw returnR("Not enough memory available in bayessurvreg2 (regresRes1)", 1);

    double *regresRes2 = &dtemp;
    int nCluster1, nCluster2;
    switch (*version){
    case 2:
    case 3:
    case 31:
      nCluster1 = bb1->nCluster();
      nCluster2 = bb2->nCluster();
      regresRes_GS2006(regresRes1, Y1, beta1, bb1->bMP(), X1, bb1->nwithinClP(), &nobs, &nCluster1);
      if (*doubly){
        regresRes2 = (double*) calloc(nobs, sizeof(double));
        if (!regresRes2) throw returnR("Not enough memory available in bayessurvreg2 (regresRes2)", 1);
        regresRes_GS2006(regresRes2, Y2, beta2, bb2->bMP(), X2, bb2->nwithinClP(), &nobs, &nCluster2);
      }
      break;

    case 32:
      nCluster1 = db->_nCluster;
      nCluster2 = db->_nCluster;
      regresRes_GS2006(regresRes1, Y1, beta1, db->_d, X1, db->_nwithinCl, &nobs, &nCluster1);
      if (*doubly){
        regresRes2 = (double*) calloc(nobs, sizeof(double));
        if (!regresRes2) throw returnR("Not enough memory available in bayessurvreg2 (regresRes2)", 1);
        regresRes_GS2006(regresRes2, Y2, beta2, db->_b, X2, db->_nwithinCl, &nobs, &nCluster2);
      }
      break;
    }


    /** Some manipulations with the design matrices  **/
    const double *XNow;
    double *XXtNow = &dtemp;

    double *XXtb1 = &dtemp;
    if (beta1->nbeta()){
      XXtb1 = (double*) calloc(nobs * beta1->lcovFixed(), sizeof(double));
      if (!XXtb1) throw returnR("Not enough memory available in bayessurvreg2 (XXtb1)", 1);

      XNow   = X1;
      XXtNow = XXtb1;
      for (i = 0; i < nobs; i++){
        for (j = 0; j < beta1->nFixed(); j++){
          for (k = j; k < beta1->nFixed(); k++){
            XXtNow[beta1->diagIFixed(j) + k - j] = XNow[beta1->indFixed(j)] * XNow[beta1->indFixed(k)];
          }
        }
        XNow += beta1->nbeta();
        XXtNow += beta1->lcovFixed();
      }
    }

    double *XXtb2 = &dtemp;
    if (*doubly && beta2->nbeta()){
      XXtb2 = (double*) calloc(nobs * beta2->lcovFixed(), sizeof(double));
      if (!XXtb2) throw returnR("Not enough memory available in bayessurvreg2 (XXtb2)", 1);

      XNow   = X2;
      XXtNow = XXtb2;
      for (i = 0; i < nobs; i++){
        for (j = 0; j < beta2->nFixed(); j++){
          for (k = j; k < beta2->nFixed(); k++){
            XXtNow[beta2->diagIFixed(j) + k - j] = XNow[beta2->indFixed(j)] * XNow[beta2->indFixed(k)];
          }
        }
        XNow += beta2->nbeta();
        XXtNow += beta2->lcovFixed();
      }
    }

    double *ZZtb1 = &dtemp;
    double *ZZtb2 = &dtemp;
    switch (*version){
    case 2:
      if (beta1->nRandom()){
        ZZtb1 = (double*) calloc(nobs * bb1->larray(), sizeof(double));
        if (!ZZtb1) throw returnR("Not enough memory available in bayessurvreg2 (ZZtb1)", 1);

        XNow    = X1;
        XXtNow  = ZZtb1;
        for (i = 0; i < nobs; i++){
          if (beta1->randomIntcpt()){
            XXtNow[0] = 1.0;
            for (k = 1; k < beta1->nRandom(); k++){
              XXtNow[k] = XNow[beta1->indbinXA(k)];
            }
          }
          for (j = beta1->randomIntcpt(); j < beta1->nRandom(); j++){
            for (k = j; k < beta1->nRandom(); k++){
              XXtNow[bb1->diagI(j) + k - j] = XNow[beta1->indbinXA(j)] * XNow[beta1->indbinXA(k)];
            }
          }
          XNow   += beta1->nbeta();
          XXtNow += bb1->larray();
        } 
      }

      if (*doubly && beta2->nRandom()){
        ZZtb2 = (double*) calloc(nobs * bb2->larray(), sizeof(double));
        if (!ZZtb2) throw returnR("Not enough memory available in bayessurvreg2 (ZZtb2)", 1);

        XNow    = X2;
        XXtNow  = ZZtb2;
        for (i = 0; i < nobs; i++){
          if (beta2->randomIntcpt()){
            XXtNow[0] = 1.0;
            for (k = 1; k < beta2->nRandom(); k++){
              XXtNow[k] = XNow[beta2->indbinXA(k)];
            }
          }
          for (j = beta2->randomIntcpt(); j < beta2->nRandom(); j++){
            for (k = j; k < beta2->nRandom(); k++){
              XXtNow[bb2->diagI(j) + k - j] = XNow[beta2->indbinXA(j)] * XNow[beta2->indbinXA(k)];
            }
          }
          XNow   += beta2->nbeta();
          XXtNow += bb2->larray();
        } 
      }
      break;

    case 3:
    case 31:
    case 32:
      break;
    }


    /** G-spline objects for the error **/
                              /* _null_weight ...... limit for the mixture component weight to be recorded in a file */
    Gspline *g_y1 = new Gspline;
    Gspline *g_y2 = new Gspline;
    if (!g_y1 || !g_y2) throw returnR("Not enough memory available in bayessurvreg2 (g_y1/g_y2)", 1);
    *g_y1 = Gspline(GsplineI1, GsplineD1);

    const double _null_weight1 = _null_mass/g_y1->total_length();        
    double _null_weight2;
    if (*doubly){
      *g_y2 = Gspline(GsplineI2, GsplineD2);
      _null_weight2 = _null_mass/g_y2->total_length();        
    }


    /** Arrays to store numbers of observations belonging to each mixture component  **/
    /** Subtract 1 from each rM (R -> C++ indeces) and check for consistency         **/
    int *mixtureN1 = (int*) calloc(g_y1->total_length(), sizeof(int));
    if (!mixtureN1) throw returnR("Not enough memory available in bayessurvreg2 (mixtureN1)", 1);
    for (i = 0; i < g_y1->total_length(); i++) mixtureN1[i] = 0;
    for (i = 0; i < *nP; i++){
      r1[i]--;
      if (r1[i] < 0 || r1[i] >= g_y1->total_length()) throw returnR("Inconsistent initial r1 supplied", 1);
      mixtureN1[r1[i]]++;
    }

    int *mixtureN2 = &itemp;
    if (*doubly){
      mixtureN2 = (int*) calloc(g_y2->total_length(), sizeof(int));
      if (!mixtureN2) throw returnR("Not enough memory available in bayessurvreg2 (mixtureN2)", 1);
      for (i = 0; i < g_y2->total_length(); i++) mixtureN2[i] = 0;
      for (i = 0; i < *nP; i++){
        r2[i]--;
        if (r2[i] < 0 || r2[i] >= g_y2->total_length()) throw returnR("Inconsistent initial r2 supplied", 1);
        mixtureN2[r2[i]]++;
      }
    }

    
    /** Arrays to store numbers of observations belonging to each mixture component: for G-spline random intercept **/
    /** Subtract 1 from each r_bM (R -> C++ indeces) and check for consistency                                     **/
    int *mixtureN_b1 = &itemp;
    int *mixtureN_b2 = &itemp;
    switch (*version){
    case 2:
    case 32:
      break;

    case 3:
    case 31:
      if (bb1->nRandom()){
        mixtureN_b1 = (int*) calloc(g_bb1->total_length(), sizeof(int));
        if (!mixtureN_b1) throw returnR("Not enough memory available in bayessurvreg2 (mixtureN_b1)", 1);
        for (i = 0; i < g_bb1->total_length(); i++) mixtureN_b1[i] = 0;
        for (i = 0; i < nCluster1; i++){
          r_b1[i]--;
          if (r_b1[i] < 0 || r_b1[i] >= g_bb1->total_length()) throw returnR("Inconsistent initial r_b1 supplied", 1);
          mixtureN_b1[r_b1[i]]++;
        }
      }
      if (bb2->nRandom()){
        mixtureN_b2 = (int*) calloc(g_bb2->total_length(), sizeof(int));
        if (!mixtureN_b2) throw returnR("Not enough memory available in bayessurvreg2 (mixtureN_b2)", 1);
        for (i = 0; i < g_bb2->total_length(); i++) mixtureN_b2[i] = 0;
        for (i = 0; i < nCluster2; i++){
          r_b2[i]--;
          if (r_b2[i] < 0 || r_b2[i] >= g_bb2->total_length()) throw returnR("Inconsistent initial r_b2 supplied", 1);
          mixtureN_b2[r_b2[i]]++;
        }
      }
      break;
    }

    /** Open files to store simulated values for writing **/
    std::string iterpath      = dir + "/iteration.sim";
    std::string betapath      = dir + "/beta.sim";        std::string betapath2  = dir + "/beta_2.sim"; 
    std::string bbpath        = dir + "/b.sim";           std::string bbpath2    = dir + "/b_2.sim";    
    std::string Dpath         = dir + "/D.sim";           std::string Dpath2     = dir + "/D_2.sim";

    std::string sigmapath     = dir + "/gspline.sim";
    std::string lambdapath    = dir + "/lambda.sim";      std::string mixmomentpath  = dir + "/mixmoment.sim";
    std::string mweightpath   = dir + "/mweight.sim";     std::string mlogweightpath = dir + "/mlogweight.sim";
    std::string mmeanpath     = dir + "/mmean.sim";       std::string Ypath          = dir + "/Y.sim";           
    std::string rpath         = dir + "/r.sim";           std::string logposterpath  = dir + "/logposter.sim";

    std::string sigmapath2    = dir + "/gspline_2.sim";
    std::string lambdapath2   = dir + "/lambda_2.sim";      std::string mixmomentpath2  = dir + "/mixmoment_2.sim";
    std::string mweightpath2  = dir + "/mweight_2.sim";     std::string mlogweightpath2 = dir + "/mlogweight_2.sim";
    std::string mmeanpath2    = dir + "/mmean_2.sim";       std::string Ypath2          = dir + "/Y_2.sim";           
    std::string rpath2        = dir + "/r_2.sim";           std::string logposterpath2  = dir + "/logposter_2.sim";

    std::string sigmapath_b   = dir + "/gspline_b.sim";
    std::string lambdapath_b  = dir + "/lambda_b.sim";      std::string mixmomentpath_b  = dir + "/mixmoment_b.sim";
    std::string mweightpath_b = dir + "/mweight_b.sim";     std::string mlogweightpath_b = dir + "/mlogweight_b.sim";
    std::string mmeanpath_b   = dir + "/mmean_b.sim";       std::string logposterpath_b  = dir + "/logposter_b.sim";
    std::string rpath_b       = dir + "/r_b.sim";           

    std::string sigmapath_b2   = dir + "/gspline_b2.sim";
    std::string lambdapath_b2  = dir + "/lambda_b2.sim";      std::string mixmomentpath_b2  = dir + "/mixmoment_b2.sim";
    std::string mweightpath_b2 = dir + "/mweight_b2.sim";     std::string mlogweightpath_b2 = dir + "/mlogweight_b2.sim";
    std::string mmeanpath_b2   = dir + "/mmean_b2.sim";       std::string logposterpath_b2  = dir + "/logposter_b2.sim";
    std::string rpath_b2       = dir + "/r_b2.sim";           

    std::string rhobpath       = dir + "/rho_b.sim";

    std::ofstream iterfile, betafile, betafile2, bbfile, bbfile2, Dfile, Dfile2;
    std::ofstream sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile;
    std::ofstream sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2;
    std::ofstream sigmafile_b, lambdafile_b, mixmomentfile_b, mweightfile_b, mlogweightfile_b, mmeanfile_b, rfile_b, logposterfile_b;
    std::ofstream sigmafile_b2, lambdafile_b2, mixmomentfile_b2, mweightfile_b2, mlogweightfile_b2, mmeanfile_b2, rfile_b2, logposterfile_b2;
    std::ofstream rhobfile;

    const bool writeSumExpa = false;
    std::string ssumexpa    = dir + "/sumexpa.sim";
    std::string ssumexpa2   = dir + "/sumexpa_2.sim";
    std::string ssumexpa_b  = dir + "/sumexpa_b.sim";
    std::string ssumexpa_b2 = dir + "/sumexpa_b2.sim";
    std::string smaxa    = dir + "/maxa.sim";
    std::string smaxa2   = dir + "/maxa_2.sim";
    std::string smaxa_b  = dir + "/maxa_b.sim";
    std::string smaxa_b2 = dir + "/maxa_b2.sim";
    std::ofstream sumexpafile, sumexpa2file, sumexpa_bfile, sumexpa_b2file;
    std::ofstream maxafile, maxa2file, maxa_bfile, maxa_b2file;
    if (writeSumExpa){
      openFile(sumexpafile, ssumexpa, 'o');
      openFile(maxafile, smaxa, 'o');

      if (*version == 3 || *version == 31){
        openFile(sumexpa_bfile, ssumexpa_b, 'o');
        openFile(maxa_bfile, smaxa_b, 'o');
      }

      if (*doubly){
        openFile(sumexpa2file, ssumexpa2, 'o');
        openFile(maxa2file, smaxa2, 'o');

        if (*version == 3 || *version == 31){
          openFile(sumexpa_b2file, ssumexpa_b2, 'o');
          openFile(maxa_b2file, smaxa_b2, 'o');
        }
      }
    }

    openFile(iterfile, iterpath, 'a');
    if (beta1->nbeta()) openFile(betafile, betapath, 'a');
    if (beta2->nbeta()) openFile(betafile2, betapath2, 'a');

    switch (*version){
    case 2:
      if (bb1->nRandom()){
        openFile(bbfile, bbpath, 'a');
        openFile(Dfile, Dpath, 'a');
      }
      if (bb2->nRandom()){
        openFile(bbfile2, bbpath2, 'a');
        openFile(Dfile2, Dpath2, 'a');
      }
      break;

    case 31:
      openFile(rhobfile, rhobpath, 'a');

    case 3:
      if (bb1->nRandom()){
        openFile(bbfile, bbpath, 'a');
        openFiles_bayesHistogram(sigmafile_b, lambdafile_b, mixmomentfile_b, mweightfile_b, mlogweightfile_b, mmeanfile_b, 
                                 Yfile, rfile_b, logposterfile_b,
	    		         sigmapath_b, lambdapath_b, mixmomentpath_b, mweightpath_b, mlogweightpath_b, mmeanpath_b, 
                                 Ypath, rpath_b, logposterpath_b,
                                 0, 'a');
      }
      if (bb2->nRandom()){
        openFile(bbfile2, bbpath2, 'a');
        openFiles_bayesHistogram(sigmafile_b2, lambdafile_b2, mixmomentfile_b2, mweightfile_b2, mlogweightfile_b2, mmeanfile_b2, 
                                 Yfile, rfile_b2, logposterfile_b2,
	    		         sigmapath_b2, lambdapath_b2, mixmomentpath_b2, mweightpath_b2, mlogweightpath_b2, mmeanpath_b2, 
                                 Ypath, rpath_b2, logposterpath_b2,
                                 0, 'a');
      }
      break;

    case 32:
      openFile(bbfile, bbpath, 'a');
      openFile(bbfile2, bbpath2, 'a');
      openFile(Dfile, Dpath, 'a');
      break;
    }

    openFiles_bayesHistogram(sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, 
                             Yfile, rfile, logposterfile,
			     sigmapath, lambdapath, mixmomentpath, mweightpath, mlogweightpath, mmeanpath, 
                             Ypath, rpath, logposterpath,
                             n1_censored, 'a');
    if (*doubly){
      openFiles_bayesHistogram(sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, 
                               Yfile2, rfile2, logposterfile2,
	  		       sigmapath2, lambdapath2, mixmomentpath2, mweightpath2, mlogweightpath2, mmeanpath2, 
                               Ypath2, rpath2, logposterpath2,
                               1, 'a');
    }


    /***** Working arrays and variables needed by error terms G-splines                              *****/
    /***** Many of these variables are doubled in the case of doubly censoring                      *****/ 
    /*** iwork[g_y->total_length()] ................. needed by update_Alloc_GS                 ***/
    /*** iwork[g_y->dim()*g_y->total_length()] ...... needed by writeToFiles_bayesHistogram     ***/
    /*** dwork[g_y->total_length()] ................. needed by update_Alloc_GS                 ***/
    /*** dwork[g_y->total_length()] ................. needed by writeToFiles_bayesHistogram     ***/
    /*** dwork[l_moments_work] ...................... needed by writeToFiles_bayesHistogram     ***/
    /***                                                                                        ***/
    /*** mu[][] ..................................... used by update_Alloc_GS to store knots    ***/
    /*** log_poster[0]               = log-likelihood                                  ***/
    /*** log_poster[1, ... 1 OR dim] = penalty                                         ***/
    /*** log_poster[2 OR dim + 1]    = sum (N_k*log(w_k))                              ***/
    /*** a_ipars[2] ................................. used by Gspline::update_alla     ***/
    /***                                                                               ***/
    int l_moments     = *dimP + ((*dimP)*(1+(*dimP)))/2;
    int l_lambda1     = g_y1->equal_lambda()? 1 : (*dimP);
    int l_lambda2     = g_y2->equal_lambda()? 1 : (*dimP);
    int l_iwork1      = g_y1->dim() * g_y1->total_length();
    int l_iwork2      = g_y2->dim() * g_y2->total_length();
    int l_dwork1      = (g_y1->total_length() > l_moments ? g_y1->total_length() : l_moments);
    int l_dwork2      = (g_y2->total_length() > l_moments ? g_y2->total_length() : l_moments);
    int l_log_poster1 = (g_y1->equal_lambda() ? 3 : 2 + g_y1->dim());
    int l_log_poster2 = (g_y2->equal_lambda() ? 3 : 2 + g_y2->dim());
    int i_logpr1      = (g_y1->equal_lambda() ? 2 : 1 + g_y1->dim());
    int i_logpr2      = (g_y2->equal_lambda() ? 2 : 1 + g_y2->dim());
    int a_ipars[2];
    a_ipars[0] = *nP;

    int *iwork1         = (int*) calloc(l_iwork1, sizeof(int));
    double *dwork1      = (double*) calloc(l_dwork1, sizeof(double));
    double **mu1        = (double**) calloc(g_y1->dim(), sizeof(double*));
    double *log_poster1 = (double*) calloc(l_log_poster1, sizeof(double));
    if (!iwork1 || !dwork1 || !mu1 || !log_poster1) throw returnR("Not enough memory available in bayessurvreg2", 1);
    for (j = 0; j < g_y1->dim(); j++){
      mu1[j] = (double*) calloc(g_y1->length(j), sizeof(double));
      if (!mu1[j]) throw returnR("Not enough memory available in bayessurvreg2", 1);
      g_y1->muP(j, mu1[j]);
    }

    int *iwork2 = &itemp;
    double *dwork2 = &dtemp;
    double **mu2 = &ddtemp; 
    double *log_poster2 = &dtemp;
    if (*doubly){
      iwork2      = (int*) calloc(l_iwork2, sizeof(int));
      dwork2      = (double*) calloc(l_dwork2, sizeof(double));
      mu2         = (double**) calloc(g_y2->dim(), sizeof(double*));
      log_poster2 = (double*) calloc(l_log_poster2, sizeof(double));
      if (!iwork2 || !dwork2 || !mu2 || !log_poster2) throw returnR("Not enough memory available in bayessurvreg2", 1);
      for (j = 0; j < g_y2->dim(); j++){
        mu2[j] = (double*) calloc(g_y2->length(j), sizeof(double));
        if (!mu2[j]) throw returnR("Not enough memory available in bayessurvreg2", 1);
        g_y2->muP(j, mu2[j]);
      }
    }

    /***** Working arrays and variables needed by  random intercept G-splines                       *****/
    /***** Many of these variables are doubled in the case of doubly censoring                      *****/ 
    int l_moments_b     = 2;
    int l_lambda1_b     = 1;
    int l_lambda2_b     = 1;
    int l_iwork1_b      = 1 * g_bb1->total_length();
    int l_iwork2_b      = 1 * g_bb2->total_length();
    int l_dwork1_b      = (g_bb1->total_length() > l_moments_b ? g_bb1->total_length() : l_moments_b);
    int l_dwork2_b      = (g_bb2->total_length() > l_moments_b ? g_bb2->total_length() : l_moments_b);
    int l_log_poster1_b = 3;
    int l_log_poster2_b = 3;
    int i_logpr1_b      = 2;
    int i_logpr2_b      = 2;
    int a_ipars1_b[2];
    int a_ipars2_b[2];

    int *iwork1_b = &itemp;
    int *iwork2_b = &itemp;
    double *dwork1_b = &dtemp;
    double *log_poster1_b = &dtemp;
    double *dwork2_b = &dtemp;
    double *log_poster2_b = &dtemp;
    double **mu1_b = &ddtemp;
    double **mu2_b = &ddtemp;

    switch (*version){
    case 2:
    case 32:
      break;

    case 3:
    case 31:
      if (bb1->nRandom()){
        iwork1_b      = (int*) calloc(l_iwork1_b, sizeof(int));
        dwork1_b      = (double*) calloc(l_dwork1_b, sizeof(double));
        mu1_b         = (double**) calloc(g_bb1->dim(), sizeof(double*));
        log_poster1_b = (double*) calloc(l_log_poster1_b, sizeof(double));
        if (!iwork1_b || !dwork1_b || !mu1_b || !log_poster1_b) throw returnR("Not enough memory available in bayessurvreg2", 1);
        for (j = 0; j < g_bb1->dim(); j++){
          mu1_b[j] = (double*) calloc(g_bb1->length(j), sizeof(double));
          if (!mu1_b[j]) throw returnR("Not enough memory available in bayessurvreg2", 1);
          g_bb1->muP(j, mu1_b[j]);
        }
        a_ipars1_b[0] = bb1->nCluster();
      }

      if (bb2->nRandom()){
        iwork2_b      = (int*) calloc(l_iwork2_b, sizeof(int));
        dwork2_b      = (double*) calloc(l_dwork2_b, sizeof(double));
        mu2_b         = (double**) calloc(g_bb2->dim(), sizeof(double*));
        log_poster2_b = (double*) calloc(l_log_poster2_b, sizeof(double));
        if (!iwork2_b || !dwork2_b || !mu2_b || !log_poster2_b) throw returnR("Not enough memory available in bayessurvreg2", 1);
        for (j = 0; j < g_bb2->dim(); j++){
          mu2_b[j] = (double*) calloc(g_bb2->length(j), sizeof(double));
          if (!mu2_b[j]) throw returnR("Not enough memory available in bayessurvreg2", 1);
          g_bb2->muP(j, mu2_b[j]);
        }
      }
      a_ipars2_b[0] = bb2->nCluster();
      break;
    }

    //Rprintf("\nERROR ONSET: initial G-spline:");
    //g_y1->print();
    //Rprintf("\nRANDOM ONSET: initial G-spline:");
    //g_bb1->print();
    //Rprintf("\nERROR EVENT: initial G-spline:");
    //g_y2->print();
    //Rprintf("\nRANDOM EVENT: initial G-spline:");
    //g_bb2->print();

    /** Main simulation: **/
      /* nstored ........ number of values currently stored in _work arrays                */
      /* write_flag ..... flag for writing to files 'a' = append                           */
      /* nullthIter ..... index for initial values (usually 0)                             */
      /* iter ........... counter of iterations                                            */
      /* iterTotal ...... total number of iterations done previously (thinnig included)    */
      /* iterTotalNow ... total number of iterations done with this function call          */
      /*                  -> this number would be used to compute acceptance rates         */  
      /* backs .......... how many times the carriage must be returned to print            */
      /*                  the next iteration number?                                       */
    int nullthIter   = *iterM;
    int iter;
    int lastIter     = nullthIter + niter;
    int iterTotal    = nullthIter*nthin;
    int iterTotalNow = 0;
    int backs        = 0;
    int writeAll     = 0;

    int *rho_acceptP = rho_accept;

    Rprintf("Iteration "); 
    switch (*version){

    /******* VERSION 2 ********/
    /******* ========= ********/    
    case 2:
      for (iter=nullthIter + 1; iter <= lastIter; iter++){
        for (int witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          update_Data_GS_regres(Y1, regresRes1, y1_left, y1_right, status1, r1, g_y1, nP);
          update_Alloc_GS(r1, mixtureN1, mu1, log_poster1 + 0, log_poster1 + i_logpr1, g_y1, regresRes1, nP, iwork1, dwork1);
          g_y1->update_alla_lambda(mixtureN1, a_ipars, &iterTotalNow);

          switch (specif[0]){
          case 1:
            g_y1->update_gamma(regresRes1, r1, nP);
            g_y1->update_sigma(regresRes1, r1, nP, &iterTotalNow);
            break;
          case 2:
            g_y1->update_Intcpt(regresRes1, r1, nP);
            g_y1->update_Scale(regresRes1, r1, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y1->dim(); j++) g_y1->muP(j, mu1[j]);

          beta1->GIBBSfixed(regresRes1, nP, X1, XXtb1, g_y1, mu1, r1);
          bb1->GIBBSupdate(regresRes1, nP, X1, ZZtb1, g_y1, mu1, r1, beta1, DD1);
          DD1->GIBBSnormalRE(bb1, beta1);
          beta1->GIBBSmeanRandom(bb1, DD1);

          if (*doubly){
            update_Data_GS_doubly(Y2, regresRes2, Y1, t2_left, t2_right, status2, r2, g_y2, nP);
            update_Alloc_GS(r2, mixtureN2, mu2, log_poster2 + 0, log_poster2 + i_logpr2, g_y2, regresRes2, nP, iwork2, dwork2);
            g_y2->update_alla_lambda(mixtureN2, a_ipars, &iterTotalNow);
            switch (specif[1]){
            case 1:
              g_y2->update_gamma(regresRes2, r2, nP);
              g_y2->update_sigma(regresRes2, r2, nP, &iterTotalNow);
              break;
            case 2:
              g_y2->update_Intcpt(regresRes2, r2, nP);
              g_y2->update_Scale(regresRes2, r2, nP, &iterTotalNow);
              break;
            }
            for (j = 0; j < g_y2->dim(); j++) g_y2->muP(j, mu2[j]);

            beta2->GIBBSfixed(regresRes2, nP, X2, XXtb2, g_y2, mu2, r2);

            bb2->GIBBSupdate(regresRes2, nP, X2, ZZtb2, g_y2, mu2, r2, beta2, DD2);
            DD2->GIBBSnormalRE(bb2, beta2);
            beta2->GIBBSmeanRandom(bb2, DD2);
          }
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
        print_iter_info(writeAll, backs, iter, nwrite, lastIter);

        if (writeSumExpa){
          writeToFile_1(g_y1->sumexpaP(), 1, sumexpafile, out_format[0], out_format[1]);
          writeToFile_1(g_y1->a_maxP(), 1, maxafile, out_format[0], out_format[1]);

          if (*doubly){
            writeToFile_1(g_y2->sumexpaP(), 1, sumexpa2file, out_format[0], out_format[1]);
            writeToFile_1(g_y2->a_maxP(), 1, maxa2file, out_format[0], out_format[1]);
          }
        }

        if ((*mainSimul) || writeAll){
          writeToFile_1(&iter, 1, iterfile, out_format[0], out_format[1]);

          if (beta1->nbeta()) writeToFile_1(beta1->betaP(), beta1->nbeta(), betafile, out_format[0], out_format[1]);
          if (beta2->nbeta()) writeToFile_1(beta2->betaP(), beta2->nbeta(), betafile2, out_format[0], out_format[1]);

          if (bb1->nRandom()) writeToFiles_random(DD1, bb1, storeb1P, &writeAll, Dfile, bbfile, out_format[0], out_format[1]);
          if (bb2->nRandom()) writeToFiles_random(DD2, bb2, storeb2P, &writeAll, Dfile2, bbfile2, out_format[0], out_format[1]);
 
          writeToFiles_bayesHistogram(g_y1, r1, Y1, log_poster1,
	      l_moments, l_lambda1, l_log_poster1, nP, storea1P, storey1P, storer1P, n1_censored, &writeAll, iwork1, dwork1,
              sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
   	      _null_weight1, out_format[0], out_format[1], check_k_effect);
          if (*doubly){
            writeToFiles_bayesHistogram(g_y2, r2, Y2, log_poster2,
	        l_moments, l_lambda2, l_log_poster2, nP, storea2P, storey2P, storer2P, 1, &writeAll, iwork2, dwork2,
                sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2,
	        _null_weight2, out_format[0], out_format[1], check_k_effect);
          }
          writeAll = 0;
        }
      }    /** end of the main cycle over iter **/
      break;


    /******* VERSION 3 ********/
    /******* ========= ********/    
    case 3:
      for (iter=nullthIter + 1; iter <= lastIter; iter++){
        for (int witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          update_Data_GS_regres(Y1, regresRes1, y1_left, y1_right, status1, r1, g_y1, nP);
          update_Alloc_GS(r1, mixtureN1, mu1, log_poster1 + 0, log_poster1 + i_logpr1, g_y1, regresRes1, nP, iwork1, dwork1);
          g_y1->update_alla_lambda(mixtureN1, a_ipars, &iterTotalNow);

          switch (specif[0]){
          case 1:
            g_y1->update_gamma(regresRes1, r1, nP);
            g_y1->update_sigma(regresRes1, r1, nP, &iterTotalNow);
            break;
          case 2:
            g_y1->update_Intcpt(regresRes1, r1, nP);
            g_y1->update_Scale(regresRes1, r1, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y1->dim(); j++) g_y1->muP(j, mu1[j]);

          beta1->GIBBSfixed(regresRes1, nP, X1, XXtb1, g_y1, mu1, r1);

          bb1->Gspl_intcpt_update(regresRes1, nP, g_bb1, mu1_b, r_b1, g_y1, mu1, r1);
          if (bb1->nRandom()){
            update_Alloc_GS(r_b1, mixtureN_b1, mu1_b, log_poster1_b + 0, log_poster1_b + i_logpr1_b, g_bb1, bb1->bMP(), 
                            &nCluster1, iwork1_b, dwork1_b);
            g_bb1->update_alla_lambda(mixtureN_b1, a_ipars1_b, &iterTotalNow);
            switch (specif_b[0]){
            case 1:
              g_bb1->update_sigma(bb1->bMP(), r_b1, &nCluster1, &iterTotalNow);
              break;
            case 2:
              g_bb1->update_Scale(bb1->bMP(), r_b1, &nCluster1, &iterTotalNow);
              break;
            }
            for (j = 0; j < g_bb1->dim(); j++) g_bb1->muP(j, mu1_b[j]);
            /* adjust_intcpt(g_y1, g_bb1, bb1); */
          }

          if (*doubly){
            update_Data_GS_doubly(Y2, regresRes2, Y1, t2_left, t2_right, status2, r2, g_y2, nP);
            update_Alloc_GS(r2, mixtureN2, mu2, log_poster2 + 0, log_poster2 + i_logpr2, g_y2, regresRes2, nP, iwork2, dwork2);
            g_y2->update_alla_lambda(mixtureN2, a_ipars, &iterTotalNow);
            switch (specif[1]){
            case 1:
              g_y2->update_gamma(regresRes2, r2, nP);
              g_y2->update_sigma(regresRes2, r2, nP, &iterTotalNow);
              break;
            case 2:
              g_y2->update_Intcpt(regresRes2, r2, nP);
              g_y2->update_Scale(regresRes2, r2, nP, &iterTotalNow);
              break;
            }
            for (j = 0; j < g_y2->dim(); j++) g_y2->muP(j, mu2[j]);

            beta2->GIBBSfixed(regresRes2, nP, X2, XXtb2, g_y2, mu2, r2);

            bb2->Gspl_intcpt_update(regresRes2, nP, g_bb2, mu2_b, r_b2, g_y2, mu2, r2);

            if (bb2->nRandom()){
              update_Alloc_GS(r_b2, mixtureN_b2, mu2_b, log_poster2_b + 0, log_poster2_b + i_logpr2_b, g_bb2, bb2->bMP(), 
                              &nCluster2, iwork2_b, dwork2_b);
              g_bb2->update_alla_lambda(mixtureN_b2, a_ipars2_b, &iterTotalNow);
              switch (specif_b[1]){
              case 1:
                g_bb2->update_sigma(bb2->bMP(), r_b2, &nCluster2, &iterTotalNow);
                break;
              case 2:
                g_bb2->update_Scale(bb2->bMP(), r_b2, &nCluster2, &iterTotalNow);
                break;
              }
              for (j = 0; j < g_bb2->dim(); j++) g_bb2->muP(j, mu2_b[j]);
              /* adjust_intcpt(g_y2, g_bb2, bb2); */
            }     
          }
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
        print_iter_info(writeAll, backs, iter, nwrite, lastIter);

        if (writeSumExpa){
          writeToFile_1(g_y1->sumexpaP(), 1, sumexpafile, out_format[0], out_format[1]);
          writeToFile_1(g_y1->a_maxP(), 1, maxafile, out_format[0], out_format[1]);

          writeToFile_1(g_bb1->sumexpaP(), 1, sumexpa_bfile, out_format[0], out_format[1]);
          writeToFile_1(g_bb1->a_maxP(), 1, maxa_bfile, out_format[0], out_format[1]);

          if (*doubly){
            writeToFile_1(g_y2->sumexpaP(), 1, sumexpa2file, out_format[0], out_format[1]);
            writeToFile_1(g_y2->a_maxP(), 1, maxa2file, out_format[0], out_format[1]);

            writeToFile_1(g_bb2->sumexpaP(), 1, sumexpa_b2file, out_format[0], out_format[1]);
            writeToFile_1(g_bb2->a_maxP(), 1, maxa_b2file, out_format[0], out_format[1]);
          }
        }

        if ((*mainSimul) || writeAll){
          writeToFile_1(&iter, 1, iterfile, out_format[0], out_format[1]);

          if (beta1->nbeta()) writeToFile_1(beta1->betaP(), beta1->nbeta(), betafile, out_format[0], out_format[1]);
          if (beta2->nbeta()) writeToFile_1(beta2->betaP(), beta2->nbeta(), betafile2, out_format[0], out_format[1]);

          if (bb1->nRandom()){
            writeToFiles_Gspl_intcpt(bb1, storeb1P, &writeAll, bbfile, out_format[0], out_format[1]);
            writeToFiles_bayesHistogram(g_bb1, r_b1, &ZERO, log_poster1_b,
                l_moments_b, l_lambda1_b, l_log_poster1_b, &nCluster1, &storea_b1P, &ZERO_INT, &storer_b1P, 0, &writeAll, iwork1_b, dwork1_b,
                sigmafile_b, lambdafile_b, mixmomentfile_b, mweightfile_b, mlogweightfile_b, mmeanfile_b, Yfile, rfile_b, logposterfile_b,
  	        _null_weight1_b, out_format[0], out_format[1], check_k_effect);
          }
          if (bb2->nRandom()){
            writeToFiles_Gspl_intcpt(bb2, storeb2P, &writeAll, bbfile2, out_format[0], out_format[1]);
            writeToFiles_bayesHistogram(g_bb2, r_b2, &ZERO, log_poster2_b,
	        l_moments_b, l_lambda2_b, l_log_poster2_b, &nCluster2, &storea_b2P, &ZERO_INT, &storer_b2P, 0, &writeAll, iwork2_b, dwork2_b,
                sigmafile_b2, lambdafile_b2, mixmomentfile_b2, mweightfile_b2, mlogweightfile_b2, mmeanfile_b2, Yfile, rfile_b2, logposterfile_b2,
	        _null_weight2_b, out_format[0], out_format[1], check_k_effect);
          }
 
          writeToFiles_bayesHistogram(g_y1, r1, Y1, log_poster1,
	      l_moments, l_lambda1, l_log_poster1, nP, storea1P, storey1P, storer1P, n1_censored, &writeAll, iwork1, dwork1,
              sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
  	    _null_weight1, out_format[0], out_format[1], check_k_effect);
          if (*doubly){
            writeToFiles_bayesHistogram(g_y2, r2, Y2, log_poster2,
	        l_moments, l_lambda2, l_log_poster2, nP, storea2P, storey2P, storer2P, 1, &writeAll, iwork2, dwork2,
                sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2,
	        _null_weight2, out_format[0], out_format[1], check_k_effect);
          }
          writeAll = 0;
        }
      }    /** end of the main cycle over iter **/
      break;


    /******* VERSION 31 ********/
    /******* ========== ********/    
    case 31:
      for (iter=nullthIter + 1; iter <= lastIter; iter++){
        for (int witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          update_Data_GS_regres(Y1, regresRes1, y1_left, y1_right, status1, r1, g_y1, nP);
          update_Alloc_GS(r1, mixtureN1, mu1, log_poster1 + 0, log_poster1 + i_logpr1, g_y1, regresRes1, nP, iwork1, dwork1);
          g_y1->update_alla_lambda(mixtureN1, a_ipars, &iterTotalNow);

          switch (specif[0]){
          case 1:
            g_y1->update_gamma(regresRes1, r1, nP);
            g_y1->update_sigma(regresRes1, r1, nP, &iterTotalNow);
            break;
          case 2:
            g_y1->update_Intcpt(regresRes1, r1, nP);
            g_y1->update_Scale(regresRes1, r1, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y1->dim(); j++) g_y1->muP(j, mu1[j]);

          beta1->GIBBSfixed(regresRes1, nP, X1, XXtb1, g_y1, mu1, r1);

          Gspl_rho_intcpt_update(bb1, bb2, rho_zb, regresRes1, regresRes2, rho_acceptP, nP, rhobI, rhobD,
                                 g_bb1, mu1_b, r_b1, g_bb2, mu2_b, r_b2, g_y1, mu1, r1, g_y2, mu2, r2);

          update_Alloc_GS(r_b1, mixtureN_b1, mu1_b, log_poster1_b + 0, log_poster1_b + i_logpr1_b, g_bb1, bb1->bMP(), 
                          &nCluster1, iwork1_b, dwork1_b);
          g_bb1->update_alla_lambda(mixtureN_b1, a_ipars1_b, &iterTotalNow);
          switch (specif_b[0]){
          case 1:
            g_bb1->update_sigma(bb1->bMP(), r_b1, &nCluster1, &iterTotalNow);
            break;
          case 2:
            g_bb1->update_Scale(bb1->bMP(), r_b1, &nCluster1, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_bb1->dim(); j++) g_bb1->muP(j, mu1_b[j]);
          /* adjust_intcpt(g_y1, g_bb1, bb1); */

          update_Data_GS_doubly(Y2, regresRes2, Y1, t2_left, t2_right, status2, r2, g_y2, nP);
          update_Alloc_GS(r2, mixtureN2, mu2, log_poster2 + 0, log_poster2 + i_logpr2, g_y2, regresRes2, nP, iwork2, dwork2);
          g_y2->update_alla_lambda(mixtureN2, a_ipars, &iterTotalNow);
          switch (specif[1]){
          case 1:
            g_y2->update_gamma(regresRes2, r2, nP);
            g_y2->update_sigma(regresRes2, r2, nP, &iterTotalNow);
            break;
          case 2:
            g_y2->update_Intcpt(regresRes2, r2, nP);
            g_y2->update_Scale(regresRes2, r2, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y2->dim(); j++) g_y2->muP(j, mu2[j]);

          beta2->GIBBSfixed(regresRes2, nP, X2, XXtb2, g_y2, mu2, r2);

          update_Alloc_GS(r_b2, mixtureN_b2, mu2_b, log_poster2_b + 0, log_poster2_b + i_logpr2_b, g_bb2, bb2->bMP(), 
                          &nCluster2, iwork2_b, dwork2_b);
          g_bb2->update_alla_lambda(mixtureN_b2, a_ipars2_b, &iterTotalNow);
          switch (specif_b[1]){
          case 1:
            g_bb2->update_sigma(bb2->bMP(), r_b2, &nCluster2, &iterTotalNow);
            break;
          case 2:
            g_bb2->update_Scale(bb2->bMP(), r_b2, &nCluster2, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_bb2->dim(); j++) g_bb2->muP(j, mu2_b[j]);
          /* adjust_intcpt(g_y2, g_bb2, bb2); */
        }    /** end of the thinning cycle  **/
        rho_acceptP++;

        /**  Write to files  **/

        print_iter_info(writeAll, backs, iter, nwrite, lastIter);

        if (writeSumExpa){
          writeToFile_1(g_y1->sumexpaP(), 1, sumexpafile, out_format[0], out_format[1]);
          writeToFile_1(g_y1->a_maxP(), 1, maxafile, out_format[0], out_format[1]);

          writeToFile_1(g_bb1->sumexpaP(), 1, sumexpa_bfile, out_format[0], out_format[1]);
          writeToFile_1(g_bb1->a_maxP(), 1, maxa_bfile, out_format[0], out_format[1]);

          writeToFile_1(g_y2->sumexpaP(), 1, sumexpa2file, out_format[0], out_format[1]);
          writeToFile_1(g_y2->a_maxP(), 1, maxa2file, out_format[0], out_format[1]);

          writeToFile_1(g_bb2->sumexpaP(), 1, sumexpa_b2file, out_format[0], out_format[1]);
          writeToFile_1(g_bb2->a_maxP(), 1, maxa_b2file, out_format[0], out_format[1]);
        }

        if ((*mainSimul) || writeAll){
          writeToFile_1(&iter, 1, iterfile, out_format[0], out_format[1]);

          if (beta1->nbeta()) writeToFile_1(beta1->betaP(), beta1->nbeta(), betafile, out_format[0], out_format[1]);
          if (beta2->nbeta()) writeToFile_1(beta2->betaP(), beta2->nbeta(), betafile2, out_format[0], out_format[1]);

          writeToFiles_Gspl_intcpt(bb1, storeb1P, &writeAll, bbfile, out_format[0], out_format[1]);
          writeToFiles_bayesHistogram(g_bb1, r_b1, &ZERO, log_poster1_b,
              l_moments_b, l_lambda1_b, l_log_poster1_b, &nCluster1, &storea_b1P, &ZERO_INT, &storer_b1P, 0, &writeAll, iwork1_b, dwork1_b,
              sigmafile_b, lambdafile_b, mixmomentfile_b, mweightfile_b, mlogweightfile_b, mmeanfile_b, Yfile, rfile_b, logposterfile_b,
  	      _null_weight1_b, out_format[0], out_format[1], check_k_effect);

          writeToFiles_Gspl_intcpt(bb2, storeb2P, &writeAll, bbfile2, out_format[0], out_format[1]);
          writeToFiles_bayesHistogram(g_bb2, r_b2, &ZERO, log_poster2_b,
	      l_moments_b, l_lambda2_b, l_log_poster2_b, &nCluster2, &storea_b2P, &ZERO_INT, &storer_b2P, 0, &writeAll, iwork2_b, dwork2_b,
              sigmafile_b2, lambdafile_b2, mixmomentfile_b2, mweightfile_b2, mlogweightfile_b2, mmeanfile_b2, Yfile, rfile_b2, logposterfile_b2,
	      _null_weight2_b, out_format[0], out_format[1], check_k_effect);

          writeToFile_1(rho_zb, 2, rhobfile, out_format[0],out_format[1]);
 
          writeToFiles_bayesHistogram(g_y1, r1, Y1, log_poster1,
	      l_moments, l_lambda1, l_log_poster1, nP, storea1P, storey1P, storer1P, n1_censored, &writeAll, iwork1, dwork1,
              sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
  	      _null_weight1, out_format[0], out_format[1], check_k_effect);

          writeToFiles_bayesHistogram(g_y2, r2, Y2, log_poster2,
	      l_moments, l_lambda2, l_log_poster2, nP, storea2P, storey2P, storer2P, 1, &writeAll, iwork2, dwork2,
              sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2,
	      _null_weight2, out_format[0], out_format[1], check_k_effect);

          writeAll = 0;
        }
      }    /** end of the main cycle over iter **/
      break;


    /******* VERSION 32 ********/
    /******* ========== ********/    
    case 32:
      for (iter=nullthIter + 1; iter <= lastIter; iter++){
        for (int witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
          iterTotal++;                        /* = (iter-1)*nthin + witer                    */
          iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

          update_Data_GS_regres(Y1, regresRes1, y1_left, y1_right, status1, r1, g_y1, nP);
          update_Alloc_GS(r1, mixtureN1, mu1, log_poster1 + 0, log_poster1 + i_logpr1, g_y1, regresRes1, nP, iwork1, dwork1);
          g_y1->update_alla_lambda(mixtureN1, a_ipars, &iterTotalNow);

          switch (specif[0]){
          case 1:
            g_y1->update_gamma(regresRes1, r1, nP);
            g_y1->update_sigma(regresRes1, r1, nP, &iterTotalNow);
            break;
          case 2:
            g_y1->update_Intcpt(regresRes1, r1, nP);
            g_y1->update_Scale(regresRes1, r1, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y1->dim(); j++) g_y1->muP(j, mu1[j]);

          beta1->GIBBSfixed(regresRes1, nP, X1, XXtb1, g_y1, mu1, r1);

          RandomEff32::update(db, regresRes1, regresRes2, nP, g_y1, mu1, r1, g_y2, mu2, r2);

          update_Data_GS_doubly(Y2, regresRes2, Y1, t2_left, t2_right, status2, r2, g_y2, nP);
          update_Alloc_GS(r2, mixtureN2, mu2, log_poster2 + 0, log_poster2 + i_logpr2, g_y2, regresRes2, nP, iwork2, dwork2);
          g_y2->update_alla_lambda(mixtureN2, a_ipars, &iterTotalNow);
          switch (specif[1]){
          case 1:
            g_y2->update_gamma(regresRes2, r2, nP);
            g_y2->update_sigma(regresRes2, r2, nP, &iterTotalNow);
            break;
          case 2:
            g_y2->update_Intcpt(regresRes2, r2, nP);
            g_y2->update_Scale(regresRes2, r2, nP, &iterTotalNow);
            break;
          }
          for (j = 0; j < g_y2->dim(); j++) g_y2->muP(j, mu2[j]);

          beta2->GIBBSfixed(regresRes2, nP, X2, XXtb2, g_y2, mu2, r2);           
        }    /** end of the thinning cycle  **/

        /**  Write to files  **/
        print_iter_info(writeAll, backs, iter, nwrite, lastIter);

        if (writeSumExpa){
          writeToFile_1(g_y1->sumexpaP(), 1, sumexpafile, out_format[0], out_format[1]);
          writeToFile_1(g_y1->a_maxP(), 1, maxafile, out_format[0], out_format[1]);

          writeToFile_1(g_y2->sumexpaP(), 1, sumexpa2file, out_format[0], out_format[1]);
          writeToFile_1(g_y2->a_maxP(), 1, maxa2file, out_format[0], out_format[1]);
        }

        if ((*mainSimul) || writeAll){
          writeToFile_1(&iter, 1, iterfile, out_format[0], out_format[1]);

          if (beta1->nbeta()) writeToFile_1(beta1->betaP(), beta1->nbeta(), betafile, out_format[0], out_format[1]);
          if (beta2->nbeta()) writeToFile_1(beta2->betaP(), beta2->nbeta(), betafile2, out_format[0], out_format[1]);

          writeToFiles_RandomEff32(db, storeb1P, storeb2P, &writeAll, Dfile, bbfile, bbfile2, out_format[0], out_format[1]);
 
          writeToFiles_bayesHistogram(g_y1, r1, Y1, log_poster1,
	      l_moments, l_lambda1, l_log_poster1, nP, storea1P, storey1P, storer1P, n1_censored, &writeAll, iwork1, dwork1,
              sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
   	      _null_weight1, out_format[0], out_format[1], check_k_effect);

          writeToFiles_bayesHistogram(g_y2, r2, Y2, log_poster2,
	      l_moments, l_lambda2, l_log_poster2, nP, storea2P, storey2P, storer2P, 1, &writeAll, iwork2, dwork2,
              sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2,
	      _null_weight2, out_format[0], out_format[1], check_k_effect);

          writeAll = 0;
        }
      }    /** end of the main cycle over iter **/
      break;

    default:
      throw returnR("Error in bayessurvreg2.cpp. Main interation loop not implemented for this version", 1);
    }

    //Rprintf("\nERROR ONSET: last G-spline:");
    //g_y1->print();
    //Rprintf("\nRANDOM ONSET: last G-spline:");
    //g_bb1->print();
    //Rprintf("\nERROR EVENT: last G-spline:");
    //g_y2->print();
    //Rprintf("\nRANDOM EVENT: last G-spline:");
    //g_bb2->print();


    if (writeSumExpa){
      sumexpafile.close();
      maxafile.close();

      if (*version == 3 || *version == 31){
        sumexpa_bfile.close();
        maxa_bfile.close();
      }

      if (*doubly){
        sumexpa2file.close();
        maxa2file.close();

        if (*version == 3 || *version == 31){
          sumexpa_b2file.close();
          maxa_b2file.close();
        }
      }
    }

    iterfile.close();
    if (beta1->nbeta()) betafile.close();
    if (beta2->nbeta()) betafile2.close();

    switch (*version){
    case 2:
      if (bb1->nRandom()){
        bbfile.close();
        Dfile.close();
      }
      if (bb2->nRandom()){
        bbfile2.close();
        Dfile2.close();
      }
      break;

    case 31:
      rhobfile.close();

    case 3:
      if (bb1->nRandom()){
        bbfile.close();
        closeFiles_bayesHistogram(sigmafile_b, lambdafile_b, mixmomentfile_b, mweightfile_b, 
                                  mlogweightfile_b, mmeanfile_b, Yfile, rfile_b, logposterfile_b, 0);
      }
      if (bb2->nRandom()){
        bbfile2.close();
        closeFiles_bayesHistogram(sigmafile_b2, lambdafile_b2, mixmomentfile_b2, mweightfile_b2, 
                                  mlogweightfile_b2, mmeanfile_b2, Yfile, rfile_b2, logposterfile_b2, 0);
      }
      break;

    case 32:
      bbfile.close();
      bbfile2.close();
      Dfile.close();
      break;
    }

    closeFiles_bayesHistogram(sigmafile, lambdafile, mixmomentfile, mweightfile, 
                              mlogweightfile, mmeanfile, Yfile, rfile, logposterfile, n1_censored);
    if (*doubly){
      closeFiles_bayesHistogram(sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, 
                                mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2, 1);
    }
    Rprintf("\n");


    /**  Do some adjustments on last sampled values before returning them back to R  **/
    *iterM = iter - 1;
    for (i = 0; i < *nP; i++) r1[i]++;                     // C++ -> R
    if (*doubly) for (i = 0; i < *nP; i++) r2[i]++;        // C++ -> R

    /** Write current G-spline values to initial arrays **/
    g_y1->Gspline2initArray(GsplineI1, GsplineD1);
    if (*doubly) g_y2->Gspline2initArray(GsplineI2, GsplineD2);

    /** Write current beta/gamma values to initial arrays **/
    if (beta1->nbeta()) beta1->BetaGamma2initArray(priorBeta1I, priorBeta1D);
    if (beta2->nbeta()) beta2->BetaGamma2initArray(priorBeta2I, priorBeta2D);

    /** Write current covariance matrix of random effects and random effects/random intercept G-spline  to initial arrays (if needed) **/
    switch (*version){
    case 2:
      if (bb1->nRandom()){
        bb1->RandomEff2initArray(priorb1I, priorb1D);
        DD1->CovMatrix2initArray(priorCovMat1I, priorCovMat1D);
      }
      if (bb2->nRandom()){
        bb2->RandomEff2initArray(priorb2I, priorb2D);
        DD2->CovMatrix2initArray(priorCovMat2I, priorCovMat2D);
      }
      break;

    case 31:
      *rhob = rho_zb[0];
    case 3:
      if (bb1->nRandom()){
        bb1->RandomEff2initArray(priorb1I, priorb1D);
        for (i = 0; i < nCluster1; i++) r_b1[i]++;                     // C++ -> R
        g_bb1->Gspline2initArray(priorCovMat1I, priorCovMat1D);
      }
      if (bb2->nRandom()){
        bb2->RandomEff2initArray(priorb2I, priorb2D);
        for (i = 0; i < nCluster2; i++) r_b2[i]++;                     // C++ -> R
        g_bb2->Gspline2initArray(priorCovMat2I, priorCovMat2D);
      }
      break;

    case 32:
      break;
    }


    /*** Cleaning ***/
    switch (*version){
    case 2:
      break;

    case 3:
    case 31:
      if (bb1->nRandom()){
        for (j = 0; j < g_bb1->dim(); j++) free(mu1_b[j]);
        free(iwork1_b);    
        free(dwork1_b);    
        free(mu1_b);   
        free(log_poster1_b);
        free(mixtureN_b1);
      }
      if (bb2->nRandom()){
        for (j = 0; j < g_bb2->dim(); j++) free(mu2_b[j]);
        free(iwork2_b);    
        free(dwork2_b);    
        free(mu2_b);   
        free(log_poster2_b);
        free(mixtureN_b2);
      }
      break;
    }


    for (j = 0; j < g_y1->dim(); j++) free(mu1[j]);
    free(iwork1);    
    free(dwork1);    
    free(mu1);   
    free(log_poster1);
    if (*doubly){
      for (j = 0; j < g_y2->dim(); j++) free(mu2[j]);
      free(iwork2);    
      free(dwork2);    
      free(mu2);   
      free(log_poster2);
    }

    free(mixtureN1);                       
    if (*doubly) free(mixtureN2);

    delete g_y1;                           
    delete g_y2;

    switch (*version){
    case 2:
      if (beta1->nRandom()) free(ZZtb1);
      if (*doubly && beta2->nRandom()) free(ZZtb2);
      break;

    case 3:
    case 31:
      break;
    }
    if (beta1->nbeta()) free(XXtb1);
    if (*doubly && beta2->nbeta()) free(XXtb2);

    free(regresRes1);
    if (*doubly) free(regresRes2);

    delete db;

    delete g_bb1;
    delete g_bb2;

    delete DD2;
    delete DD1;

    delete bb1;
    delete bb2;

    delete beta1;
    delete beta2;

    free(ddtemp);
    PutRNGstate();
    return;
  }
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }
 }  /** end of the function bayessurvreg2 **/


}  /** end of extern "C" **/


//
/***** adjust_intcpt                          *****/
//
// * used in the case of G-spline random intercept model
//
// g_eps ....... G-spline for the random error term
// g_b ......... G-spline for the random intercept
// bb .......... individual random intercepts
//
void
adjust_intcpt(Gspline* g_eps,  Gspline* g_b,  RandomEff* bb)
{
  static double Eb;

  /*** 1) compute expectation of the random intercept G-spline  ***/
  g_b->mean_univariate(&Eb);

  /*** 2) subtract Eb from all individual random intercepts     ***/
  bb->adjust_intcpt(&Eb);

  /*** 3) add Eb to error (epsilon) G-spline intercept  ***/
  g_eps->adjust_intcpt(&Eb);

  /*** 4) subtract Eb from random intercept G-spline intercept  ***/
  Eb *= (-1);
  g_b->adjust_intcpt(&Eb);

  return;
}


void
print_iter_info(int &writeAll,  int &backs,  const int &iter,  const int &nwrite,  const int &lastIter)
{
  static int i;

  if (!(iter % nwrite) || iter == lastIter){
    writeAll = 1;
    for (i = 0; i < backs; i++) Rprintf("\b");
    Rprintf("%d", iter);
    backs = int(log10(double(iter))) + 1;
  }

  return;
}
