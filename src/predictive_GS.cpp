// Function to compute predictive survivor, hazard, cumulative hazard curves and densities
//   for specific values of covariates and models with G-splines
//   * in the case of a bivariate model, marginal functions are always computed
//   * possibly computed also quantiles of these quantities
//
//   

//
// 22/01/2005: start working on it
// 24/01/2005: 'predictive_GS'
//             'evalPredFuns'
// 31/01/2005: possibility for normal random effects added to 'predictive_GS'
// 07/02/2005: possibility for G-spline random intercept added to 'predictive_GS'
// 20/04/2005: bug in 'evalPredFuns' causing SegFault in some cases fixed
// 13/12/2006: version 32 implemented (prediction in the model with doubly censored data and bivariate normal
//             random intercepts)
//

#include "predictive_GS.h"

extern "C"{

/*** predictive_GS                         ***/
/*** ====================================  ***/
//

// 
// PARAMETERS:
//
// averDens[sum(ngrid)] ......... computed predictive density for each observation (no memory has to be supplied 
//                                if predictive density is not asked)
// averS[sum(ngrid)] ............ computed predictive survivor function for each observation (no memory has to be supplied 
//                                if predictive survivor function is not asked)
// averHaz[sum(ngrid)] .......... computed predictive hazard function for each observation (no memory has to be supplied 
//                                if predictive hazard function is not asked)
// averCumHaz[sum(ngrid)] ....... computed predictive hazard function for each observation (no memory has to be supplied 
//                                if predictive cumulative hazard function is not asked)
//
// valDens[M_now*sum(ngrid)] .... values of predictive density at each iteration
//           * only an array of length sum(ngrid) has to be supplied if onlyAver == 1
//           * if onlyAver == 0, it is sorted in the following way:
//             0, ..., M_now_max-1             = values for observation 0, grid[0]
//             M_now_max, ..., 2*(M_now_max)-1 = values for observation 0, grid[1]
//             etc.
//
// valS[M_now*sum(ngrid)] ....... values of predictive density at each iteration
// valHaz[M_now*sum(ngrid)] ..... values of hazards at each iteration
// valCumHaz[M_now*sum(ngrid)] .. values of cumulative hazards at each iteration
//
// quantDens[nquant*sum(ngrid)] .... quantiles of predictive density at each iteration
//            * it does not have to be supplied if onlyAver == 1
//            * sorted in the following way:
//              0, ..., ngrid[0]-1           = quantile[0] for observation 0
//              ngrid[0], ... 2*ngrid[0] - 1 = quantile[1] for observation 0
//              etc.
// quantS[nquant*sum(ngrid)] ....... quantiles of predictive density at each iteration
// quantHaz[nquant*sum(ngrid)] ..... quantiles of hazards at each iteration
// quantCumHaz[nquant*sum(ngrid)] .. quantiles of cumulative hazards at each iteration
//
// dimsP[ ] ..................... dimensionality parameters
// X[nobs*length(beta)] ......... covariates (given in ROW major order, i.e. first all covariates for the first subject etc.)
// obsdims[nobs] ................ indicators of a dimension for each observation
//                                * if dim == 1, this should be a vector filled with zeros only
//                                * if dim == 2, this should be a vector of 0 and 1
// M_now ........................ 
//         on INPUT:  maximal sample size used here
//         on OUTPUT: current sample size used to get average and quantiles here (after taking into account 'skip' and 'by')
//
// dirP ......................... directory where the sample is stored
// extensP ...................... additional extension by file names (usually "_2" for doubly censored data)
// extens_adjP .................. additional extension by file name where to take G-spline for the random intercept
//                                (if version == 3), usually "_b" or "_b2"

// GsplI[ ] ..................... needed G-spline parameters
//   GsplI[0]          = dim (1 or 2)
//   GsplI[1]          = total_length of the G-spline
//   GsplI[2, 2+dim-1] = [K1, K2]
   
// objBetaI[ ] .................. integer parameters for betaGamma constructor
// objBetaD[ ] .................. double parameters for betaGamma constructor

// objbI[ ] ..................... integer parameters for RandomEff constructor (needed only if there are some random effects)
//                                if version = 32: integer argument for RandomEff32::RE initializer
// objbD[ ] ..................... double parameters for RandomEff constructor (needed only if there are some random effects)
//                                *only space is needed, it can be filled by whatever
//                                *if version = 32: space of length nCluster 
// b_GsplI[] .................... needed G-spline parameters for random intercept (if its distribution is a G-spline)
//    b_GsplI[0]          = dim (here only 1 is allowed)
//    b_GsplI[1]          = total_length

// gridA[sum(ngrid)] ............ grids to compute predictive quantities for each observation
// loggridA[sum(ngrid)] ......... logarithm of the grid
// ngrid[nobs] .................. lengths of grids for each observation
// onlyAver ..................... 0/1: compute only predictive quantities or return values as well?
// predictP[4] .................. 0/1 indicating which predictive quantities are to be computed
//   predictP[0] ... densities?
//   predictP[1] ... survivor functions?
//   predictP[2] ... hazards?
//   predictP[3] ... cumulative hazards?
// M ............................ McMC sample size (total, 'skip' and 'by' iterations included)
//                                * M should be <= number of rows in *.sim files
//                                * here: it is an index of the last iteration used to compute the average
// skip ......................... how many rows are to be skipped at the beginning of the sample
// by ........................... only every 'by' G-spline will be taken into account
// nwrite ....................... frequency of informing the user about the progress
// version ...................... arbitrary or 32
//                                if = 32, then model for doubly-interval censored data is assumed with G-spline errors
//                                and bivariate normal random intercepts in the onset and time-to-event parts of the model
// Onset ........................ only used by version = 32
//                                equal to 1 if we are predicting the onset 
//                                equak to 0 if we are predicting the event
// errP ......................... error flag (0 on output if everything OK)
//
void
predictive_GS(double *averDens,         double *averS,           double *averHaz,      double *averCumHaz,
              double *valDens,          double *valS,            double *valHaz,       double *valCumHaz,
              double *quantDens,        double *quantS,          double *quantHaz,     double *quantCumHaz,
              const int *dimsP,         const double *X,         const int *obsdims,
              int *M_now,               char **dirP,             char **extensP,       char **extens_adjP,
              const int *GsplI,
              const int *objBetaI,      const double *objBetaD,
              const int *objbI,         const double *objbD,
              const int *b_GsplI,
              const double *gridA,      const double *loggridA,  const int *ngrid,
              double *probsA,           const int *nquant,       int *onlyAver,        const int *predictP,
              const int *M,             const int *skip,         const int *by,        const int *nwrite,
              const int *version,       const int *Onset,        int *errP)
{
  try{
    GetRNGstate();
    double dtemp;
    int itemp;

    int i, j, ix;
    double tmpd;

    *errP = 0;
    string dir = *dirP;
    string extens = *extensP;
    string extens_adj = *extens_adjP;

    /*** Dimensionality parameters ***/
    const int *nobs     = dimsP;
    const int *ncluster = dimsP + 1;
    const int *nwithin  = dimsP + 2;
    const int M_now_max = *M_now;

    /*** What to predict? ***/
    const int *predDens   = predictP + 0;
    const int *predS      = predictP + 1;
    const int *predHaz    = predictP + 2;
    const int *predCumHaz = predictP + 3;

    /*** Quantiles ***/
    if (*nquant <= 0) *onlyAver = 1;

    /*** Needed G-spline parameters      ***/
    const int *dim          = GsplI + 0;
    const int *total_length = GsplI + 1;
    const int *GsplK        = GsplI + 2;          /* K1 (and K2) */
    int *Glength            = (int*) calloc(*dim, sizeof(int));
    if (!Glength) throw returnR("Not enough memory available in predictive_GS (Glength)", 1);
    for (j = 0; j < *dim; j++) Glength[j] = 2*GsplK[j] + 1;

    /*** Check obsdims ***/
    for (i = 0; i < *nobs; i++){
      if (obsdims[i] < 0 || obsdims[i] >= *dim) throw returnR("Error: Inconsistent 'obsdims' parameter supplied to predictive_GS", 1);
    }

    /*** Check grid  and log-grid ***/
    int sum_ngrid = 0;
    for (i = 0; i < *nobs; i++) sum_ngrid += ngrid[i];

    /*** Object for regression parameters ***/
    BetaGamma* beta = new BetaGamma;
    if (!beta) throw returnR("Not enough memory available in predictive_GS (beta)", 1);
    *beta = BetaGamma(objBetaI, objBetaD);

    /*** Object for random effects      ***/
    RandomEff *bb        = new RandomEff;
    RandomEff32::RE *db  = new RandomEff32::RE;
    bool reff_NORMAL = true;                      /** BUT NOT version = 32 !**/

    /*** Objects for bivariate normal random effects in version = 32 ***/
    double *dval, *bval, *dbval;
    double D32[7] = {1, 0, 1,  2,  1, 0, 1};                /** parD argument for RandomEff32::RE initializer  (filled arbitrary) **/
    if (*version == 32){
      reff_NORMAL = false;
      dval = (double*) calloc(objbI[2], sizeof(double));        // objbI[2] = nCluster
      bval = (double*) calloc(objbI[2], sizeof(double));        // objbI[2] = nCluster
      RandomEff32::init(db, dval, bval, D32, objbI, objbI);
      if (*Onset) dbval = dval;
      else        dbval = bval;
    }
    else{
      if (beta->nRandom()){
        *bb = RandomEff(objbI, objbD);
        if (bb->type_prior() == Gspline_) reff_NORMAL = false;
      }
    }

    
    /*** Object for covariance matrix of random effects                              ***/
    /*** or arrays for G-spline parameters definig distribution of random effects    ***/
    CovMatrix *DD = new CovMatrix;
    const int nD = (beta->nRandom() * (beta->nRandom() + 1)) / 2;

    const int *dim_b          = b_GsplI + 0;
    const int *total_length_b = b_GsplI + 1;
    int k_effect_b;
    int *rM_b = &itemp;
    double *cum_w_b = &dtemp;
    double *sig_scale_b = &dtemp;
    double *prop_mu_b = &dtemp;

    if (*version != 32){
      if (beta->nRandom()){
        if (reff_NORMAL){
          int DDparmI[2];
          DDparmI[0] = beta->nRandom();
          DDparmI[1] = InvWishart;                                       /** it does not matter what is filled here **/
          double *DDparmD = (double*) calloc(2*nD + 1, sizeof(double));
          if (!DDparmD) throw returnR("Not enough memory available in predictive_GS (DDparmD)", 1);
          for (j = 0; j < beta->nRandom(); j++){                         /** initial cov matrix and scale matrix equal to identity  **/
            ix = (j * (2*beta->nRandom() - j + 1))/2;                    /** again, it does not matter what is filled here          **/
            DDparmD[ix] = DDparmD[nD + 1 + ix] = 1.0;                    /** initial matrix must only be positive definite          **/
            for (i = j+1; i < beta->nRandom(); i++){                     /** to pass the CovMatrix constructor                      **/
              DDparmD[ix + i - j] = DDparmD[nD + 1 + ix + i - j] = 0.0;
            }
          }
          DDparmD[nD] = beta->nRandom() + 2;                              /** 'prior degrees of freedom', it does not matter what   **/
          *DD = CovMatrix(DDparmI, DDparmD);
          free(DDparmD);
        }
        else{                          /** G-spline random effects **/
          cum_w_b     = (double*) calloc(*total_length_b, sizeof(double));
          prop_mu_b   = (double*) calloc(*total_length_b, sizeof(double));
          sig_scale_b = (double*)  calloc(*dim_b, sizeof(double));
          rM_b        = (int*) calloc(*ncluster, sizeof(int));
          if (!cum_w_b || !prop_mu_b || !sig_scale_b) throw returnR("Not enough memory available in predictive_GS (cum_w_b/sig_scale_b)", 1);
          if (!rM_b) throw returnR("Not enough memory available in predictive_GS (rM_b)", 1);
        }
      }
    }  /** end of if (*version != 32) **/

    /*** Space for linear predictors    ***/
    double *linPred = (double*) calloc(*nobs, sizeof(double));
    if (!linPred) throw returnR("Not enough memory available in predictive_GS (linPred)", 1);
    for (i = 0; i < *nobs; i++) linPred[i] = 0.0;

    /*** Allocate memory for needed quantities from simulated G-splines ***/
    int k_effect;
    double *sigma         = (double*)  calloc(*dim, sizeof(double));
    double *gamma         = (double*)  calloc(*dim, sizeof(double));
    double *delta         = (double*)  calloc(*dim, sizeof(double));
    double *intcpt        = (double*)  calloc(*dim, sizeof(double));
    double *scale         = (double*)  calloc(*dim, sizeof(double));
    double *delta_sig     = (double*)  calloc(*dim, sizeof(double));
    double *inv_sig_scale = (double*)  calloc(*dim, sizeof(double));
    if (!sigma || !gamma || !delta || !inv_sig_scale || !intcpt || !scale || !delta_sig) 
      throw returnR("Not enough memory available in predictive_GS (sigma/gamma/delta/intcpt/scale/delta_sig/inv_sig_scale)", 1);

    double **w_marg           = (double**) calloc(*dim, sizeof(double*));
    double **mu_sig_marg      = (double**) calloc(*dim, sizeof(double*));
    if (!w_marg || !mu_sig_marg)
      throw returnR("Not enough memory available in predictive_GS (w_marg/sc_mu_marg)", 1);
    for (j = 0; j < *dim; j++){
      w_marg[j]      = (double*) calloc(Glength[j], sizeof(double));
      mu_sig_marg[j] = (double*) calloc(Glength[j], sizeof(double));
      if (!w_marg[j] || !mu_sig_marg[j]) throw returnR("Not enough memory available in predictive_GS (w_marg[j]/mu_sig_marg[j])", 1);
    }

    /*** Open files with simulated G-splines ***/
    std::string kpath     = dir + "/mixmoment" + extens + ".sim";    
    std::string wpath     = dir + "/mweight" + extens + ".sim";
    std::string mupath    = dir + "/mmean" + extens + ".sim";       
    std::string sigmapath = dir + "/gspline" + extens + ".sim";
    std::ifstream kfile, wfile, mufile, sigmafile;
    openGsplineFiles(kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath, *skip + 1);   /* skip also header */

    /*** Open files with simulated remaining quantities ***/
    std::string betapath = dir + "/beta" + extens + ".sim";
    std::ifstream betafile;

    std::string Dpath    = dir + "/D" + extens + ".sim";
    std::ifstream Dfile;

    std::string D32path  = dir + "/D" + ".sim";
    std::ifstream D32file;

    std::string kpath_b     = dir + "/mixmoment" + extens_adj + ".sim";    
    std::string wpath_b     = dir + "/mweight" + extens_adj + ".sim";
    std::string mupath_b    = dir + "/mmean" + extens_adj + ".sim";       
    std::string sigmapath_b = dir + "/gspline" + extens_adj + ".sim";
    std::ifstream kfile_b, wfile_b, mufile_b, sigmafile_b;

    openRegresFiles(betafile, Dfile, betapath, Dpath, *skip + 1, beta->nbeta(), beta->nRandom(), reff_NORMAL);    /* skip also header */
    if (*version == 32){     
      openD32File(D32file, D32path, *skip + 1);     /* skip also header */
    }
    else{    
      if (beta->nRandom() && !reff_NORMAL){
        openGsplineFiles(kfile_b, wfile_b, mufile_b, sigmafile_b, kpath_b, wpath_b, mupath_b, sigmapath_b, *skip + 1);
      }
    }

    /*** Reset averages ***/
    resetAverage(averDens, nobs, ngrid, predDens);
    resetAverage(averS, nobs, ngrid, predS);
    resetAverage(averHaz, nobs, ngrid, predHaz);
    resetAverage(averCumHaz, nobs, ngrid, predCumHaz);

    /*** Loop over McMC iterations ***/
    double *vvDens   = valDens;
    double *vvS      = valS;
    double *vvHaz    = valHaz;
    double *vvCumHaz = valCumHaz;
    const int *shift_pointer_inEval = (*onlyAver ? &ONE_INT : &M_now_max);
    
    if (*skip >= *M) throw returnR("More McMC iterations should be skipped than available", 1);    
    readGsplineFromFiles2(&k_effect, w_marg, mu_sig_marg, gamma, sigma, delta, intcpt, scale, delta_sig, 0, *skip, 
                          *dim, *total_length, GsplK,
                          kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
    readRegresFromFiles(beta, DD, 0, *skip, betafile, Dfile, betapath, Dpath, reff_NORMAL);
    if (*version == 32){
      readDfromFile(db, 0, *skip, D32file, D32path);
      predict_db(db);
      linPred_GS(linPred, beta, dbval, X, nwithin, nobs, ncluster);
    }
    else{
      if (beta->nRandom()){
        if (reff_NORMAL){
          bb->predictNormalRE(beta, DD);
        }
        else{
          readGsplineFromFiles3(&k_effect_b, cum_w_b, prop_mu_b, sig_scale_b, 0, *skip, *dim_b, *total_length_b,
                                kfile_b, wfile_b, mufile_b, sigmafile_b, kpath_b, wpath_b, mupath_b, sigmapath_b);
          bb->predictGspl_intcpt(&k_effect_b, cum_w_b, prop_mu_b, sig_scale_b, rM_b);
        }
      }
      linPred_GS(linPred, beta, bb->bMP(), X, nwithin, nobs, ncluster);
    }
    evalPredFuns(averDens, averS, averHaz, averCumHaz, vvDens, vvS, vvHaz, vvCumHaz, obsdims, nobs, ngrid, gridA, loggridA,
                 linPred, dim, Glength, w_marg, mu_sig_marg, intcpt, sigma, scale, inv_sig_scale, predictP, &_zero_weight,
                 shift_pointer_inEval);

    *M_now = 1;

    int by_1 = *by - 1;
    int jump_value = (*onlyAver ? 0 : 1);
    int backs = 0;
    Rprintf("Iteration ");
    for (int iter = *skip + 1 + (*by); iter <= *M; iter += (*by)){
      if (*M_now >= M_now_max) throw returnR("Error: Higher sample size would be used than indicated", 1);
      if (*predDens)   vvDens += jump_value;
      if (*predS)      vvS += jump_value;
      if (*predHaz)    vvHaz += jump_value;
      if (*predCumHaz) vvCumHaz += jump_value;
      readGsplineFromFiles2(&k_effect, w_marg, mu_sig_marg, gamma, sigma, delta, intcpt, scale, delta_sig, by_1, iter, 
                            *dim, *total_length, GsplK,
                            kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
      readRegresFromFiles(beta, DD, by_1, iter, betafile, Dfile, betapath, Dpath, reff_NORMAL);
      if (*version == 32){
        readDfromFile(db, by_1, iter, D32file, D32path);
        predict_db(db);
        linPred_GS(linPred, beta, dbval, X, nwithin, nobs, ncluster);
      }
      else{
        if (beta->nRandom()){
          if (reff_NORMAL){
            bb->predictNormalRE(beta, DD);
          }
          else{
            readGsplineFromFiles3(&k_effect_b, cum_w_b, prop_mu_b, sig_scale_b, by_1, iter, *dim_b, *total_length_b,
                                  kfile_b, wfile_b, mufile_b, sigmafile_b, kpath_b, wpath_b, mupath_b, sigmapath_b);
            bb->predictGspl_intcpt(&k_effect_b, cum_w_b, prop_mu_b, sig_scale_b, rM_b);
          }
        }
        linPred_GS(linPred, beta, bb->bMP(), X, nwithin, nobs, ncluster);
      }
      evalPredFuns(averDens, averS, averHaz, averCumHaz, 
                   vvDens, vvS, vvHaz, vvCumHaz, 
                   obsdims, nobs, ngrid, gridA, loggridA,
                   linPred, dim, Glength, w_marg, mu_sig_marg, intcpt, sigma, scale, inv_sig_scale, predictP, &_zero_weight,
                   shift_pointer_inEval);
      (*M_now)++;

      if (!(iter % (*nwrite)) || iter == *M){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
      }
    }    /** end of the while over iterations **/
    Rprintf("\n");
         
    /*** Close files with simulated G-splines and regression quantities ***/
    closeGsplineFiles(kfile, wfile, mufile, sigmafile);
    closeRegresFiles(betafile, Dfile, beta->nbeta(), beta->nRandom(), reff_NORMAL);
    if (*version == 32){
      D32file.close();
    }
    else{
      if (beta->nRandom() && !reff_NORMAL) closeGsplineFiles(kfile_b, wfile_b, mufile_b, sigmafile_b);
    }

    /*** McMC averages ***/
    cumsum2average(averDens, M_now, nobs, ngrid, predDens);
    cumsum2average(averS, M_now, nobs, ngrid, predS);
    cumsum2average(averHaz, M_now, nobs, ngrid, predHaz);
    cumsum2average(averCumHaz, M_now, nobs, ngrid, predCumHaz);

    /*** Indeces of quantile values in sampled chain (indexing starting from 0)   ***/
    // indquant1, indquant2 ..... quantile = q*sample[indquant1] + (1-q)sample[indquant2]
    //
    int *indquant1 = &itemp;
    int *indquant2 = &itemp;
    if (!(*onlyAver)){
      indquant1  = (int*) calloc(*nquant, sizeof(int));
      indquant2  = (int*) calloc(*nquant, sizeof(int));
      if (!indquant1 || !indquant2) throw returnR("Error Not enough memory available in predictive_GS (indquant1/indquant2)", 1);
      for (i = 0; i < *nquant; i++){
        if (probsA[i] < 0 || probsA[i] > 1) throw returnR("Error: Incorrect probs values supplied.", 1);
        if (probsA[i] <= 0) indquant1[i] = indquant2[i] = 0;
        else{
          if (probsA[i] >= 1) indquant1[i] = indquant2[i] = *M_now - 1;
          else{
            tmpd = probsA[i] * double(*M_now);
            if (fabs(tmpd - floor(tmpd + 1e-8)) < 1e-8){
              indquant1[i] = int(floor(tmpd)) - 1;
              indquant2[i] = int(floor(tmpd));
            }
            else{
              indquant1[i] = indquant2[i] = int(floor(tmpd));
            }
          }
	}
      }
      Rprintf("\nComputing quantiles.");
      value2quantile(valDens, quantDens, probsA, indquant1, indquant2, nquant, M_now, nobs, ngrid, predDens, shift_pointer_inEval);
      value2quantile(valS, quantS, probsA, indquant1, indquant2, nquant, M_now, nobs, ngrid, predS, shift_pointer_inEval);
      value2quantile(valHaz, quantHaz, probsA, indquant1, indquant2, nquant, M_now, nobs, ngrid, predHaz, shift_pointer_inEval);
      value2quantile(valCumHaz, quantCumHaz, probsA, indquant1, indquant2, nquant, M_now, nobs, ngrid, predCumHaz, shift_pointer_inEval);
    }

    PutRNGstate();

    /*** Cleaning ***/
    if (!(*onlyAver)){
      free(indquant1);
      free(indquant2);
    }

    for (j = 0; j < *dim; j++){
      free(w_marg[j]);   
      free(mu_sig_marg[j]);
    }    
    free(w_marg);    
    free(mu_sig_marg);
    free(sigma);     
    free(gamma);   
    free(delta);     
    free(inv_sig_scale);
    free(intcpt);    
    free(scale);   
    free(delta_sig);
    free(Glength);

    free(linPred);

    delete DD;
   
    if (*version == 32){
      free(bval);
      free(dval);
    }
    else{
      if (beta->nRandom()){
        if (reff_NORMAL){
	  //        delete DD;
        }
        else{
          free(sig_scale_b);
          free(prop_mu_b);
          free(cum_w_b);
          free(rM_b);
        }
      }
    }
    delete db;
    delete bb;
    delete beta;

    return; 
  }
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }
}

}    /** end of extern "C" **/


// =============================================================================
// ******* evalPredFuns 
// =============================================================================
//
// obsdims[nobs] ................ indicators of a dimension for each observation
//                                * if dim == 1, this should be a vector filled with zeros only
//                                * if dim == 2, this should be a vector of 0 and 1
// nobs ......................... number of observations
// ngrid[nobs] .................. length of the grid for each observation
// gridA[sum(ngrid)] ............ ! already adjusted for possible t_0
// loggridA[sum(ngrid)] ......... log(gridA), also already adjusted for possible t_0
// linPred[nobs] ................ linear predictors (x'beta + z'b) for required covariate combinations
//                      * linear predictor without possible intercept!
// dim .......................... dimension of the G-spline
// Glength[dim] ................. length of the G-spline in each dimension
// zero_weight .................. limit to ignore a mixture component
// shift_pointer ................ a value to be used to shift the pointer in val* arrays to go to the following grid point
//     usually: if (onlyAver) shift_pointer = 1
//              else          shift_pointer = M_now_max
//
void
evalPredFuns(double* averDens,        double* averS,               double* averHaz,            double* averCumHaz,
             double* valDens,         double* valS,                double* valHaz,             double* valCumHaz,
             const int* obsdims,      const int* nobs,             const int* ngrid,           
             const double* gridA,     const double* loggridA,
             const double* linPred,   const int* dim,              const int* Glength,
             double** const w_marg,   double** const mu_sig_marg,  const double* intcpt,
             const double* sigma,     const double* scale,         double* inv_sig_scale,
             const int* predictP,     const double* zero_weight,   const int* shift_pointer)
{
  static int i, j, ix;
  static double stres1, stres2;

  /** Compute some transformations of sigma's **/
  for (j = 0; j < *dim; j++) inv_sig_scale[j] = 1/(sigma[j]*scale[j]);

  const int shift_pointer_Dens = predictP[0] ? (*shift_pointer) : 0;
  const int shift_pointer_S = predictP[1] ? (*shift_pointer) : 0;

  double* avDens   = averDens;
  double* avS      = averS;
  double* avHaz    = averHaz;
  double* avCumHaz = averCumHaz;
  double* vvDens   = valDens;
  double* vvS      = valS;
  double* vvHaz    = valHaz;
  double* vvCumHaz = valCumHaz;
  const double* grid     = gridA;
  const double* loggrid  = loggridA;
  const double* lin_pr   = linPred;
  const int* obs_dim     = obsdims;
  for (i = 0; i < *nobs; i++){
    for (ix = 0; ix < ngrid[i]; ix++){
      stres1 = inv_sig_scale[*obs_dim] * ((*loggrid) - intcpt[*obs_dim] - (*lin_pr));      

      /** Predictive density **/
      if (predictP[0] || predictP[2]){

        *vvDens = 0.0;
        for(j = 0; j < Glength[*obs_dim]; j++){
          if (w_marg[*obs_dim][j] <= *zero_weight) continue;
          stres2 = stres1 - mu_sig_marg[*obs_dim][j];
          if (stres2 <= -_stres0 || stres2 >= _stres0) continue;
          *vvDens += w_marg[*obs_dim][j] * dnorm(stres2, 0, 1, 0);
        }
        *vvDens *= inv_sig_scale[*obs_dim]/(*grid);

        if (predictP[0]){
          *avDens += (*vvDens);
          avDens++;
        }
        vvDens += shift_pointer_Dens;
      }

      /** Predictive survivor function **/
      if (predictP[1] || predictP[2] || predictP[3]){

        *vvS = 0.0;
        for(j = 0; j < Glength[*obs_dim]; j++){
          if (w_marg[*obs_dim][j] <= *zero_weight) continue;
          stres2 = stres1 - mu_sig_marg[*obs_dim][j];
          if (stres2 <= -_stres0) *vvS += w_marg[*obs_dim][j];
          else{
            if (stres2 < _stres0){
              *vvS += w_marg[*obs_dim][j] * pnorm(stres2, 0, 1, 0, 0);
	    }
          }
        }
        
        if (predictP[1]){
          *avS += (*vvS);
          avS++;
        }
        vvS += shift_pointer_S;
      }

      /** Predictive hazard            **/
      if (predictP[2]){
        if (*(vvS-shift_pointer_S) <= ZERO) *vvHaz = FLT_MAX;
        else                                *vvHaz = (*(vvDens-shift_pointer_Dens))/(*(vvS-shift_pointer_S)); 
        if (*vvHaz >= FLT_MAX) *vvHaz = FLT_MAX;
        *avHaz += (*vvHaz);

        avHaz++;
        vvHaz += (*shift_pointer);
      }

      /** Predictive cumulative hazard  **/
      if (predictP[3]){
        if (*(vvS-shift_pointer_S) <= ZERO) *vvCumHaz = FLT_MAX;
        else                                *vvCumHaz = -log(*(vvS-shift_pointer_S));
        if (*vvCumHaz <= 0)       *vvCumHaz *= (-1);
        if (*vvCumHaz >= FLT_MAX) *vvCumHaz = FLT_MAX;
        *avCumHaz += (*vvCumHaz);
 
        avCumHaz++;
        vvCumHaz += (*shift_pointer);
      }
          
      grid++;
      loggrid++;
    }
    lin_pr++;
    obs_dim++;
  }

  return;
}
