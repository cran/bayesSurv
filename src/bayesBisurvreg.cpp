// Main function to run McMC for possibly bivariate AFT model
//   without random effects
// * allow also for doubly censoring
//
// 28/12/2004: start working on it (ve vlaku Veseli n.L. - Praha)
// 12/01/2005: more or less finished
// 25/04/2005: option mainSimul added to write simulated values during burn-in only after nwrite iteration
//

#include "bayesBisurvreg.h"

extern "C" {

// n    = sample size (number of (bi)variate observational vectors)
// dim  = dimension of the sample (1 or 2)
// p1   = number of regressors (intercept excluded) for the event time (or onset time in the case of doubly-censored data)
// p2   = number of regressors (intercept excluded) for the event time in the case of doubly-censored data

// dirP ................... path to the directory where to store simulated values
// dimsP[] ................ dimension parameters
//    dimsP[0] = sample size (n)
//    dimsP[1] = is doubly-censored? (0 = no, 1 = yes) 
//
// X1[dim*n*p1] ........... covariates for the event time (or onset time in the case of doubly-censored data)
// X2[dim*n*p2] ........... covariates for the event time in the case of doubly-censored data
//    !!! CHANGE COMPARED TO EARLIER VERSIONS OF BAYESSURVREG (bayessurvreg1): X1 and X2 are stored in ROW major order, i.e.
//        X[0, ..., p-1]  ... covariates for the first observation
//        X[p, ..., 2p-1] ... covariates for the second observation etc.
//
// y1_left[dim*n] ......... log-event time or log-onset time - lower limit
// y1_right[dim*n] ........ log-event time or log-onset time - upper limit
// status1[dim*n] ......... censoring indicator vector for event time or onset time
//
// t2_left[dim*n] ......... event-time in the case of doubly-censored data - lower limit
// t2_right[dim*n] ........ event-time in the case of doubly-censored data - upper limit
// status2[dim*n] ......... censoring indicator vector for event time in the case of doubly-censored data
//
// Y1[dim*n] .............. (augmented) log-event time or log-onset time
// Y2[dim*n] .............. (augmented) log-time-to-event (i.e. log(T2 - T1))  in the case of doubly censored data
//
// r1[n] .................. allocation indicators for event times or onset times  (R indeces, i.e. from 1)
// r2[n] .................. allocation indicators for event times in the case of doubly-censored data (R indeces, i.e. from 1)
// specif[2] .............. G-spline specification for onset and time-to-event time (see bayesHistogram.cpp for explanation)
// GsplineI1[] ............ integer G-spline defining parameters for distribution of log-event times or log-onset times
// GsplineD1[] ............ double G-spline defining parameters for distribution of log-event times or log-onset times
// GsplineI2[] ............ integer G-spline defining parameters for distribution of log-onset times in the case of doubly-censored data
// GsplineD2[] ............ double G-spline defining parameters for distribution of log-onset times in the case of doubly-censored data
//
// priorBeta1I[] .......... integer parameters defining prior for beta1 (regression parameter for the log-event time or log-onset time)
// priorBeta1D[] .......... double parameters defining inits and prior for beta1
//
// priorBeta2I[] .......... integer parameters defining prior for beta2 (regres. par. for the log-event time in the case of doubly cens.)
// priorBeta2D[] .......... double parameters defining inits and prior for beta2
//
// iterM .................. index of the first iteration
// nsimulP[3].............. length of simulation, thinning etc.
// storeP[6] .............. parameters defining what should be stored
// mainSimul[1] ........... 0/1, if 0 then all parameters are written to files only every nwrite iteration
// errP ................... error flag
//
void
bayesBisurvreg(char** dirP,             const int* dimsP,         const double* X1,     const double* X2,
               const double* y1_left,   const double* y1_right,   const int* status1,
               const double* t2_left,   const double* t2_right,   const int* status2,               
               double* Y1,              double* Y2,
               int* r1,                 int* r2,                  const int* specif,
               int* GsplineI1,          double* GsplineD1,
               int* GsplineI2,          double* GsplineD2,
               int* priorBeta1I,        double* priorBeta1D,
               int* priorBeta2I,        double* priorBeta2D,
               int* iterM,              int* nsimulP,             int* storeP,
               const int* mainSimul,    int* errP)
{
  try{
    int out_format[2] = {6, 1};                /** precision and width for output **/

    double dtemp;
    int itemp;
    double* ddtemp = (double*) malloc(sizeof(double));

    int i, j, k;
    GetRNGstate();

    *errP = 0;
    std::string dir = *dirP;
 
    /** Numbers of iterations etc.  **/
    int niter = nsimulP[0]; 
    int nthin = nsimulP[1];
    int nwrite = nsimulP[2];

    /** Mandatory storing of simulated values **/
    int* storea1P       = storeP;
    int* storey1P       = storeP + 1;
    int* storer1P       = storeP + 2;
    int* storea2P       = storeP + 3;
    int* storey2P       = storeP + 4;
    int* storer2P       = storeP + 5;

    /** Dimensionality parameters  **/        
    const int* nP     = dimsP;                                      /* number of observational uni- or bivariate vectors  */
    const int* doubly = dimsP + 1;                                  /* 0/1 .... is doubly censored?                       */
    const int* dimP   = GsplineI1 + 0;                              /* dimension of the response (1 or 2)                 */
    const int nobs    = (*nP) * (*dimP);                            /* number of observations                             */
    if (!(*doubly)) *storea2P = *storey2P = *storer2P = 0;
    if (*doubly){
      if (GsplineI2[0] != (*dimP)) throw returnR("Error: Inconsistent G-spline dimensions", 1);
    }


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

    /** Object for beta1           **/
    BetaGamma* beta1 = new BetaGamma;
    if (!beta1) throw returnR("Not enough memory available in bayesBisurvreg (beta1)", 1);
    if (priorBeta1I[0]){
      *beta1 = BetaGamma(priorBeta1I, priorBeta1D);
    }

    /** Object for beta2           **/
    BetaGamma* beta2 = new BetaGamma;
    if (!beta2) throw returnR("Not enough memory available in bayesBisurvreg (beta2)", 1);
    if (*doubly && priorBeta2I[0]){
      *beta2 = BetaGamma(priorBeta2I, priorBeta2D);
    }

    /** Regression residuals       **/
    double* regresRes1 = (double*) calloc(nobs, sizeof(double));
    if (!regresRes1) throw returnR("Not enough memory available in bayesBisurvreg (regresRes1)", 1);
    regresRes_GS(regresRes1, Y1, beta1, &ZERO, X1, &ZERO_INT, &nobs, &ZERO_INT);

    double* regresRes2 = &dtemp;
    if (*doubly){
      regresRes2 = (double*) calloc(nobs, sizeof(double));
      if (!regresRes2) throw returnR("Not enough memory available in bayesBisurvreg (regresRes2)", 1);
      regresRes_GS(regresRes2, Y2, beta2, &ZERO, X2, &ZERO_INT, &nobs, &ZERO_INT);
    }


    /** Some manipulations with the design matrices  **/
    const double*  XNow;
    double* XXtNow = &dtemp;

    double* XXtb1 = &dtemp;
    if (beta1->nbeta()){
      XXtb1 = (double*) calloc(nobs * beta1->lcovFixed(), sizeof(double));
      if (!XXtb1) throw returnR("Not enough memory available in bayesBisurvreg (XXtb1)", 1);

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

    double* XXtb2 = &dtemp;
    if (*doubly && beta2->nbeta()){
      XXtb2 = (double*) calloc(nobs * beta2->lcovFixed(), sizeof(double));
      if (!XXtb2) throw returnR("Not enough memory available in bayesBisurvreg (XXtb2)", 1);

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

    /** G-spline objects **/
                              /* _null_weight ...... limit for the mixture component weight to be recorded in a file */
    Gspline* g_y1 = new Gspline;
    Gspline* g_y2 = new Gspline;
    if (!g_y1 || !g_y2) throw returnR("Not enough memory available in bayesBisurvreg (g_y1/g_y2)", 1);
    *g_y1 = Gspline(GsplineI1, GsplineD1);

    const double _null_weight1 = _null_mass/g_y1->total_length();        
    double _null_weight2;
    if (*doubly){
      *g_y2 = Gspline(GsplineI2, GsplineD2);
      _null_weight2 = _null_mass/g_y2->total_length();        
    }

    /** Arrays to store numbers of observations belonging to each mixture component  **/
    /** Subtract 1 from each rM (R -> C++ indeces) and check for consistency         **/
    int* mixtureN1 = (int*) calloc(g_y1->total_length(), sizeof(int));
    if (!mixtureN1) throw returnR("Not enough memory available in bayesBisurvreg (mixtureN1)", 1);
    for (i = 0; i < g_y1->total_length(); i++) mixtureN1[i] = 0;
    for (i = 0; i < *nP; i++){
      r1[i]--;
      if (r1[i] < 0 || r1[i] >= g_y1->total_length()) throw returnR("Inconsistent initial r1 supplied", 1);
      mixtureN1[r1[i]]++;
    }

    int* mixtureN2 = &itemp;
    if (*doubly){
      mixtureN2 = (int*) calloc(g_y2->total_length(), sizeof(int));
      if (!mixtureN2) throw returnR("Not enough memory available in bayesBisurvreg (mixtureN2)", 1);
      for (i = 0; i < g_y2->total_length(); i++) mixtureN2[i] = 0;
      for (i = 0; i < *nP; i++){
        r2[i]--;
        if (r2[i] < 0 || r2[i] >= g_y2->total_length()) throw returnR("Inconsistent initial r2 supplied", 1);
        mixtureN2[r2[i]]++;
      }
    }


    /** Open files to store simulated values for writing **/
    std::string iterpath      = dir + "/iteration.sim";
    std::string betapath      = dir + "/beta.sim"; 
    std::string betapath2     = dir + "/beta_2.sim"; 

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

    std::ofstream iterfile, betafile, betafile2;
    std::ofstream sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile;
    std::ofstream sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2;

    openFile(iterfile, iterpath, 'a');
    if (beta1->nbeta()) openFile(betafile, betapath, 'a');
    if (beta2->nbeta()) openFile(betafile2, betapath2, 'a');
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

    /***** Working arrays and variables                                                             *****/
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

    int* iwork1         = (int*) calloc(l_iwork1, sizeof(int));
    double* dwork1      = (double*) calloc(l_dwork1, sizeof(double));
    double** mu1        = (double**) calloc(g_y1->dim(), sizeof(double*));
    double* log_poster1 = (double*) calloc(l_log_poster1, sizeof(double));
    if (!iwork1 || !dwork1 || !mu1 || !log_poster1) throw returnR("Not enough memory available in bayesBisurvreg", 1);
    for (j = 0; j < g_y1->dim(); j++){
      mu1[j] = (double*) calloc(g_y1->length(j), sizeof(double));
      if (!mu1[j]) throw returnR("Not enough memory available in bayesBisurvreg", 1);
    }

    int* iwork2 = &itemp;
    double* dwork2 = &dtemp;
    double** mu2 = &ddtemp; 
    double* log_poster2 = &dtemp;
    if (*doubly){
      iwork2      = (int*) calloc(l_iwork2, sizeof(int));
      dwork2      = (double*) calloc(l_dwork2, sizeof(double));
      mu2         = (double**) calloc(g_y2->dim(), sizeof(double*));
      log_poster2 = (double*) calloc(l_log_poster2, sizeof(double));
      if (!iwork2 || !dwork2 || !mu2 || !log_poster2) throw returnR("Not enough memory available in bayesBisurvreg", 1);
      for (j = 0; j < g_y2->dim(); j++){
        mu2[j] = (double*) calloc(g_y2->length(j), sizeof(double));
        if (!mu2[j]) throw returnR("Not enough memory available in bayesBisurvreg", 1);
      }
    }

    int a_ipars[2];
    a_ipars[0] = *nP;


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
    Rprintf("Iteration ");    
    for (iter = nullthIter + 1; iter <= lastIter; iter++){
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
        }
      }    /** end of the thinning cycle  **/


      /**  Write to files  **/
      if (!(iter % nwrite) || iter == lastIter){
        writeAll = 1;
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
      }

      if ((*mainSimul) || writeAll){
        writeToFile_1(&iter, 1, iterfile, out_format[0], out_format[1]);      
        if (beta1->nbeta()) writeToFile_1(beta1->betaP(), beta1->nbeta(), betafile, out_format[0], out_format[1]);
        if (beta2->nbeta()) writeToFile_1(beta2->betaP(), beta2->nbeta(), betafile2, out_format[0], out_format[1]);
        writeToFiles_bayesHistogram(g_y1, r1, Y1, log_poster1,
	    l_moments, l_lambda1, l_log_poster1, nP, storea1P, storey1P, storer1P, n1_censored, &writeAll, iwork1, dwork1,
            sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
            _null_weight1, out_format[0], out_format[1]);
        if (*doubly){
          writeToFiles_bayesHistogram(g_y2, r2, Y2, log_poster2,
	      l_moments, l_lambda2, l_log_poster2, nP, storea2P, storey2P, storer2P, 1, &writeAll, iwork2, dwork2,
              sigmafile2, lambdafile2, mixmomentfile2, mweightfile2, mlogweightfile2, mmeanfile2, Yfile2, rfile2, logposterfile2,
              _null_weight2, out_format[0], out_format[1]);
        }
        writeAll = 0;
      }
    }    /** end of the main cycle over iter **/

    iterfile.close();
    if (beta1->nbeta()) betafile.close();
    if (beta2->nbeta()) betafile2.close();
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

    /*** Cleaning ***/
    for (j = 0; j < g_y1->dim(); j++) free(mu1[j]);
    free(iwork1);    free(dwork1);    free(mu1);   free(log_poster1);
    if (*doubly){
      for (j = 0; j < g_y2->dim(); j++) free(mu2[j]);
      free(iwork2);    free(dwork2);    free(mu2);   free(log_poster2);
    }

    free(mixtureN1);                       if (*doubly) free(mixtureN2);
    delete g_y1;                           delete g_y2;
    if (beta1->nbeta()) free(XXtb1);       if (*doubly && beta2->nbeta()) free(XXtb2);
    free(regresRes1);                      if (*doubly) free(regresRes2);
    delete beta1;                          delete beta2;

    free(ddtemp);

    PutRNGstate();
    return;
  }    
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }
}

}    /*** end of extern "C"  ***/
