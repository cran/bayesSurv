// Bayesian estimate of a uni- or bivariate density using G-splines
//  * main function to be called from R
//
// 21/10/2004: 'bayesHistogram'
// 03/11/2004: 'bayesHistogram' rewritten such that sampled values are directly written to files
// 07/06/2005: during burn-in not all sampled values are writen to files
//
// ==============================================================================================
//
#include "bayesHistogram.h"

extern "C" {

using namespace std;

// dimsP ......... dimsP[0] = sample size (N)
// y_left ........ long vector with observations y[i, j] = y_left[i*_dim + j], 
//         where i = 0, ..., N-1, j = 0, ... _dim-1,
//         i.e. observational vector for each person is assumed to form a column of 
//         a matrix _dim x N (N = sample size)
//       * these are either observed data or left or right or lower limit of interval 
//         in the case of censoring
// y_right ....... long vector with upper limits of interval censored observations
//       * ignored if the observation is not interval censored
//       * it does not have to be supplied if there are no interval censored observations
// status ........ long vector with censoring indicators for each observations 
//         (1 = observed, 0 = right censored, 2 = left censored, 3 = interval censored)
// rM[nP] ........ initial component labels for each observation
//         R indeces (1, 2, ...) on input and on return
//         C++ indeces (0, 1, ...) during the computation
// YsM ........... initial sampled latent responses
// iterM ......... index of the initial iteration
// specif ........ type of the specification of the model (1 or 2)
//
// storeP[2] ..... flags whether to store Ys, r
//
void
bayesHistogram(char** dirP,           int* dimsP,
               double* y_left,        double* y_right,      int* status,
               int* rM,               double* YsM,
               int* iterM,
               const int* specif,     int* GsplineI,        double* GsplineD,
               int* nsimulP,          int* storeP,
               const int* mainSimul,  int* errP)
{
  try{
    int out_format[2] = {6, 1};     /** precision and width for output **/

    int i, j;
    GetRNGstate();

    *errP = 0;
    string dir = *dirP;
 
    /** Numbers of iterations etc.  **/
    int niter = nsimulP[0]; 
    int nthin = nsimulP[1];
    int nwrite = nsimulP[2];

    /** Mandatory storing of simulated values **/
    int* storeaP       = storeP;
    int* storeyP       = storeP + 1;
    int* storerP       = storeP + 2;
    
    /** Check data for consistency and find out how many censored observations we have **/
    int* nP = dimsP;
    int* dimP = GsplineI + 0;
    int n_interval = 0;
    int n_censored = 0;
    for (i = 0; i < (*nP) * (*dimP); i++){
      switch (status[i]){
      case 1:
        break;
      case 0:
      case 2:
        n_censored++;
        break;
      case 3:
        n_censored++;
        n_interval++;
        if (y_left[i] > y_right[i]) throw returnR("Incorrect data supplied, lower limit of the interval is higher than the upper limit", 1);
        break;
      default:
        throw returnR("Incorrect status indicator supplied", 1);
      }
    }
    if (!n_censored) *storeyP = 0;
   
    /** Check initial augmented data for consistency **/
    for (i = 0; i < (*nP) * (*dimP); i++){
      switch (status[i]){
      case 1: YsM[i] = y_left[i]; break;
      case 0: if (YsM[i] < y_left[i])                        throw returnR("Inconsistent augmented data supplied", 1); break;
      case 2: if (YsM[i] > y_left[i])                        throw returnR("Inconsistent augmented data supplied", 1); break;
      case 3: if (YsM[i] < y_left[i] || YsM[i] > y_right[i]) throw returnR("Inconsistent augmented data supplied", 1); break;
      }
    }

    /** Storage space for the G-spline **/
                                          /* _null_weight ...... limit for the mixture component weight to be recorded in a file */
    Gspline* g_y = new Gspline;
    *g_y = Gspline(GsplineI, GsplineD);
    const double _null_weight = _null_mass/g_y->total_length();

    /** Array to store numbers of observations belonging to each mixture component  **/
    /** Subtract 1 from each rM (R -> C++ indeces) and check for consistency        **/
    int* mixtureNM = new int[g_y->total_length()];
    if (mixtureNM == NULL) throw returnR("C++ Error: Could not allocate a memory for the working space", 1);
    for (i = 0; i < g_y->total_length(); i++) mixtureNM[i] = 0;
    for (i = 0; i < *nP; i++){
        rM[i]--;
        if (rM[i] < 0 || rM[i] >= g_y->total_length()) throw returnR("Inconsistent initial rM supplied", 1);
        mixtureNM[rM[i]]++;
    }

    /** Open files to store simulated values for writing **/
    std::string iterpath      = dir + "/iteration.sim";   

    std::string sigmapath      = dir + "/gspline.sim";
    std::string lambdapath    = dir + "/lambda.sim";      std::string mixmomentpath  = dir + "/mixmoment.sim";
    std::string mweightpath   = dir + "/mweight.sim";     std::string mlogweightpath = dir + "/mlogweight.sim";
    std::string mmeanpath     = dir + "/mmean.sim";       std::string Ypath          = dir + "/Y.sim";           
    std::string rpath         = dir + "/r.sim";           std::string logposterpath  = dir + "/logposter.sim";

    std::ofstream iterfile;
    std::ofstream  sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile;

    openFile(iterfile, iterpath, 'a');
    openFiles_bayesHistogram(sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, 
                             Yfile, rfile, logposterfile,
			     sigmapath, lambdapath, mixmomentpath, mweightpath, mlogweightpath, mmeanpath, 
                             Ypath, rpath, logposterpath,
                             n_censored, 'a');

//std::string mixtureNMpath  = dir + "/mixtureNM.sim";
//std::ofstream mixtureNMfile;
//mixtureNMfile.open(mixtureNMpath.c_str(), std::ios::out | std::ios::trunc);
//mixtureNMfile.close();
//mixtureNMfile.open(mixtureNMpath.c_str(), std::ios::out | std::ios::app);

    /***** Working arrays and variables                                                             *****/
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
    int l_moments    = *dimP + ((*dimP)*(1+(*dimP)))/2;
    int l_lambda     = g_y->equal_lambda()? 1 : (*dimP);
    int l_iwork      = g_y->dim() * g_y->total_length();
    int l_dwork      = (g_y->total_length() > l_moments ? g_y->total_length() : l_moments);
    int l_log_poster = (g_y->equal_lambda() ? 3 : 2 + g_y->dim());
    int i_logpr      = (g_y->equal_lambda() ? 2 : 1 + g_y->dim());
    int* iwork         = new int[l_iwork];
    double* dwork      = new double[l_dwork];
    double** mu        = new double*[g_y->dim()];
    double* log_poster = new double[l_log_poster];
    if (iwork == NULL || dwork == NULL || mu == NULL || log_poster == NULL) throw returnR("Could not allocate a working memory", 1);
    for (j = 0; j < g_y->dim(); j++){
      mu[j] = new double[g_y->length(j)];
      if (mu[j] == NULL) throw returnR("Could not allocate a working memory", 1);
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
    int nullthIter = *iterM;
    int iter;
    int lastIter = nullthIter + niter;
    int iterTotal = nullthIter*nthin;
    int iterTotalNow = 0;
    int backs = 0;
    int writeAll = 0;
    Rprintf("Iteration ");    
    for (iter=nullthIter + 1; iter <= lastIter; iter++){
      for (int witer = 1; witer <= nthin; witer++){  /**  thinning cycle  **/
        iterTotal++;                        /* = (iter-1)*nthin + witer                    */
        iterTotalNow++;                     /* = (iter-1)*nthin + witer - nullthIter*nthin */

        update_Data_GS(YsM, y_left, y_right, status, rM, g_y, nP, &n_censored);
        g_y->update_alla(mixtureNM, a_ipars, &iterTotalNow);
        switch (*specif){
        case 1:
          g_y->update_gamma(YsM, rM, nP);
          g_y->update_sigma(YsM, rM, nP, &iterTotalNow);
          break;
        case 2:
          g_y->update_Intcpt(YsM, rM, nP);
          g_y->update_Scale(YsM, rM, nP, &iterTotalNow);
          break;
        }
        g_y->update_lambda();
        update_Alloc_GS(rM, mixtureNM, mu, log_poster + 0, log_poster + i_logpr, g_y, YsM, nP, iwork, dwork);

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
        writeToFiles_bayesHistogram(g_y, rM, YsM, log_poster,
	    l_moments, l_lambda, l_log_poster, nP, storeaP, storeyP, storerP, n_censored, &writeAll, iwork, dwork,
            sigmafile, lambdafile, mixmomentfile, mweightfile, mlogweightfile, mmeanfile, Yfile, rfile, logposterfile,
            _null_weight, out_format[0], out_format[1]);
        writeAll = 0;
      }

//writeToFile_1(mixtureNM, g_y->total_length(), mixtureNMfile, 6, 1);
      
    }    /** end of the main cycle over iter **/
    iterfile.close();
    closeFiles_bayesHistogram(sigmafile, lambdafile, mixmomentfile, mweightfile, 
                              mlogweightfile, mmeanfile, Yfile, rfile, logposterfile, n_censored);
    Rprintf("\n");
//mixtureNMfile.close();


    /**  Do some adjustments on last sampled values before returning them back to R  **/
    *iterM = iter - 1;
    for (i = 0; i < *nP; i++) rM[i]++;        // C++ -> R

    /** Write current G-spline values to initial arrays **/
    g_y->Gspline2initArray(GsplineI, GsplineD);

    PutRNGstate();

    /** Cleaning **/
    for (j = 0; j < g_y->dim(); j++) delete [] mu[j];
    delete [] mu;                      
    delete [] log_poster;     delete g_y;    delete [] mixtureNM;
    delete [] iwork;          delete [] dwork;

    return;
  }  /** end of try **/
  catch(returnR rr){
    *errP = rr.errflag();
    PutRNGstate();
    return;
  }
}    /** end of function bayesHistogram **/

}    /** end of extern "C"  **/
