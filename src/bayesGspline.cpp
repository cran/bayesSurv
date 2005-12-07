// Compute McMC average of sampled G-splines,
//   possibly return also values of each G-spline
//   

// 01/11/2004: 'bayesGspline'
// 04/11/2004: 'evalGspline'
//             'printReadGspline'
// 13/12/2004: addapted to allow general mean and scale of the G-spline (2nd specification)
// 22/01/2005: added possibility for file names with additional extension (used for doubly-censored data)
// 02/02/2005: added possibility to adjust the G-spline intercept (version 30 and 31)
//             * useful for evaluating fitted G-spline when there is G-spline random intercept in the model
// 16/05/2005: computation of standardized G-splines allowed
//
#include "bayesGspline.h"

extern "C" {

using namespace std;

// 
// PARAMETERS:
// 
// average[nx1 OR nx1*nx2] ........ McMC average of the sampled G-spline evaluated in a (rectangular) grid of values 
//                         * sorted in the folowing way:
//                           if (nx2 == 0): eval = {g(x1[0]), ..., g(x1[nx1-1])}
//                           if (nx2 > 0):  eval = {g(x1[0], x2[0]), ..., g(x1[nx1-1], x2[0]), 
//                                                  g(x1[0], x2[1]), ..., g(x1[nx1-1], x2[1]), 
//                                                   ...................................
//                                                  g(x1[0], x2[nx2-1]), ..., g(x1[nx1-1], x2[nx2-1])}
// value[ngrid * (M-skip)/by] ... values of the G-spline evaluated at each iteration
//                                * if onlyAver == true, it is sufficient to supply an array of length ngrid
// M_now ........................ current sample size used to get average here (after taking into account 'skip' and 'by')
// onlyAver ..................... 1/0 to indicate whether only McMC average is to be returned or all values as well
// dirP ......................... directory where the sample is stored
// extensP ...................... additional extension by file names (usually "_2" for doubly censored data)
// extens_adjP .................. additional extension by file name where to take adjustment intercept
//                                (if version == 30 or 31), usually "_b" or "_b2"
// x1[nx1] ...................... grid for the first dimension
// x2[nx2] ...................... grid for the second dimension 
// total_length ................. total_length of the G-spline (number of components)
// M ............................ McMC sample size (total, 'skip' and 'by' iterations included)
//                                * M should be <= number of rows in *.sim files
//                                * here: it is an index of the last iteration used to compute the average
// skip ......................... how many rows are to be skipped at the beginning of the sample
// by ........................... only every 'by' G-spline will be taken into account
// nwrite ....................... frequency of informing the user about the progress
// nx1 .......................... length of the grid for the first dimension
// nx2 .......................... length of the grid for the second dimension
// version ...................... 0 if used after running 'bayesHistogram' or 'bayesBisurvreg' or 'bayessurvreg2'
//                                30 if used after running 'bayessurvreg3' and error density is to be computed
//                                31 if used after running 'bayessurvreg3' and density of the random intercept is to be computed
//                                => intercept of the G-spline has to be adjusted then
//                   for error:         new intercept = sampled intercept + mean(random intcpt G-spline)
//                   for random intcpt: new intercept = 0 - mean(random intcpt G-spline)
//                                                      (0 = sampled (fixed) intercept)
// standard ..................... 0/1, 1 if I want to compute standardized (zero mean, unit variance) density
//                                currently implemented only for version >= 30 (i.e. for univariate bayessurvreg3)
// errP ......................... error flag
//
void
bayesGspline(double* average,          double* value,   int* M_now,          const int* onlyAver,
             char** dirP,              char** extensP,  char** extens_adjP,  
             double* x1,               double* x2,   
             const int* total_length,  const int* M,    const int* skip,     const int* by,         const int* nwrite,
             const int* nx1,           int* nx2,        const int* version,  int* standard,
             int* errP)
{
  try{
    bool test = false;
    int ix, i, j;
    double* pvalue;

    *errP = 0;
    string dir        = *dirP;
    string extens     = *extensP;
    string extens_adj = *extens_adjP;
    int dim   = (*nx2 > 0 ? 2 : 1);
    int ngrid = (*nx2 > 0 ? (*nx1)*(*nx2) : (*nx1));
    if (*nx2 == 0) *nx2 = 1;
    double* x[2];
    switch (dim){
    case 1: x[0] = x1; break;
    case 2: x[0] = x1; x[1] = x2; break;
    default: throw returnR("bayesGspline function not implemented for dim > 2", 1);
    }
    if (*version < 30) *standard = 0;

    /* Open files with simulated G-splines and skip rows at the beginning of each file that are to be skipped */
    int k_effect;
    double* w                 = (double*)  calloc(*total_length, sizeof(double));
    int** ind_mu              = (int**)    calloc(dim, sizeof(int*));
    double** mu               = (double**) calloc(dim, sizeof(double*));
    double* sigma             = (double*)  calloc(dim, sizeof(double));
    double* gamma             = (double*)  calloc(dim, sizeof(double));
    double* delta             = (double*)  calloc(dim, sizeof(double));
    double* intcpt            = (double*)  calloc(dim, sizeof(double));
    double* scale             = (double*)  calloc(dim, sizeof(double));
    double* min_half_inv_sig2 = (double*)  calloc(dim, sizeof(double));
    if (!w || !ind_mu || !mu || !sigma || !gamma || !delta || !min_half_inv_sig2 || !intcpt || !scale) 
      throw returnR("Not enough memory available in bayesGspline (w/ind_mu/mu/sigma/gamma/delta/intcpt/scale/min_half_inv_sigma2)", 1);
    for (j = 0; j < dim; j++){
      ind_mu[j] = (int*)    calloc(*total_length, sizeof(int));
      mu[j]     = (double*) calloc(*total_length, sizeof(double));
      if (!ind_mu[j] || !mu[j]) throw returnR("Not enough memory available in bayesGspline (ind_mu[j]/mu[j])", 1);
    }
    std::string kpath = dir + "/mixmoment" + extens + ".sim";    
    std::string wpath = dir + "/mweight" + extens + ".sim";
    std::string mupath = dir + "/mmean" + extens + ".sim";       
    std::string sigmapath = dir + "/gspline" + extens + ".sim";
    std::ifstream kfile, wfile, mufile, sigmafile;
    openGsplineFiles(kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath, *skip + 1);   /* skip also header */

    /* Open file with adjustment intercept (for version 30 or 31) */
    std::string int_adjpath = dir + "/mixmoment" + extens_adj + ".sim";
    std::ifstream int_adjfile;
    if (*version >= 30) open_File_toRead(int_adjfile, int_adjpath, *skip + 1);

    /* Open file with mean and variance of the current g-spline (for version 30 or 31) */
    std::string mixmomentpath = dir + "/mixmoment" + extens + ".sim";
    std::ifstream mixmomentfile;
    double E_gx = 0;
    double sd_gx = 1;
    if (*standard) open_File_toRead(mixmomentfile, mixmomentpath, *skip + 1);

    /* Reset averages */
    for (ix = 0; ix < ngrid; ix++){
      average[ix] = 0.0;
    }

    /* Loop over McMC iterations */
    if (*skip >= *M) throw returnR("More McMC iterations should be skipped than available", 1);
    pvalue = value;
    readGsplineFromFiles(&k_effect, w, ind_mu, mu, gamma, sigma, delta, intcpt, scale, 0, *skip, dim, *total_length, 
                         kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
    if (*standard){
      readMean_and_Scale(&E_gx, &sd_gx, 0, *skip, dim, mixmomentfile, mixmomentpath);
    }
    if (*version >= 30){
      adjust_intercept(intcpt, version, &E_gx, 0, *skip, int_adjfile, int_adjpath);
    }
    evalGspline(average, pvalue, nx1, nx2, x, &dim, &k_effect, w, mu, intcpt, sigma, scale, min_half_inv_sig2, standard, &E_gx, &sd_gx);
    *M_now = 1;
    if (test) printReadGspline(*skip + 1, dim, k_effect, w, mu, intcpt, sigma, scale);

    int by_1 = *by - 1;
    int jump_value = (*onlyAver ? 0 : ngrid);
    int backs = 0;
    Rprintf("Iteration ");
    for (int iter = *skip + 1 + (*by); iter <= *M; iter += (*by)){      
      pvalue += jump_value;
      readGsplineFromFiles(&k_effect, w, ind_mu, mu, gamma, sigma, delta, intcpt, scale, by_1, iter, dim, *total_length, 
                           kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
      if (*standard){
        readMean_and_Scale(&E_gx, &sd_gx, by_1, iter, dim, mixmomentfile, mixmomentpath);
      }
      if (*version >= 30){
        adjust_intercept(intcpt, version, &E_gx, by_1, iter, int_adjfile, int_adjpath);
      }
      evalGspline(average, pvalue, nx1, nx2, x, &dim, &k_effect, w, mu, intcpt, sigma, scale, min_half_inv_sig2, standard, &E_gx, &sd_gx);
      (*M_now)++;
      if (test) printReadGspline(iter, dim, k_effect, w, mu, intcpt, sigma, scale);

      if (!(iter % (*nwrite)) || iter == *M){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
      }
    }    /** end of the while over iterations **/
    Rprintf("\n");

    /* Close files with simulated G-splines */
    closeGsplineFiles(kfile, wfile, mufile, sigmafile);
    if (*version >= 30) int_adjfile.close();
    if (*standard) mixmomentfile.close();

    /* McMC averages */
    for (ix = 0; ix < ngrid; ix++){
      average[ix] /= (*M_now);
    }

    /* Cleaning */
    for (j = 0; j < dim; j++){
      free(ind_mu[j]);   free(mu[j]);
    }    
    free(ind_mu);    free(mu);
    free(w);         free(sigma);   free(gamma);   free(delta);   free(min_half_inv_sig2);
    free(intcpt);    free(scale);
    
    return;
  }
  catch(returnR rr){
    *errP = rr.errflag();
    return;
  }
}

}   /** end of extern "C" **/


// ======================================================================================
// ******* evalGspline
// ======================================================================================
//
// average[ngrid] ......... on INPUT: something,
//                          on OUTPUT: something + value[ngrid]
//                          !!! to get average, this must be at the end divided by the sample size !!!
// value[ngrid] ........... on OUTPUT: G-spline evaluated in a grid given by 'x' 
//
void
evalGspline(double* average,      double* value,  
            const int* nx1,       const int* nx2,       double** x,
            const int* dim,       const int* k_effect,  
            const double* w,      double** mu,          const double* intcpt,  
            const double* sigma,  const double* scale,  double* min_half_inv_sig2,
            const int* standard,  const double* E_gx,   const double* sd_gx)
{
  //  static int ix;
  static int j, k, col, row;
  static double inv_sigma, inv_sigmas, inv_scale, x_min_mu, val_k_ix, tmp;
  static double inv_sd_gx;

  /** Compute some transformations of sigma's **/
  inv_sigmas = 1.0;
  for (j = 0; j < *dim; j++){
    inv_sigma = 1/sigma[j];
    inv_scale = 1/scale[j];
    tmp = inv_sigma * inv_scale;
    inv_sigmas *= invSQRT_TWO_PI * tmp;              // 1/(sqrt(2*pi) * sigma * scale)
    min_half_inv_sig2[j] = -0.5 * tmp * tmp;
  }

  /** Compute G-spline values **/
  double* val_ix = value;
  double* aver_ix = average;
  if (*standard){
    if (*dim > 1) throw returnR("C++ Error: evalGspline not implemented for dim > 1", 1);

    inv_sigmas *= (*sd_gx);                        // sd_gx/(sqrt(2*pi) * sigma * scale)
    min_half_inv_sig2[0] *= (*sd_gx) * (*sd_gx);   // -0.5 * sd_gx^2 /(sigma^2 * scale^2)
    inv_sd_gx = 1/(*sd_gx);

    for (row = 0; row < *nx1; row++){
      *val_ix = 0.0;
      for (k = 0; k < *k_effect; k++){
        val_k_ix = 0.0;
        x_min_mu = x[0][row] + inv_sd_gx*((*E_gx) -  intcpt[0] - scale[0]*mu[0][k]);
        val_k_ix += min_half_inv_sig2[0] * x_min_mu * x_min_mu;
        val_k_ix = (val_k_ix < -3*_emax ? 0.0 : w[k]*exp(val_k_ix));
        *val_ix += val_k_ix;
      }
      *val_ix *= inv_sigmas;
      *aver_ix += (*val_ix);
      val_ix++;
      aver_ix++;
    }
  }
  else{
    for (col = 0; col < *nx2; col++){
      for (row = 0; row < *nx1; row++){
  //      ix = col*(*nx1) + row;
        *val_ix = 0.0;
        for (k = 0; k < *k_effect; k++){
          val_k_ix = 0.0;
          for (j = 0; j < *dim; j++){
            x_min_mu = x[j][j==0 ? row : col] - intcpt[j] - scale[j]*mu[j][k];
            val_k_ix += min_half_inv_sig2[j] * x_min_mu * x_min_mu;
          }            
          val_k_ix = (val_k_ix < -3*_emax ? 0.0 : w[k]*exp(val_k_ix));
          *val_ix += val_k_ix;
        }
        *val_ix *= inv_sigmas;
        *aver_ix += (*val_ix);
        val_ix++;
        aver_ix++;
      }
    }
  }

  return;
}


// ======================================================================================
// ******* printReadGspline
// ======================================================================================
void
printReadGspline(const int& iter,       const int& dim,       const int& k_effect,  const double* w,  double** mu,  
                 const double* intcpt,  const double* sigma,  const double* scale)
{
  int i, j;

  Rprintf("G-spline %d: ", iter);
  Rprintf("  k = %d,\n", k_effect);
  Rprintf("   w = "); for (i = 0; i < k_effect; i++) Rprintf("%f, ", w[i]); Rprintf("\n");
  for (j = 0; j < dim; j++){
    Rprintf("   mu%d = ", j+1); for (i = 0; i < k_effect; i++) Rprintf("%f, ", mu[j][i]); Rprintf("\n");
  }
  Rprintf("        sigma = "); for (j = 0; j < dim; j++) Rprintf("%f, ", sigma[j]); Rprintf("\n");
  Rprintf("    intercept = "); for (j = 0; j < dim; j++) Rprintf("%f, ", intcpt[j]); Rprintf("\n");
  Rprintf("        scale = "); for (j = 0; j < dim; j++) Rprintf("%f, ", scale[j]); Rprintf("\n");
  Rprintf("**************************************\n");
  return;
}
