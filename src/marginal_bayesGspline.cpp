// Compute McMC average of sampled marginal G-splines (in the bivariate case),
//   possibly return also values of each marginal G-spline
//   
// 27/11/2005: 'marginal_bayesGspline'
//             'marginal_evalGspline'             
//
#include "marginal_bayesGspline.h"

extern "C" {

using namespace std;

// 
// PARAMETERS:
// 
// average1[nx1] ................ McMC average of the sampled marginal G-spline (margin=1) evaluated in a grid of values 
// average2[nx2] ................ McMC average of the sampled marginal G-spline (margin=2) evaluated in a grid of values
// value1[nx1 * (M-skip)/by] .... values of the marginal G-spline (margin=1) evaluated at each iteration
//                                * if onlyAver == true, it is sufficient to supply an array of length nx1
// value2[nx2 * (M-skip)/by] .... values of the marginal G-spline (margin=2) evaluated at each iteration
//                                * if onlyAver == true, it is sufficient to supply an array of length nx2
// M_now ........................ current sample size used to get average here (after taking into account 'skip' and 'by')
// onlyAver ..................... 1/0 to indicate whether only McMC average is to be returned or all values as well
// dirP ......................... directory where the sample is stored
// extensP ...................... additional extension by file names (usually "_2" for doubly censored data)
// x1[nx1] ...................... grid for the first dimension
// x2[nx2] ...................... grid for the second dimension 
// KK[2] ........................ numbers of knots on each side of the reference knot
// M ............................ McMC sample size (total, 'skip' and 'by' iterations included)
//                                * M should be <= number of rows in *.sim files
//                                * here: it is an index of the last iteration used to compute the average
// skip ......................... how many rows are to be skipped at the beginning of the sample
// by ........................... only every 'by' G-spline will be taken into account
// nwrite ....................... frequency of informing the user about the progress
// nx1 .......................... length of the grid for the first margin
// nx2 .......................... length of the grid for the second margin
// errP ......................... error flag
//
void
marginal_bayesGspline
   (double* average1,   double* average2,      double* value1,       double* value2,
    int* M_now,         const int* onlyAver,   char** dirP,          char** extensP,
    const double* x1,   const double* x2,
    const int* KK,      
    const int* M,       const int* skip,       const int* by,        const int* nwrite,
    const int* nx1,     const int* nx2,        int* errP)
{
  try{
    const int dim = 2;
    const int length0 = 2*KK[0] + 1;
    const int length1 = 2*KK[1] + 1;
    const int total_length = length0 * length1;
    const int standard = 0;
    double* pvalue1;
    double* pvalue2;
    int i, j, ix;

    *errP = 0;
    string dir        = *dirP;
    string extens     = *extensP;

    /* Open files with simulated G-splines and skip rows at the beginning of each file that are to be skipped */
    int k_effect;
    double* w_temp            = (double*)  calloc(total_length, sizeof(double));
    double** w                = (double**) calloc(dim, sizeof(double*));
    double** mu               = (double**) calloc(dim, sizeof(double*));
    double* sigma             = (double*)  calloc(dim, sizeof(double));
    double* gamma             = (double*)  calloc(dim, sizeof(double));
    double* delta             = (double*)  calloc(dim, sizeof(double));
    double* intcpt            = (double*)  calloc(dim, sizeof(double));
    double* scale             = (double*)  calloc(dim, sizeof(double));
    double* inv_sigmas        = (double*)  calloc(dim, sizeof(double));
    double* min_half_inv_sig2 = (double*)  calloc(dim, sizeof(double));
    if (!w_temp || !w || !mu || !sigma || !gamma || !delta || !inv_sigmas || !min_half_inv_sig2 || !intcpt || !scale) 
      throw returnR("Not enough memory available in marginal_bayesGspline (w_temp/w/mu/sigma/gamma/delta/intcpt/scale/inv_sigmas/min_half_inv_sigma2)", 1);
    w[0] = (double*) calloc(length0, sizeof(double));
    w[1] = (double*) calloc(length1, sizeof(double));
    mu[0] = (double*) calloc(length0, sizeof(double));
    mu[1] = (double*) calloc(length1, sizeof(double));
    if (!w[0] || !w[1] || !mu[0] || !mu[1]) throw returnR("Not enough memory available in marginal_bayesGspline (w[j]/mu[j])", 1);
    std::string kpath = dir + "/mixmoment" + extens + ".sim";    
    std::string wpath = dir + "/mweight" + extens + ".sim";
    std::string mupath = dir + "/mmean" + extens + ".sim";       
    std::string sigmapath = dir + "/gspline" + extens + ".sim";
    std::ifstream kfile, wfile, mufile, sigmafile;
    openGsplineFiles(kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath, *skip + 1);   /* skip also header */

    /* Reset averages */
    for (ix = 0; ix < *nx1; ix++) average1[ix] = 0.0;
    for (ix = 0; ix < *nx2; ix++) average2[ix] = 0.0;

    /* Loop over McMC iterations */
    if (*skip >= *M) throw returnR("More McMC iterations should be skipped than available", 1);
    pvalue1 = value1;
    pvalue2 = value2;
    readGsplineFromFiles_forMarginal(w_temp, w, mu, gamma, sigma, delta, intcpt, scale, KK, 0, *skip, total_length, 
                                     kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
    marginal_evalGspline(average1, average2, pvalue1, pvalue2, &length0, &length1, nx1, nx2, x1, x2, 
                         w, mu, intcpt, sigma, scale, inv_sigmas, min_half_inv_sig2);
    *M_now = 1;

    int by_1 = *by - 1;
    int jump_value1 = (*onlyAver ? 0 : (*nx1));
    int jump_value2 = (*onlyAver ? 0 : (*nx2));
    int backs = 0;
    Rprintf("Iteration ");
    for (int iter = *skip + 1 + (*by); iter <= *M; iter += (*by)){
      pvalue1 += jump_value1;
      pvalue2 += jump_value2;
      readGsplineFromFiles_forMarginal(w_temp, w, mu, gamma, sigma, delta, intcpt, scale, KK, by_1, iter, total_length, 
                                       kfile, wfile, mufile, sigmafile, kpath, wpath, mupath, sigmapath);
      marginal_evalGspline(average1, average2, pvalue1, pvalue2, &length0, &length1, nx1, nx2, x1, x2, 
                           w, mu, intcpt, sigma, scale, inv_sigmas, min_half_inv_sig2);
      (*M_now)++;

      if (!(iter % (*nwrite)) || iter == *M){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
      }
    }    /** end of the while over iterations **/
    Rprintf("\n");

    /* Close files with simulated G-splines */
    closeGsplineFiles(kfile, wfile, mufile, sigmafile);

    /* McMC averages */
    for (ix = 0; ix < *nx1; ix++) average1[ix] /= (*M_now);
    for (ix = 0; ix < *nx2; ix++) average2[ix] /= (*M_now);

    /* Cleaning */
    for (j = 0; j < dim; j++){
      free(w[j]);
      free(mu[j]);
    }    
    free(w);                  free(mu);
    free(w_temp);             free(sigma);   
    free(gamma);              free(delta);   
    free(min_half_inv_sig2);  free(inv_sigmas);
    free(intcpt);             free(scale);
    
    return;
  }
  catch(returnR rr){
    *errP = rr.errflag();
    return;
  }
}

}  /** end of extern "C"  **/

// ======================================================================================
// ******* marginal_evalGspline
// ======================================================================================
//
// average1[nx1] ......... on INPUT: something,
//                         on OUTPUT: something + value1[nx1]
//                          !!! to get average, this must be at the end divided by the sample size !!!
// average2[nx2] ......... on INPUT: something,
//                         on OUTPUT: something + value2[nx2]
//                          !!! to get average, this must be at the end divided by the sample size !!!
// value1[nx1] ........... on OUTPUT: marginal G-spline (margin=1) evaluated in a grid given by 'x1' 
// value2[nx2] ........... on OUTPUT: marginal G-spline (margin=2) evaluated in a grid given by 'x2' 
// nx1 ................... length of the grid for the first margin
// nx2 ................... length of the grid for the second margin
// x1[nx1] ............... grid for the first margin
// x2[nx2] ............... grid for the second margin
//
void
marginal_evalGspline
   (double* average1,     double* average2,     double* value1,        double* value2,  
    const int* length0,   const int* length1,
    const int* nx1,       const int* nx2,       const double* x1,      const double* x2,
    double** w,           double** mu,          const double* intcpt,  
    const double* sigma,  const double* scale,  double* inv_sigmas,    double* min_half_inv_sig2)
{
  const int dim = 2;
  static int j, k, ii, dd;
  static double inv_sigma, inv_scale, x_min_mu, val_k_ix, tmp;
  double* val_ix;
  double* aver_ix;
  const double* x_ix;

  /** Compute some transformations of sigma's **/
  for (j = 0; j < dim; j++){
    inv_sigma = 1/sigma[j];
    inv_scale = 1/scale[j];
    tmp = inv_sigma * inv_scale;
    inv_sigmas[j] = invSQRT_TWO_PI * tmp;              // 1/(sqrt(2*pi) * sigma * scale)
    min_half_inv_sig2[j] = -0.5 * tmp * tmp;
  }

  /** Compute G-spline values for margin 1 **/
  dd = 0;
  val_ix = value1;
  aver_ix = average1;
  x_ix = x1;
  for (ii = 0; ii < *nx1; ii++){
    *val_ix = 0.0;    
    for (k = 0; k < *length0; k++){
      x_min_mu = *x_ix - intcpt[dd] - scale[dd]*mu[dd][k];
      val_k_ix = min_half_inv_sig2[dd] * x_min_mu * x_min_mu;
      val_k_ix = (val_k_ix < -3*_emax ? 0.0 : w[dd][k]*exp(val_k_ix));
      *val_ix += val_k_ix;
    }
    *val_ix *= inv_sigmas[dd];
    *aver_ix += (*val_ix);
    val_ix++;
    aver_ix++;
    x_ix++;
  }

  /** Compute G-spline values for margin 2 **/
  dd = 1;
  val_ix = value2;
  aver_ix = average2;
  x_ix = x2;
  for (ii = 0; ii < *nx2; ii++){
    *val_ix = 0.0;    
    for (k = 0; k < *length1; k++){
      x_min_mu = *x_ix - intcpt[dd] - scale[dd]*mu[dd][k];
      val_k_ix = min_half_inv_sig2[dd] * x_min_mu * x_min_mu;
      val_k_ix = (val_k_ix < -3*_emax ? 0.0 : w[dd][k]*exp(val_k_ix));
      *val_ix += val_k_ix;
    }
    *val_ix *= inv_sigmas[dd];
    *aver_ix += (*val_ix);
    val_ix++;
    aver_ix++;
    x_ix++;
  }

  return;
}


