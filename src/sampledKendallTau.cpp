// Compute Kendalls's tau at each iteration of MCMC 
// and return a sample from it
//
// * only for a bivariate version!!! -> to be checked in R, not here
//
// 26/11/2005: 'sampledKendallTau'
//             'evalKendallTau'
//
#include "sampledKendallTau.h"

extern "C" {

// 
// PARAMETERS:
// 
// Tau[M_now] ................... computed values of Tau at each iteration
// M_now ........................ current sample size used here (after taking into account 'skip' and 'by')
// dirP ......................... directory where the sample is stored
// extensP ...................... additional extension by file names (usually "_2" for doubly censored data)
// KK[2] ........................ numbers of knots on each side of the reference knot
// Phi0[(2*KK[0]+1)^2] .......... values of Phi((mu[0,i] - mu[0,j])/(sqrt(2)*sigma0)) (in COLUMN major order)
// Phi1[(2*KK[1]+1)^2] .......... values of Phi((mu[1,i] - mu[1,j])/(sqrt(2)*sigma0)) (in COLUMN major order)
// M ............................ McMC sample size (total, 'skip' and 'by' iterations included)
//                                * M should be <= number of rows in *.sim files
//                                * here: it is an index of the last iteration used to compute the average
// skip ......................... how many rows are to be skipped at the beginning of the sample
// by ........................... only every 'by' G-spline will be taken into account
// nwrite ....................... frequency of informing the user about the progress
// errP ......................... error flag
//
void
sampledKendallTau(double* Tau,           int* M_now,
                  char** dirP,           char** extensP,
                  const int* KK,
                  const double* Phi0,    const double* Phi1,   
                  const int* M,          const int* skip,     const int* by,         const int* nwrite,
                  int* errP)
{
  try{
    double* pTau = Tau;
  
    const int dim = 2;
    const int length0 = 2*KK[0] + 1;
    const int length1 = 2*KK[1] + 1;
    const int total_length = length0 * length1;
    //    bool test = false;
    int i, j, k, l, pp;

    *errP = 0;
    string dir        = *dirP;
    string extens     = *extensP;

    /* Open files with simulated G-splines and skip rows at the beginning of each file that are to be skipped */
    int k_effect;
    double* w                 = (double*)  calloc(total_length, sizeof(double));
    int** ind_mu              = (int**)    calloc(dim, sizeof(int*));
    if (!w || !ind_mu) 
      throw returnR("Not enough memory available in sampledKendallTau (w/ind_mu)", 1);
    for (j = 0; j < dim; j++){
      ind_mu[j] = (int*)    calloc(total_length, sizeof(int));
      if (!ind_mu[j]) throw returnR("Not enough memory available in sampledKendallTau (ind_mu[j])", 1);
    }
    std::string kpath = dir + "/mixmoment" + extens + ".sim";    
    std::string wpath = dir + "/mweight" + extens + ".sim";
    std::string mupath = dir + "/mmean" + extens + ".sim";
    std::ifstream kfile, wfile, mufile;
    openGsplineFiles_forTau(kfile, wfile, mufile, kpath, wpath, mupath, *skip + 1);   /* skip also header */

    /* Rearange Phis to have them as matrices   */
    double** mPhi0 = (double**) calloc(length0, sizeof(double*));
    double** mPhi1 = (double**) calloc(length1, sizeof(double*));
    if (!mPhi0 || !mPhi1) throw returnR("Not enough memory available in sampledKendallTau (mPhi0/mPhi1)", 1);
    for (i = 0; i < length0; i++){
      mPhi0[i] = (double*) calloc(length0, sizeof(double));
      if (!mPhi0[i]) throw returnR("Not enough memory available in sampledKendallTau (mPhi0[i])", 1);
    }
    for (j = 0; j < length1; j++){
      mPhi1[j] = (double*) calloc(length1, sizeof(double));
      if (!mPhi1[j]) throw returnR("Not enough memory available in sampledKendallTau (mPhi1[j])", 1);
    }
    pp = 0;
    for (k = 0; k < length0; k++){
      for (i = 0; i < length0; i++){
        mPhi0[i][k] = Phi0[pp];
        pp++;
      }
    }
    pp = 0;
    for (l = 0; l < length1; l++){
      for (j = 0; j < length1; j++){
        mPhi1[j][l] = Phi1[pp];
        pp++;
      }
    }

    /* Compute Phi((mu[0,i] - mu[0,k])/(sqrt(2)*sigma0)) * Phi((mu[1,j] - mu[1,l])/(sqrt(2)*sigma1))   */
    double**** PhiPhi = (double****) calloc(length0, sizeof(double***));
    if (!PhiPhi) throw returnR("Not enough memory available in sampledKendallTau (PhiPhi)", 1);
    for (i = 0; i < length0; i++){
      PhiPhi[i] = (double***) calloc(length1, sizeof(double**));
      if (!PhiPhi[i]) throw returnR("Not enough memory available in sampledKendallTau (PhiPhi[i])", 1);
      for (j = 0; j < length1; j++){
        PhiPhi[i][j] = (double**) calloc(length0, sizeof(double*));
        if (!PhiPhi[i][j]) throw returnR("Not enough memory available in sampledKendallTau (PhiPhi[i][j])", 1);
        for (k = 0; k < length0; k++){
          PhiPhi[i][j][k] = (double*) calloc(length1, sizeof(double));
          if (!PhiPhi[i][j][k]) throw returnR("Not enough memory available in sampledKendallTau (PhiPhi[i][j][k])", 1);
          for (l = 0; l < length1; l++){
	    PhiPhi[i][j][k][l] = mPhi0[i][k] * mPhi1[j][l];
          }
        }
      }
    }

    /* Loop over McMC iterations */
    if (*skip >= *M) throw returnR("More McMC iterations should be skipped than available", 1);
    readGsplineFromFiles_forTau(&k_effect, w, ind_mu, 0, *skip, dim, KK, total_length, kfile, wfile, mufile, kpath, wpath, mupath);
    evalKendallTau(pTau, &dim, &k_effect, w, ind_mu, PhiPhi);

    *M_now = 1;
    int by_1 = *by - 1;
    int backs = 0;
    Rprintf("Iteration ");
    for (int iter = *skip + 1 + (*by); iter <= *M; iter += (*by)){
      pTau++;
      readGsplineFromFiles_forTau(&k_effect, w, ind_mu, by_1, iter, dim, KK, total_length, kfile, wfile, mufile, kpath, wpath, mupath);
      evalKendallTau(pTau, &dim, &k_effect, w, ind_mu, PhiPhi);

      (*M_now)++;
      if (!(iter % (*nwrite)) || iter == *M){
        for (i = 0; i < backs; i++) Rprintf("\b");
        Rprintf("%d", iter);
        backs = int(log10(double(iter))) + 1;
      }
    }    /** end of the while over iterations **/
    Rprintf("\n");

    /* Close files with simulated G-splines */
    kfile.close();
    wfile.close();
    mufile.close();

    /* Cleaning */
    for (i = 0; i < length0; i++){
      for (j = 0; j < length1; j++){
        for (k = 0; k < length0; k++){
          free(PhiPhi[i][j][k]);
        }
        free(PhiPhi[i][j]);
      }
      free(PhiPhi[i]);
    }
    free(PhiPhi);

    for (j = 0; j < dim; j++){
      free(ind_mu[j]);
    }    
    free(ind_mu);
    free(w);


    return;
  }
  catch(returnR rr){
    *errP = rr.errflag();
    return;
  }
}  /** end of function 'sampledKendallTau'  **/


}  /** end of extern "C" **/


// ======================================================================================
// ******* evalKendallTau
// ======================================================================================
void
evalKendallTau(double* value,  const int* dim,  const int* k_effect,  const double* w,  int** ind_mu,  double**** PhiPhi)
{
  static int pp, qq;
  static int* ip;
  static int* jp;
  static int* kp; 
  static int* lp;
  static const double* wp;
  static const double* wq;
  static double w_w;

  if (*dim != 2) throw returnR("Function 'evalKendallTau' implemented only for dim = 2", 1);

  *value = 0.0;
  ip = ind_mu[0];
  jp = ind_mu[1]; 
  wp = w;
  for (pp = 0; pp < *k_effect; pp++){
    w_w = (*wp) * (*wp);
    *value += w_w * PhiPhi[*ip][*jp][*ip][*jp];
    kp = ip + 1;
    lp = jp + 1;
    wq = wp + 1;
    for (qq = pp+1; qq < *k_effect; qq++){
      w_w = (*wp) * (*wq);
      *value += w_w * PhiPhi[*ip][*jp][*kp][*lp];
      *value += w_w * PhiPhi[*kp][*lp][*ip][*jp];
      kp++;
      lp++;
      wq++;
    }
    ip++;
    jp++;
    wp++;
  }
  *value *= 4;
  *value -= 1;

  return;
}

