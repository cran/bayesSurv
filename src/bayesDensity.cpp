// Compute McMC averages of sampled densities
//  * both conditional on k and overall
//  * both standardized and unstandardized

// 04/03/2004: start working on it
//

#include "bayessurvreg.h"

extern "C" {

using namespace std;

//
// PARAMETERS:
//
// aver .......... McMC averages of densities, for each k in one column after returning back to R
//                                             k = 0 <=> unconditional density                         (ngrid x (1 + kmax))
// staver ........ McMC averages of standardized densities                                             (nstgrid x (1 + kmax))
// intercept ..... sample of intercepts (computed from the sampled densities)                          (M)
// scale ......... sample of scales (computed from the sampled densities)                              (M)
// Mk ............ numbers of densities with given number of components,                               (1 + kmax)
//                                             Mk[0] := M                           
// dirP .......... directory where the sample is stored
// grid .......... grid of values where to evaluate unstandardized density                             (ngrid)
// stgrid ........ grid of values where to evaluate standardized density                               (nstgrid)
// kmax .......... maximal number of components given to the McMC                                      (1)
// M ............. sample size                                                                         (1)
// skip .......... how many rows are to be skipped at the beginning of the sample                      (1)
// ngrid ......... length of the grid                                                                  (1)
// nstgrid ....... length of the stgrid                                                                (1)
// type .......... what to compute: 0 = both standardized and unstandardized
//                                  1 = only unstandardized
//                                  2 = only standardized
//
void
bayesDensity(double* aver,      double* staver,      double* intercept,     double* scale,      int* Mk,
             char** dirP,       const double* grid,  const double* stgrid,
             const int* kmax,   const int* M,        const int* skip,
             const int* ngrid,  const int* nstgrid,
             const int* type,   int* errP)
{
  try{
    *errP = 0;
    string dir = *dirP;

    // lrow ..... length of one mixture in the array
    int lrow = 3*(*kmax) + 1;           

    int i, j;

    // Reset averages and numbers of sampled components
    for (j = 0; j <= *kmax; j++){
      Mk[j] = 0;
      if (*type <= 1)               for (i = 0; i < *ngrid; i++)   aver[j*(*ngrid) + i] = 0.0;
      if (*type == 0 || *type == 2) for (i = 0; i < *nstgrid; i++) staver[j*(*nstgrid) + i] = 0.0;
    }
    Mk[0] = *M;

    // Read sampled densities
    double* sample = new double[(*M)*(1 + 3*(*kmax))];
    if (sample == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    readMixtureFromFiles(sample, *M, *kmax, dir, "/mixmoment.sim", "/mweight.sim", "/mmean.sim", "/mvariance.sim", 1 + (*skip));    

    // Compute
    double* w;     double* mu;     double* sig2;
    double* sig = new double[*kmax];
    if (sig == NULL) throw returnR("C++ Error: Could not allocate a memory for a working space, buy more memory...", 1);
    double f, stmu, stsig;
    int iter, k;
    Rprintf("Computing predictive densities. \n");
    for (iter = 0; iter < *M - (*skip); iter++){

      // Sampled values
      k = int(sample[iter*lrow]);
      w = sample + iter*lrow + 1;
      mu = sample + iter*lrow + 1 + (*kmax);
      sig2 = sample + iter*lrow + 1 + 2*(*kmax);
      Mk[k]++;

      // Compute intercept and scale
      intercept[iter] = 0.0;
      scale[iter] = 0.0;
      for (j = 0; j < k; j++){
        sig[j] = sqrt(sig2[j]);
        intercept[iter] += w[j]*mu[j];
        scale[iter] += w[j]*(mu[j]*mu[j] + sig2[j]);
      }
      scale[iter] -= intercept[iter]*intercept[iter];
      if (scale[iter] <= 0.0) scale[iter] = 0.0;
      scale[iter] = sqrt(scale[iter]);

      // Compute values of a sampled density in grid and add it to cumsum densities
      if (*type <= 1){
        for (i = 0; i < *ngrid; i++){
          f = 0.0;
          for (j = 0; j < k; j++){
            f += w[j] * dnorm(grid[i], mu[j], sig[j], 0);
          }
          aver[k*(*ngrid) + i] += f;    // conditional density
          aver[i] += f;                 // overall density
        }
      }

      if (*type == 0 || *type == 2){
        for (i = 0; i < *nstgrid; i++){
          f = 0.0;
          for (j = 0; j < k; j++){
            stmu = (mu[j] - intercept[iter])/scale[iter];
            stsig = sig[j]/scale[iter];
            f += w[j] * dnorm(stgrid[i], stmu, stsig, 0);
          }
          staver[k*(*nstgrid) + i] += f;    // conditional density
          staver[i] += f;                   // overall density
        }
      }

    }    // end of the loop over iterations

    // Compute averages
    for (j = 0; j <= *kmax; j++){
      if (Mk[j] == 0) continue;
      if (*type <= 1)               for (i = 0; i < *ngrid; i++)   aver[j*(*ngrid) + i] /= Mk[j];
      if (*type == 0 || *type == 2) for (i = 0; i < *nstgrid; i++) staver[j*(*nstgrid) + i] /= Mk[j];
    }

    delete [] sig;
    delete [] sample;
    return;
  }  // end of try
  catch(returnR rr){
    *errP = rr.errflag();
    return;
  }
}    // end of function bayesDensity

}    // end of extern "C"
