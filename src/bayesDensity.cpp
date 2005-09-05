// Compute McMC averages of sampled densities
//  * both conditional on k and overall
//  * both standardized and unstandardized

// 04/03/2004: start working on it
//

#include "bayesDensity.h"

extern "C" {

using namespace std;

//
// PARAMETERS:
//
// aver .......... McMC averages of densities, for each k in one column after returning back to R
//                                             k = 0 <=> unconditional density                         (ngrid x (1 + kmax))
// staver ........ McMC averages of standardized densities                                             (nstgrid x (1 + kmax))
// centaver ...... McMC averages of centered densities
// intercept ..... sample of intercepts (computed from the sampled densities)                          (M)
// scale ......... sample of scales (computed from the sampled densities)                              (M)
// Mk ............ numbers of densities with given number of components,                               (1 + kmax)
//                                             Mk[0] := total sample size used
// dirP .......... directory where the sample is stored
// grid .......... grid of values where to evaluate unstandardized density                             (ngrid)
// stgrid ........ grid of values where to evaluate standardized density                               (nstgrid)
// kmax .......... maximal number of components given to the McMC                                      (1)
// M ............. available McMC sample size                                                          (1)
// skip .......... how many rows are to be skipped at the beginning of the sample                      (1)
// by ............ possible additional thinnning
// ngrid ......... length of the grid (set to 0 if do not want to compute it)                          (1)
// nstgrid ....... length of the stgrid (set to 0 if do not want to compute it)                        (1)
// ncentgrid ..... length of the centgrid (set to 0 if do not want to compute it)                      (1)
//
void
bayesDensity(double* aver,        double* staver,        double* centaver,
             double* intercept,   double* scale,         int* Mk,
             char** dirP,       
             const double* grid,  const double* stgrid,  const double* centgrid,
             const int* kmax,     const int* M,          const int* skip,         const int* by,
             const int* ngrid,    const int* nstgrid,    const int* ncentgrid,
             int* errP)
{
  try{
    *errP = 0;
    string dir = *dirP;

    // lrow ..... length of one mixture in the array
    int lrow = 3*(*kmax) + 1;           

    int i, j;

    // Reset averages and numbers of sampled components
    if (*skip >= *M) throw returnR("C++ Error: Too many iterations are to be skipped", 1);
    if (*by <= 0) throw returnR("C++ Error: by parameter must be positive", 1);
    if (*skip < 0) throw returnR("C++ Error: skip parameter must not be negative", 1);
    for (j = 0; j <= *kmax; j++){
      Mk[j] = 0;
      for (i = 0; i < *ngrid; i++)   aver[j*(*ngrid) + i] = 0.0;
      for (i = 0; i < *nstgrid; i++) staver[j*(*nstgrid) + i] = 0.0;
      for (i = 0; i < *ncentgrid; i++) centaver[j*(*ncentgrid) + i] = 0.0;
    }
    Mk[0] = 1 + (*M - (*skip) - 1) / (*by);

    // Read sampled densities
    int nread;
    double* sample = (double*) calloc(Mk[0]*lrow, sizeof(double));
    if (!sample) throw returnR("C++ Error: Could not allocate a memory for a working space 'sample'", 1);
    readMixtureFromFiles(sample, &nread, *M, *skip, *by, *kmax, dir, "/mixmoment.sim", "/mweight.sim", "/mmean.sim", "/mvariance.sim");
    if (nread != Mk[0]) throw returnR("Different MCMC sample sizes indicated by mixture files", 1);

    // Compute
    double* sig = (double*) calloc(*kmax, sizeof(double));
    if (!sig) throw returnR("C++ Error: Could not allocate a memory for a working space 'sig'", 1);
    double f, stmu, stsig, centmu;
    int iter, k;
    double* proef = sample;
    double* w     = proef + 1;     
    double* mu    = w + (*kmax);     
    double* sig2  = mu + (*kmax);
    Rprintf("Computing predictive densities. \n");
    for (iter = 0; iter < Mk[0]; iter++){

      // Sampled values
      k = int(*proef);
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
      for (i = 0; i < *ngrid; i++){
        f = 0.0;
        for (j = 0; j < k; j++){
          f += w[j] * dnorm(grid[i], mu[j], sig[j], 0);
        }
        aver[k*(*ngrid) + i] += f;    // conditional density
        aver[i] += f;                 // overall density
      }

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

      for (i = 0; i < *ncentgrid; i++){
        f = 0.0;
        for (j = 0; j < k; j++){
          centmu = (mu[j] - intercept[iter]);
          f += w[j] * dnorm(centgrid[i], centmu, sig[j], 0);
        }
        centaver[k*(*ncentgrid) + i] += f;    // conditional density
        centaver[i] += f;                     // overall density
      }

      proef += lrow;
      w     = proef + 1;     
      mu    = w + (*kmax);     
      sig2  = mu + (*kmax);
    }    // end of the loop over iterations

    // Compute averages
    for (j = 0; j <= *kmax; j++){
      if (Mk[j] == 0) continue;
      for (i = 0; i < *ngrid; i++)   aver[j*(*ngrid) + i] /= Mk[j];
      for (i = 0; i < *nstgrid; i++) staver[j*(*nstgrid) + i] /= Mk[j];
      for (i = 0; i < *ncentgrid; i++) centaver[j*(*ncentgrid) + i] /= Mk[j];
    }

    free(sig);
    free(sample);
    return;
  }  // end of try
  catch(returnR rr){
    *errP = rr.errflag();
    return;
  }
}    // end of function bayesDensity

}    // end of extern "C"
