// Routines to handle a proposal vector for the reversible moves

// 11/12/2003: start working on it
// 16/02/2004: completely changed

#include "bayessurvreg.h"

extern "C"{

using namespace std;


// logd*** ........... logarithm of a density of a distribution used to generate
//                     a canonical seed vector u
//                     * it may depend on additional parameters 'priorParmu'

// r*** .............. function to generate a random seed vector u
//                     * it may depend on additional parameters 'priorParmu'

// trans*** .......... function to transform a canonical seed vector u to a vector v
//                     which is then used to propose steps to a higher dimensional space,
//                     i.e. trans = v(u)
//                     * this transformation may depend on some additional parameters
//                       given by 'transParmu'

// invtrans*** ....... inverse of the transformation to get back a canonical seed, 
//                     i.e. invtrans = u(v)

// logJtrans*** ...... functions to give a log(Jacobian of a transformation), i.e.
//                     logJtrans = log(dv(u)/du)
//                     * it is a function of u, however sometimes it is simpler to express
//                       it as a function of v 
//                       -> that is why I give both u and v as parameters to logJtrans
//                          and use these that are more appropriate

// DIMENSIONS:
//   u ........... always a vector of length 3
//   v ........... always a vector of length 3
//   priorParmu .. dimension varies accross different proposals
//   transParmu .. dimension varies accross different proposals

// =============================================================================
//
// Basic density and random generator for a canonical seed vector u
//
// =============================================================================
// 

// ********** logdUnif *********
//
// Log-density of three independent standard uniform variates
// and a corresponding random generator
//
double
logdUnif(const double* u,  const double* priorParmu)
{
  return 0.0;
}

// *********** rUnif ************
//
// Sample independently three times from a standard uniform distribution
//
void
rUnif(double* u,  const double* priorParmu)
{
  for (int i = 0; i < 3; i++){
    u[i] = runif(0, 1);
  }
  return;
}


// =============================================================================
//
// Basic transformation (identity)
//
// =============================================================================

// ********** transId **********
//
void
transId(double* v, const double* u, const double* transParmu)
{
  for (int i = 0; i < 3; i++){
    v[i] = u[i];
  }
  return;
}


// ********** logJtransId *********
//
double
logJtransId(const double* u, const double* v, const double* transParmu)
{
  return 0.0;
}


// =============================================================================
//
// Functions to implement an original proposal for split-combine move 
// of Richardson and Green (1997)
// 
// =============================================================================
//
// Transformation of a canonical seed such that the resulting v vector
// has marginally Beta components, provided that the original canonical seed vector u
// has marginally standard uniform components
//
// transParmu ..... 6-dimensional vector giving Beta distributions of the transformed vector v,
//                  v[i] ~ Beta(transParmu[2i], transParmu[2i+1]), i = 0, 1, 2

// ********** transBeBeBe **********
//
void
transBeBeBe(double* v,  const double* u,  const double* transParmu)
{
  for (int i = 0; i < 3; i++){
    v[i] = qbeta(u[i], transParmu[2*i], transParmu[2*i + 1], 1, 0);
    if (v[i] <= NORM_ZERO) v[i] = NORM_ZERO;
    else                   if (v[i] >=  1 - NORM_ZERO) v[i] = 1 - NORM_ZERO;    
  }
  return;  
}

// ********* invtransBeBeBe *********
//
void
invtransBeBeBe(double* u,  const double* v,  const double* transParmu)
{
  for (int i = 0; i < 3; i++){
    u[i] = pbeta(v[i], transParmu[2*i], transParmu[2*i + 1], 1, 0);
  }
  return;
}

// ********* logJtransBeBeBe *********
//
double
logJtransBeBeBe(const double* u, const double* v, const double* transParmu)
{
  double logJ = 0.0;
  for (int i = 0; i < 3; i++){
    logJ -= dbeta(v[i], transParmu[2*i], transParmu[2*i + 1], 1);
  }
  return logJ;
}


// ==================================================================================
//
// Functions to implement a proposal for split-combine move
// of Brooks, Giudici and Roberts (2003), pp. 33
//
// ==================================================================================
//
// * (v[0], v[1], v[2]) ~ Beta(a[0], a[1]) * Beta(a[2], a[3]) * Beta(a[4], a[5]), where
//   a[i] = transParmu[i], i = 0, ..., 5
// * v[1] further transformed, v[1]' = |2*v[1] - 1|, i.e.
//   (v[0], v[1]', v[2]) = (v[0], |2*v[1] - 1|, v[2])
//   

// ********** transBrooks **********
//
void
transBrooks(double* v, const double* u, const double* transParmu)
{
  transBeBeBe(v, u, transParmu);
  v[1] = fabs(2*v[1] - 1);
  if (v[1] <= NORM_ZERO) v[1] = NORM_ZERO;
  else                   if (v[1] >=  1 - NORM_ZERO) v[1] = 1 - NORM_ZERO;    

  return;
}

// ********* invtransBrooks *********
void
invtransBrooks(double* u, const double* v, const double* transParmu)
{
  u[0] = pbeta(v[0], transParmu[0], transParmu[1], 1, 0);
  u[2] = pbeta(v[2], transParmu[4], transParmu[5], 1, 0);
  double unif = runif(0, 1);
  if (unif < 0.5) u[1] = -0.5*(v[1] - 1);
  else            u[1] = 0.5*(v[1] + 1);
  u[1] = pbeta(u[1], transParmu[2], transParmu[3], 1, 0);

  return;
}

// ******** logJtransBrooks *********
double
logJtransBrooks(const double* u, const double* v, const double* transParmu)
{

  double v1 = 0.5*(v[1] + 1);         // it does not matter whether I take 0.5*(v[1] + 1) or
                                      // -0.5*(v[1] - 1) due to symmetry of Beta(a, a) distribution
  double logJ = -dbeta(v1, transParmu[2], transParmu[3], 1) + LOG_2;
  logJ -= dbeta(v[0], transParmu[0], transParmu[1], 1);
  logJ -= dbeta(v[2], transParmu[4], transParmu[5], 1);
  
  return logJ;
}


// =============================================================================
//
// Functions to implement an original proposal for birth-death move 
// of Richardson and Green (1997)
// 
// =============================================================================
//
// Transformation of a canonical seed vector u such that the resulting vector v
// has marginally the following distributions
// v[0] ~ Beta(a, k)       (mixture weight)
// v[1] ~ N(xi, kappa)     (mixture mean)
// v[2] ~ Gamma(zeta, eta) (mixture inverse-variance), where
//
// (a, k, xi, sqrt(kappa), zeta, eta) = transParmu
//

// ******** transBeNG *********
//
void
transBeNG(double* v, const double* u, const double* transParmu)
{
  v[0] = qbeta(u[0], transParmu[0], transParmu[1], 1, 0);
  v[1] = qnorm(u[1], transParmu[2], transParmu[3], 1, 0);
  v[2] = qgamma(u[2], transParmu[4], 1/transParmu[5], 1, 0);

  if (v[0] <= NORM_ZERO) v[0] = NORM_ZERO;
  else                   if (v[0] >=  1 - NORM_ZERO) v[0] = 1 - NORM_ZERO;    

  if (v[1] <= -FLT_MAX) v[1] = -FLT_MAX;
  else                  if (v[1] >= FLT_MAX) v[1] = FLT_MAX;    

  if (v[2] <= NORM_ZERO) v[2] = NORM_ZERO;
  else                   if (v[2] >=  FLT_MAX) v[2] = FLT_MAX;    

  return;
}

// ********* invtransBeNG ********
//
void
invtransBeNG(double* u, const double* v, const double* transParmu)
{
  u[0] = pbeta(v[0], transParmu[0], transParmu[1], 1, 0);
  u[1] = pnorm(v[1], transParmu[2], transParmu[3], 1, 0);
  u[2] = pgamma(v[2], transParmu[4], 1/transParmu[5], 1, 0);
  return;  
}

// ******** logJtransBeNG *********
double
logJtransBeNG(const double* u, const double* v, const double* transParmu)
{
  double logJ = -dbeta(v[0], transParmu[0], transParmu[1], 1);
  logJ -= dnorm(v[1], transParmu[2], transParmu[3], 1);
  logJ -= dgamma(v[2], transParmu[4], 1/transParmu[5], 1);
  return logJ;
}
 

}   // end of extern "C"
