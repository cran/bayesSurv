// Functions to generate a sequence of possibly dependent vectors u
//  with possibly dependent components
//  where a stationary distribution of the vector u
//  is a (multivariate) standard uniform distribution.

// Using the method of a moody ring.
// see: Brooks, S.P., Giudici, P., and Roberts, G.O. (2003).
//      Efficient construction of reversible jump Markov chain
//      Monte Carlo proposal distributions (with Discussion).
//      JRSS B, vol. 65, pp. 3-55.
//      Section 6.4.

// I will write these functions such that they can be used also
// directly in R, i.e. with all parameters being pointers.

// 04/11/2003: start working on it

#include "bayessurvreg.h"

extern "C"{

using namespace std;

// ********** moodyRing **********
//
// PARAMETERS:
//
// moodP ........... mood parameter from the previous iteration: C(t-1)
//                   (ignored if corrP = false)
// timeDepP ........ epsilon \in [-0.5, 0.5] which determines dependence between consecutive vectors u
//                   epsilon = 0 <=> perfect dependence
// componentDepP ... delta \in [-0.5, 0.5] which determines dependence between the components of the vector u
//                   (ignored if corrP = false)
//                   delta = 0    <=> perfect dependence
//                   delta = 0.5  <=> corrP = false
// nuP ............. length of the u vector
// uA .............. u vector stored in array
//                   (if corrP = false, this gives vector u from the previous iteration)
//                   (ignored on input if corrP = true) 
// corrP ........... logical true = moody ring for correlated AV's is required
//                           false = moody ring for independent AV's is required
// callFromR ....... logical true = if this function is called directly from R (random seed must be got and put)
//
// RETURN:
//
// moodP ........... updated mood (if corrP = true)
//                   or -1 if there were some problems
// uA .............. vector u for the next iteration
//
void
moodyRing(double* uA,              double* moodP,             
          const double* timeDepP,  const double* componentDepP,     
          const int* nuP,          const int* corrP,
          const int* callFromR)
{
  try{
    if (*callFromR) GetRNGstate();

    int i;  
    double czt;

    double epsilon = fabs(*timeDepP);
    double delta = fabs(*componentDepP);

    // Input checks
    if (epsilon > 0.5){
      returnR error("C++ Error: epsilon for moody ring higher than 0.5", 99);
      throw error;
    }

    if (delta > 0.5){
      returnR error("C++ Error: delta for moody ring higher than 0.5", 99);
      throw error;
    }

    if (*moodP < 0 || *moodP > 1){
      *moodP = *moodP - floor(*moodP);
    }

    // Moody ring for correlated AV's
    if (*corrP){     
      double cwt = *moodP + runif(-epsilon, epsilon);
      *moodP = cwt - floor(cwt);      // new mood parameter (lying between 0 and 1)
      for (i = 0; i < *nuP; i++){
        czt = *moodP + runif(-delta, delta);
        uA[i] = czt - floor(czt);
      }
    }

    // Moody ring for independent AV's
    else{
      for (i = 0; i < *nuP; i++){
        czt = uA[i] + runif(-epsilon, epsilon);
        uA[i] = czt - floor(czt);
      }   
    }

    if (*callFromR) PutRNGstate();
    return;

  }  // end of try
  catch(returnR){
    *moodP = -1;
    if (*callFromR){
      PutRNGstate();
      return;
    }
    throw;
  }  

}   // end of function moodyRing


// ********** corr_moodyRing **********
// * generate a sequence of n dependent vectors whose components are also dependent
//   using a moody ring
//
// PARAMETERS:
//
// initmoodP ........... initial mood parameter C(0) (lying between 0 and 1)
// timeDepP ............ epsilon \in [-0.5, 0.5] which determines dependence between consecutive vectors u
// componentDepP ....... delta \in [-0.5, 0.5] which determines dependence between the components of the vector u
// nuP ................. length of the u vector
// nP .................. number of u vectors that should be generated
// callFromR ........... logical true = if this function is called directly from R (random seed must be got and put)
//
// RETURN:
//
// uA .................. generated u vectors stored in array,
//                       a pointer to an array of length (*nuP) * (*nP) must be supplied
// moodA ............... mood parameters used to generate the u vectors,
//                       a pointer to an array of length *nP must be supplied
// initmoodP ........... last used mood parameter or -1 if there were some problems
//
void
corr_moodyRing(double* uA,              double* moodA,                double* initmoodP,
               const double* timeDepP,  const double* componentDepP,  const int* nuP,
               const int* nP,           const int* callFromR)
{
  if (*callFromR) GetRNGstate();

  for (int i = 0; i < *nP; i++){
    moodyRing(uA + i*(*nuP), initmoodP, timeDepP, componentDepP, nuP, &ONE_INT, &ZERO_INT);
    if (*initmoodP < 0){
      if (*callFromR) PutRNGstate();
      return;
    }
    moodA[i] = *initmoodP;    
  }

  if (*callFromR) PutRNGstate();
  return;

}   // end of function corr_moodyRing


// ********** indep_moodyRing **********
// * generate a sequence of n dependent vectors whose components are independent
//   using a moody ring
//
// PARAMETERS:
//
// inituA .............. initial u vector (array of length *nuP)
// timeDepP ............ epsilon \in [-0.5, 0.5] which determines dependence between consecutive vectors u
// nuP ................. length of the u vector
// nP .................. number of u vectors that should be generated
// callFromR ........... logical true = if this function is called directly from R (random seed must be got and put)
//
// RETURN:
//
// uA .................. generated u vectors stored in array,
//                       a pointer to an array of length (*nuP) * (*nP) must be supplied
// inituA .............. -1 in its nulth component if there were some problems
//
//
void
indep_moodyRing(double* uA,               double* inituA,
                const double* timeDepP,   const int* nuP,
                const int* nP,            const int* callFromR)
{
  try{
    if (*callFromR) GetRNGstate();

    double epsilon = fabs(*timeDepP);

    // Input checks
    if (epsilon > 0.5){
      returnR error("C++ Error: epsilon for moody ring higher than 0.5", 99);
      throw error;
    }

    int i, j;
    double czt;
    
    // First iteration   
    for (j = 0; j < *nuP; j++){
      czt = inituA[j] + runif(-epsilon, epsilon);
      uA[j] = czt - floor(czt);
    }

    // All other iterations
    for (i = 1; i < *nP; i++){
      for (j = 0; j < *nuP; j++){
        czt = uA[(i-1)*(*nuP) + j] + runif(-epsilon, epsilon);
        uA[i*(*nuP) + j] = czt - floor(czt);
      }
    }   

    if (*callFromR) PutRNGstate();
    return;

  }   // end of try
  catch(returnR){
    inituA[0] = -1;
    if (*callFromR){
      PutRNGstate();
      return;
    }
    throw;
  }  

}   // end of function indep_moodyRing 
              

}   // end of extern "C"
