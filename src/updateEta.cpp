// Function to update the hyperparameter for the inverse variance (eta)
//   Gibbs step is used

// 28/11/2003: start woking on it: 28/11/2003

// This step is based on the following prior assumptions:
//    sigma2 ~ Inv-Gamma(zeta, eta)
//       eta ~ Gamma(g, h)

#include "bayessurvreg.h"

using namespace std;

void
updateEta(double* etaP,
          const int* kP,        const double* invsigma2M,
          const double* zetaP,  const double* gP,           const double* hP)
{
  double scale = *hP;
  for (int j = 0; j < *kP; j++) scale += invsigma2M[j];
  scale = 1/scale;
  if (scale <= SCALE_ZERO) scale = SCALE_ZERO;
  double shape = *gP + (*kP)*(*zetaP);
  *etaP = rgamma(shape, scale);
  return;
  
}   // end of function updateEta
