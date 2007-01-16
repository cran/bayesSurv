/*** classRandomEff_rhoNorm.h ***/

#ifndef _CLASS_RANDOM_EFF_RHO_NORM_H_
#define _CLASS_RANDOM_EFF_RHO_NORM_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "AK_BasicFun.h"

#include "classRandomEff.h"
#include "Mvtdist3.h"
#include "rhoNorm.h"
#include "Gspline.h"

void
Gspl_rho_intcpt_update(RandomEff *d,            RandomEff *b,              double *rho_zb,
                       double *regResOnset,     double *regResTime,        int *rho_accept,
                       const int *nP,           const int *rho_algor,      const double *rho_scaleL,
                       const Gspline *gg_d,     double** const mu_d,       const int *rM_d,
                       const Gspline *gg_b,     double** const mu_b,       const int *rM_b,
                       const Gspline *gg_zeta,  double** const mu_zeta,    const int *rM_zeta,
                       const Gspline *gg_eps,   double** const mu_eps,     const int *rM_eps);


#endif

