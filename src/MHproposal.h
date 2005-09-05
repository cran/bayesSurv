#ifndef _MH_PROPOSAL_H_
#define _MH_PROPOSAL_H_

#include <R.h>
#include <Rmath.h>

#include "AK_Error.h"
#include "mvtdist.h"

void
AMproposal(double* proppar,              double* chcovpar,     
           const double* covpar,         const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif,   
           const double* eps,            const double* sdNum,        const double* tolerChol);

void
AMadapt(double* covpar,     double* meanpar,     const double* par,
        const int* indUpd,  const int* nupdate,  const int* diagI,   const int* iter,
        const double* eps,  const double* sdNum);

void
MHproposal(double* proppar,              
           const double* chcovpar,       const double* par,          const int* indUpd,
           const int* npar,              const int* nupdate,         const int* diagI,             
           const double* halfRangeUnif,  const double* weightUnif);

#endif
