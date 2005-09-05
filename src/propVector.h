#ifndef _PROP_VECTOR_H_
#define _PROP_VECTOR_H_

#include <R.h>
#include <Rmath.h>

#include "constants.h"

extern "C"{

double
logdUnif(const double* u,  const double* priorParmu);

void
rUnif(double* u,  const double* priorParmu);


void
transId(double* v, const double* u, const double* transParmu);

double
logJtransId(const double* u, const double* v, const double* transParmu);


void
transBeBeBe(double* v,  const double* u,  const double* transParmu);

void
invtransBeBeBe(double* u,  const double* v,  const double* transParmu);

double
logJtransBeBeBe(const double* u, const double* v, const double* transParmu);


void
transBrooks(double* v, const double* u, const double* transParmu);

void
invtransBrooks(double* u, const double* v, const double* transParmu);

double
logJtransBrooks(const double* u, const double* v, const double* transParmu);

void
transBeNG(double* v, const double* u, const double* transParmu);

void
invtransBeNG(double* u, const double* v, const double* transParmu);

double
logJtransBeNG(const double* u, const double* v, const double* transParmu);

}   // end of extern "C"

#endif


