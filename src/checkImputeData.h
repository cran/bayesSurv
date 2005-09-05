#ifndef _CHECK_IMPUTE_DATA_H_
#define _CHECK_IMPUTE_DATA_H_

#include <R.h>

extern "C"{

const double _zero_ = 1e-10;

void
midimputeData(int* err,              double* t_impute,       const int* nP,
              const double* t_left,  const double* t_right,  const int* status);

void
midimputeDataDoubly(int* err,               double* t1_impute,       double* t2_impute,   const int* nP,
                    const double* t1_left,  const double* t1_right,  const int* status1,
                    const double* t2_left,  const double* t2_right,  const int* status2);

}

#endif
