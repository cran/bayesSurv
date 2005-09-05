// Functions to check (doubly) censored data for consistency
//  and impute mid-points as initial event times
//
// 18/01/2005: 'midimputeData'
//             'midimputeDataDoubly'
//

#include "checkImputeData.h"
 
extern "C"{

/*** Simple mid-point imputation for censored data     ***/
/*** It performs some checks for consistency as well   ***/
void
midimputeData(int* err,              double* t_impute,       const int* nP,
              const double* t_left,  const double* t_right,  const int* status)
{
  *err = 0;
  int i;

  const int* stat = status;
  const double* t_L = t_left;
  const double* t_R = t_right;
  double* tImp = t_impute;

  for (i = 0; i < *nP; i++){
    switch (*stat){
    case 1:
    case 0:
      *tImp = *t_L;
      break;

    case 2:
      *tImp = 0.5*(*t_L);
      break;

    case 3:
      if (*t_R < *t_L || fabs(*t_R - *t_L) < _zero_){
        REprintf("\nError: time[%d] = (%f, %f],\n", i, *t_L, *t_L);
        *err = 3;
        return;
      }
      *tImp = 0.5*((*t_L) + (*t_R));
      break;

    default:
      REprintf("\nError: unknown censoring indicator\n");
      *err = 9;
      return;
    }
    stat++;
    t_L++;
    t_R++;
    tImp++;
  }

  return;
}  /** end of function imputeData **/


/*** Simple mid-point imputation for doubly-censored data     ***/
/*** It performs some checks for consistency as well          ***/
/***                                                          ***/
// t2_impute are imputed times-to-event (i.e. differences between event and onset times)
//
void
midimputeDataDoubly(int* err,               double* t1_impute,       double* t2_impute,   const int* nP,
                    const double* t1_left,  const double* t1_right,  const int* status1,
                    const double* t2_left,  const double* t2_right,  const int* status2)
{
  *err = 0;
  int i;

  const int* stat1 = status1;
  const int* stat2 = status2;
  const double* t1_L = t1_left;
  const double* t1_R = t1_right;
  const double* t2_L = t2_left;
  const double* t2_R = t2_right;
  double* t1Imp = t1_impute;
  double* t2Imp = t2_impute;

  for (i = 0; i < *nP; i++){
    switch (*stat1){

    /*** Onset exactly observed ***/
    case 1:
      *t1Imp = *t1_L;
      switch (*stat2){
      case 1:
        if (*t1_L > *t2_L || fabs(*t2_L - (*t1_L)) < _zero_){
          REprintf("\nError: onset[%d] = %f, event[%d] = %f\n", i, *t1_L, i, *t2_L);
          *err = 11;
          return;
        }
        *t2Imp = *t2_L - (*t1Imp);
        break;

      case 0:
        if (*t1_L > *t2_L){
          REprintf("\nError: onset[%d] = %f, event[%d] = %f+\n", i, *t1_L, i, *t2_L);
          *err = 10;
          return;
        }
        if (fabs(*t2_L - (*t1_L)) < _zero_) *t2Imp = 1000*_zero_;
        else                                *t2Imp = *t2_L - (*t1Imp);
        break;

      case 2:
        REprintf("\nError: onset[%d] = %f, event[%d] = %f-\n", i, *t1_L, i, *t2_L);
        *err = 12;
        return;

      case 3:
        if (*t1_L > *t2_L){
          REprintf("\nError: onset[%d] = %f, event[%d] = (%f, %f]\n", i, *t1_L, i, *t2_L, *t2_R);
          *err = 131;
          return;
        }
        if (*t2_R < *t2_L || fabs(*t2_R - *t2_L) < _zero_){
          REprintf("\nError: event[%d] = (%f, %f]\n", i, *t2_L, *t2_R);
          *err = 132;
          return;
        }
        *t2Imp = 0.5*(*t2_L + *t2_R) - *t1Imp;
        break;

      default:
        REprintf("\nError: unknown censoring indicator\n");
        *err = 9;
        return;
      }
      break;

    /*** Onset right censored ***/
    case 0:
      *t1Imp = *t1_L;
      switch (*stat2){
      case 1:
        REprintf("\nError: onset[%d] = %f+, event[%d] = %f\n", i, *t1_L, i, *t2_L);
        *err = 41;
        return;
  
      case 0:
        if (!(fabs(*t1_L - *t2_L) < _zero_)){
          REprintf("\nError: onset[%d] = %f+, event[%d] = %f+\n", i, *t1_L, i, *t2_L);
          *err = 40;
          return;
        }
        *t2Imp = 1000*_zero_;
        break;

      case 2:
        REprintf("\nError: onset[%d] = %f+, event[%d] = %f-\n", i, *t1_L, i, *t2_L);
        *err = 42;
        return;

      case 3:
        REprintf("\nError: onset[%d] = %f+, event[%d] = (%f, %f]\n", i, *t1_L, i, *t2_L, *t2_R);
        *err = 43;
        return;

      default:
        REprintf("\nError: unknown censoring indicator\n");
        *err = 9;
        return;
      }
      break;

    /*** Onset left censored ***/
    case 2:
      *t1Imp = 0.5*(*t1_L);
      switch (*stat2){
      case 1:
        if (*t2_L < *t1_L){
          REprintf("\nError: onset[%d] = %f-, event[%d] = %f\n", i, *t1_L, i, *t2_L);
          *err = 21;
          return;
        }
        *t2Imp = *t2_L - (*t1Imp);
        break;

      case 0:
        if (*t2_L < *t1_L){
          REprintf("\nError: onset[%d] = %f-, event[%d] = %f+\n", i, *t1_L, i, *t2_L);
          *err = 0;
          return;
        }
        *t2Imp = *t2_L - (*t1Imp);
        break;

      case 2:
        if (fabs(*t1_L - *t2_L) > _zero_){
          REprintf("\nError: onset[%d] = %f-, event[%d] = %f-\n", i, *t1_L, i, *t2_L);
          *err = 22;
          return;
        }
        *t2Imp = *t2_L - (*t1Imp);
        break;

      case 3:
        if (*t1_L > *t2_L){
          REprintf("\nError: onset[%d] = %f-, event[%d] = (%f, %f]\n", i, *t1_L, i, *t2_L, *t2_R);
          *err = 231;
          return;
        }
        if (*t2_R < *t2_L || fabs(*t2_R - *t2_L) < _zero_){
          REprintf("\nError: event[%d] = (%f, %f]\n", i, *t2_L, *t2_R);
          *err = 232;
          return;
        }
        *t2Imp = 0.5*(*t2_L + *t2_R) - *t1Imp;
        break;

      default:
        REprintf("\nError: unknown censoring indicator\n");
        *err = 9;
        return;
      }
      break;

    /*** Onset interval censored ***/
    case 3:
      if (*t1_R < *t1_L || fabs(*t1_R - *t1_L) < _zero_){
        REprintf("\nError: onset[%d] = (%f, %f]\n", i, *t1_L, *t1_R);
        *err = 3;
        return;
      }
      *t1Imp = 0.5*(*t1_L + *t1_R);
      switch (*stat2){
      case 1:
        if (*t2_L < *t1_R){
          REprintf("\nError: onset[%d] = (%f, %f], event[%d] = %f\n", i, *t1_L, *t1_R, *t2_L);
          *err = 31;
          return;
        }
        *t2Imp = *t2_L - *t1Imp;
        break;

      case 0:
        if (*t2_L < *t1_R){
          REprintf("\nError: onset[%d] = (%f, %f], event[%d] = %f+\n", i, *t1_L, *t1_R, *t2_L);
          *err = 30;
          return;
        }
        *t2Imp = *t2_L - *t1Imp;
        break;

      case 2:
        REprintf("\nError: onset[%d] = (%f, %f], event[%d] = %f-\n", i, *t1_L, *t1_R, i, *t2_L);
        *err = 32;
        return;

      case 3:
        if (fabs(*t2_L - *t1_L) < _zero_ && fabs(*t2_R - *t1_R) < _zero_){
          *t2Imp = *t2_R - *t1Imp;
        }
        else{
          if (*t2_R < *t2_L || fabs(*t2_R - *t2_L) < _zero_){
            REprintf("\nError: event[%d] = (%f, %f]\n", i, *t2_L, *t2_R);
            *err = 332;
            return;
          }
          if (*t2_L < *t1_R){
            REprintf("\nError: onset[%d] = (%f, %f], event[%d] = (%f, %f]\n", i, *t1_L, *t1_R, i, *t2_L, *t2_R);
            *err = 331;
            return;
          }
          *t2Imp = 0.5*(*t2_L + *t2_R) - *t1Imp;                    
        }
        break;

      default:
        REprintf("\nError: unknown censoring indicator\n");
        *err = 9;
        return;
      }
      break;
    }

    stat1++;
    stat2++;
    t1_L++;
    t1_R++;
    t2_L++;
    t2_R++;
    t1Imp++;
    t2Imp++;
  }

  return;
}    /** end of function checkData **/

}
