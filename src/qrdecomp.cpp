// QR decomposition of a matrix
//
// Rewritten from Fortran routines
// 
// 03/02/2004: start working on it 
//
#include "qrdecomp.h"

extern "C"{

using namespace std;

// Main routine originally from /usr/lib/R/src/appl/dqrdc2.f:
//
//
//    dqrdc2 uses householder transformations to compute the qr
//    factorization of an n by p matrix x.  a limited column
//     pivoting strategy based on the 2-norms of the reduced columns
//     moves columns with near-zero norm to the right-hand edge of
//     the x matrix.  this strategy means that sequential one
//     degree-of-freedom effects can be computed in a natural way.
//
//     i am very nervous about modifying linpack code in this way.
//     if you are a computational linear algebra guru and you really
//     understand how to solve this problem please feel free to
//     suggest improvements to this code.
//
//     Another change was to compute the rank.
//
//     on entry
//
//        x       double precision(ldx,p), where ldx .ge. n.
//                x contains the matrix whose decomposition is to be
//                computed.
//
//       ldx     integer.
//               ldx is the leading dimension of the array x.  --> this was removed from the arguments in C++ version
//
//        n       integer.
//                n is the number of rows of the matrix x.
//
//        p       integer.
//                p is the number of columns of the matrix x.
//
//        tol     double precision
//                tol is the nonnegative tolerance used to
//                determine the subset of the columns of x
//                included in the solution.
//
//        jpvt    integer(p).
//                integers which are swapped in the same way as the
//                the columns of x during pivoting.  on entry these
//                should be set equal to the column indices of the
//                columns of the x matrix (typically 1 to p).
//
//        work    double precision(p,2).
//                work is a work array.    --> this was removed from the arguments in C++ version
//
//    on return
//
//        x       x contains in its upper triangle the upper
//                triangular matrix r of the qr factorization.
//                below its diagonal x contains information from
//                which the orthogonal part of the decomposition
//                can be recovered.  note that if pivoting has
//                been requested, the decomposition is not that
//                of the original matrix x but that of x
//                with its columns permuted as described by jpvt.
//
//        k       integer.
//                k contains the number of columns of x judged
//                to be linearly independent.
//
//        qraux   double precision(p).
//                qraux contains further information required to recover
//                the orthogonal part of the decomposition.
//
//        jpvt    jpvt(k) contains the index of the column of the
//               original matrix that has been interchanged into
//                the k-th column.
//
//     this version dated 22 august 1995
//     ross ihaka
//
//     bug fixes 29 September 1999 BDR (p > n case, inaccurate ranks)
//
//
//     dqrdc uses the following functions and subprograms.
//
//     blas daxpy,ddot,dscal,dnrm2
//     fortran dabs,dmax1,min0,dsqrt
//
//   rewritten to C++ 3 February 2004
//   Arnost Komarek
//
// =================================================================================
void 
dqrdc2CPP(double* x, const int* n, const int* p, const double* tol, int* k, double* qraux, int* jpvt)
{
  double* work = (double*) calloc((*p)*2, sizeof(double));
  if (!work) throw returnR("Could not allocate working space for dqrdc2CPP", 1);
  int i, j, l, lp1, lup;
  double tt,ttt;
  double nrmxl,t;

//
//     compute the norms of the columns of x.
//
  for (j = 0; j < *p; j++){
    qraux[j] = dnrm2CPP(*n, x + j*(*n), 1);
    work[j] = qraux[j];
    work[(*p)+j] = qraux[j];
    if(work[(*p)+j] == 0.0) work[(*p)+j] = 1.0;
  }

//
//     perform the householder reduction of x.
//
  lup = (*n < *p ? *n : *p);
  *k = *p + 1;
  for (l = 0; l < lup; l++){     // DO 200
//
//     previous version only cycled l to lup
//
//     cycle the columns from l to p left-to-right until one
//     with non-negligible norm is located.  a column is considered
//     to have become negligible if its norm has fallen below
//     tol times its original norm.  the check for l .le. k
//     avoids infinite cycling.
//
    while (l < *k - 1 && qraux[l] < work[(*p) + l]*(*tol)){
      lp1 = l+1;
      for (i = 0; i < *n; i++){
        t = x[l*(*n) + i];
        for (j = lp1; j < (*p); j++){
          x[(j-1)*(*n) + i] = x[j*(*n) + i];
        }
        x[(*p-1)*(*n)+i] = t;
      }
      i = jpvt[l];
      t = qraux[l];
      tt = work[l];
      ttt = work[(*p)+l];
      for (j = lp1; j < *p; j++){
        jpvt[j-1] = jpvt[j];
        qraux[j-1] = qraux[j];
        work[j-1] = work[j];
        work[(*p)+j-1] = work[(*p)+j];
      }    
      jpvt[*p-1] = i;
      qraux[*p-1] = t;
      work[*p-1] = tt;
      work[(*p)+(*p)-1] = ttt;
      *k = *k - 1;
    }
    if (l == *n - 1) continue;          // continue DO 200

//
//          compute the householder transformation for column l.
//
    nrmxl = dnrm2CPP(*n-l, x + (*n)*l + l, 1);

    if (nrmxl == 0.0) continue;    // continue DO 200
    if (x[(*n)*l + l] != 0.0) nrmxl = (x[(*n)*l + l] > 0.0 ? nrmxl : -nrmxl);
    dscalCPP(*n-l, 1.0/nrmxl, x + (*n)*l + l, 1);
    x[(*n)*l + l] = 1.0 + x[(*n)*l + l];
//
//              apply the transformation to the remaining columns,
//              updating the norms.
//
    lp1 = l + 1;

    if (*p - 1 < lp1){
//
//              save the transformation.
//
      qraux[l] = x[(*n)*l + l];
      x[(*n)*l + l] = -nrmxl; 
      continue;                  // continue DO 200
    }

    for (j = lp1; j < *p; j++){    // DO 160
      t = -ddotCPP(*n-l, x+(*n)*l+l, 1, x+(*n)*j+l, 1)/ x[(*n)*l+l];
      daxpyCPP(*n-l, t, x+(*n)*l+l, 1, x+(*n)*j+l, 1);
      if (qraux[j] == 0.0) continue;      // continue DO 160
      tt = 1.0 - (fabs(x[(*n)*j+l])/qraux[j])*(fabs(x[(*n)*j+l])/qraux[j]);
      tt = (tt > 0.0 ? tt: 0.0);
      t = tt;
//
// modified 9/99 by BDR. Re-compute norms if there is large reduction
// The tolerance here is on the squared norm
// In this version we need accurate norms, so re-compute often.
//  work(j,1) is only updated in one case: looks like a bug -- no longer used
//
//                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j,1))**2
//                     if (tt .eq. 1.0d0) go to 130
  
      if (fabs(t) < 1e-6){
        qraux[j] = dnrm2CPP(*n-l-1, x+(*n)*j+l+1, 1); 
        work[j] = qraux[j];
      }
      else{
        qraux[j] = qraux[j]*sqrt(t);
      }
    }    // end of DO 160

//
//              save the transformation.
//
    qraux[l] = x[(*n)*l + l];
    x[(*n)*l + l] = -nrmxl;
  }  // end of DO 200

  *k = (*k-1 < *n ? *k-1 : *n);

  free(work);
  return;
}    // end of the function dqrdc2




// Helping routines originally from /usr/lib/R/src/appl/blas.f:
//
// =============================================================
//
//     forms the dot product of two vectors.
//     uses unrolled loops for increments equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03, Arnost Komarek
//     n ............ how many components are to be multiplied and summed
//     dx, dy ....... vectors
//     incx, incy ... increments

// arrays: dx, dy
// =============================================================
double
ddotCPP(const int n, double* dx, const int incx, double* dy, const int incy)
{
   double dtemp;
   int i, ix, iy, m, mp1;

   double ddot = 0.0;
   dtemp = 0.0;
   if(n <= 0) return ddot;
   if(!(incx == 1 && incy == 1)){
//
//        code for unequal increments or equal increments
//          not equal to 1
//
      ix = 1;
      iy = 1;
      if(incx < 0) ix = (-n+1)*incx + 1;
      if(incy < 0) iy = (-n+1)*incy + 1;
      for (i = 1; i <= n; i++){
        dtemp += dx[ix-1]*dy[iy-1];
        ix += incx;
        iy += incy;
      }
      ddot = dtemp;
      return ddot;
   }
   else{
//
//        code for both increments equal to 1
//
//
//        clean-up loop
//
      m = n % 5;
      if(m != 0)
         for (i = 1; i <= m; i++)
           dtemp += dx[i-1]*dy[i-1];
      if(n >= 5){
         mp1 = m + 1;
         for (i = mp1; i <= n; i+=5)
           dtemp += dx[i - 1]*dy[i - 1] + dx[i]*dy[i] +
                dx[i + 1]*dy[i + 1] + dx[i + 2]*dy[i + 2] + dx[i + 3]*dy[i + 3];
      }
      ddot = dtemp;
      return ddot;
   }
}

// =============================================================
//
//     constant times a vector plus a vector.
//     uses unrolled loops for increments equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03, Arnost Komarek
//     Input
//          n ....... how many components are to be multiplied and summed
//          da ...... scalar
//          dx ...... vector
//          incx .... increment for the vector dx
//          dy ...... second vector
//          incy .... increment for the second vector
//
//     Output
//          dy
// =============================================================
void
daxpyCPP(const int n, const double da, double* dx, const int incx, double* dy, const int incy)
{
      int i, ix, iy, m, mp1;

      if (n <= 0) return;
      if (da == 0.0) return;
      if(!(incx == 1 && incy ==1)){
//
//        code for unequal increments or equal increments
//          not equal to 1
//
         ix = 1;
         iy = 1;
         if(incx < 0) ix = (-n+1)*incx + 1;
         if(incy < 0) iy = (-n+1)*incy + 1;
         for (i = 1; i <= n; i++){
           dy[iy - 1] += da*dx[ix - 1];
           ix += incx;
           iy += incy;
         }
         return;
      }
      else{
//
//        code for both increments equal to 1
//
//
//        clean-up loop
//
         m = n % 4;
         if(m != 0)
            for (i = 1; i <= m; i++)
               dy[i - 1] += da*dx[i - 1];
         if(n < 4) return;
         mp1 = m + 1;
         for (i = mp1; i<= n; i+=4){
           dy[i - 1] += da*dx[i - 1];
           dy[i] += da*dx[i];
           dy[i + 1] += da*dx[i + 1];
           dy[i + 2] += da*dx[i + 2];
         }
         return;
      }
}


// =============================================================
//
//     scales a vector by a constant.
//     uses unrolled loops for increment equal to one.
//     jack dongarra, linpack, 3/11/78.
//     modified 3/93 to return if incx .le. 0.
//     modified 12/3/93, array(1) declarations changed to array(*)
//     rewritten to C++ 5/5/03, Arnost Komarek
//
// =============================================================
void
dscalCPP(const int n, const double da, double* dx, const int incx)
{
      int i, m, mp1, nincx;

      if(n <= 0 || incx <= 0) return;
      if(incx != 1){
//
//        code for increment not equal to 1
//
         nincx = n*incx;
         for (i = 1; i <= nincx; i+=incx)
            dx[i - 1] *= da;
         return;
      }
      else{
//
//        code for increment equal to 1
//
//
//        clean-up loop
//

         m = n % 5;
         if(m != 0 )
            for (i = 1; i <= m; i++)
              dx[i - 1] *= da;
         if(n < 5) return;
         mp1 = m + 1;
         for (i = mp1; i<= n; i+=5){
           dx[i - 1] *= da;
           dx[i]  *= da;
           dx[i + 1] *= da;
           dx[i + 2] *= da;
           dx[i + 3] *= da;
         }
         return;
      }
}


// ===================================================================
//
//  DNRM2 returns the euclidean norm of a vector via the function
//  name, so that
//
//     DNRM2 := sqrt( x'*x )
//
//  -- This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to DLASSQ.
//     Sven Hammarling, Nag Ltd.
// 
//     rewritten to C++ 2/2/04
//     Arnost Komarek
// 
// ==================================================================
double 
dnrm2CPP(const int n, const double* x, const int incx)
{
  const double ONE = 1.0;
  const double ZERO = 0.0;
  int IX;
  double ABSXI, NORM, SCALE, SSQ;

  if (n < 1 || incx < 1)
    NORM = ZERO;
  else 
    if (n == 1)
      NORM = fabs(x[0]);
    else{
      SCALE = ZERO;
      SSQ = ONE;
            // The following loop is equivalent to this call to the LAPACK
            //  auxiliary routine:
            //  CALL DLASSQ( N, x, incx, SCALE, SSQ )
      for (IX = 1; IX <= 1 + ( n - 1 )*incx; IX += incx){
        if (x[IX-1] != ZERO){
          ABSXI = fabs(x[IX-1]);
          if(SCALE < ABSXI){
	    SSQ   = ONE   + SSQ*( SCALE/ABSXI )*( SCALE/ABSXI );
	    SCALE = ABSXI;
	  }
          else
            SSQ   = SSQ   +     ( ABSXI/SCALE )*( ABSXI/SCALE );
        }
      }
      NORM = SCALE * sqrt(SSQ);    
    }

  return NORM;
}

}   // end of extern "C"
