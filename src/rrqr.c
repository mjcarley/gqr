/* Copyright (C) 2020 by Michael Carley */

/**********************************************************************
 *
 * This file is part of gqr.
 *
 * gqr is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gqr is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gqr.  If not, see <http://www.gnu.org/licenses/>.
 *
 **********************************************************************/

/*
  rank revealing QR factorizations
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <blaswrap.h>

#include "gqr.h"
#include "gqr-private.h"

gint dgeqp3_(gint *m, gint *n, gdouble *A, gint *lda, gint *jpvt,
	     gdouble *tau, gdouble *work, gint *lwork, gint *info) ;

gint dlarf_(gchar *side, gint *m, gint *n, gdouble *v, gint *incv,
	    gdouble *tau, gdouble *C, gint *ldc, gdouble *work) ;

gint dtrtrs_(gchar *uplo, gchar *trans, gchar *diag,
	     gint *n, gint *nrhs, gdouble *a, gint *lda, gdouble *b,
	     gint *ldb, gint *info) ;
gint dtrtri_(gchar *uplo, gchar *diag, gint *n, gdouble *a,
	     gint *lda, gint *info) ;
void drotg_(gdouble *da, gdouble *db, gdouble *c, gdouble *s) ;
void drot_(gint *N, gdouble *dx, gint *incx, gdouble *dy, gint *incy,
	   gdouble *C, gdouble *S) ;

#define matrix_index(_m,_n,_i,_j) ((_i) + (_j)*(_m))

gint column_swap(gdouble *a, gint m, gint n, gint j, gint k)

{
  gint one = 1 ;
  
  blaswrap_dswap(m,
		 &(a[matrix_index(m,n,0,j)]), one,
		 &(a[matrix_index(m,n,0,k)]), one) ;
  return 0 ;
}

gint column_shuffle(gdouble *a, gint m, gint n, gint j0, gint lda)

{
  gint i, j ;
  gdouble tmp ;

  for ( i = 0 ; i < m ; i ++ ) {
    tmp = a[matrix_index(lda,n,i,j0)] ;
    for ( j = j0 ; j < n ; j ++ ) {
      a[matrix_index(lda,n,i,j)] = a[matrix_index(lda,n,i,j+1)] ;
    }
    a[matrix_index(lda,n,i,n)] = tmp ;
  }

  return 0 ;
}  

gint row_shuffle(gdouble *a, gint m, gint n, gint i0, gint i1, gint lda)

{
  gint i, j ;
  gdouble tmp ;

  for ( j = 0 ; j < n ; j ++ ) {
    tmp = a[matrix_index(m,n,i0,j)] ;
    for ( i = i0 ; i < i1 ; i ++ ) {
      a[matrix_index(m,n,i,j)] = a[matrix_index(lda,n,i+1,j)] ;
    }
    a[matrix_index(m,n,i1,j)] = tmp ;
  }

  return 0 ;
}  

gint array_shuffle_int(gint *a, gint n, gint i0)

{
  gint tmp, i ;
  
  tmp = a[i0] ;
  
  for ( i = i0 ; i < n ; i ++ ) a[i] = a[i+1] ;
  a[n] = tmp ;
    
  return 0 ;
}

gint array_shuffle_double(gdouble *a, gint n, gint i0)

{
  gint i ;
  gdouble tmp ;
  
  tmp = a[i0] ;
  
  for ( i = i0 ; i < n ; i ++ ) a[i] = a[i+1] ;
  a[n] = tmp ;
    
  return 0 ;
}

gint gqr_srrqr(gdouble *A, gint m, gint n, gdouble f, gdouble tol,
	       gdouble *Q, gdouble *R, gint *pvt, gint *rank, gint *ldr,
	       gdouble *work, gint lwork)

/*
 * algorithm of Gu and Eisenstat, based on the Matlab code of Xu Xing
 *
 */
  
{
  gdouble tau[8192], *AB, *gamma, omega[1024], tmp ;
  gdouble ta, tb, Gs, Gc ;
  gint info, i, j, k, mab, nab, n12, ii, itmp, iter, one = 1 ;
  gint nr, mr ;
  gint cp ;
  cp = 0 ;
  
  /* lwork = 8192 ; */
  /*do the pivoting QR factorization*/
  memset(pvt, 0, MAX(m,n)*sizeof(gint)) ;
  memset(tau, 0, 8192*sizeof(gdouble)) ;
  dgeqp3_(&m, &n, A, &m, pvt, tau, work, &lwork, &info) ;

  /*subtract one from pvt to convert to C indexing*/
  for ( i = 0 ; i < n ; i ++ ) pvt[i] -= 1 ;
  
  /*find the rank*/
  for ( k = 1 ; k < MIN(m,n) ; k ++ ) {
    if ( fabs(A[matrix_index(m,n,k,k)]) < tol ) break ;
  }

  /*size R*/
  /* mr = n ; nr = n ; */
  mr = k+3 ; nr = n ;
  
  /* copy the useful bit of A into R */
  for ( j = 0 ; j < nr ; j ++ ) {
    for ( i = 0 ; i <= MIN(mr,j) ; i ++ ) {
      R[matrix_index(mr, nr, i, j)] = A[matrix_index(m,n,i,j)] ;
    }
  }

  /*fill Q with the matrix generated by the reflectors*/
  memset(Q, 0, k*m*sizeof(gdouble)) ;
  for ( j = 0 ; j < k ; j ++ ) {
    Q[matrix_index(m, k, j, j)] = 1.0 ;
  }
  for ( i = k-1 ; i >= 0 ; i -- ) {
    gint nv, nc ;
    A[matrix_index(m,n,i,i)] = 1.0 ;
    nv = m - i ; nc = k - i ;
    dlarf_("L", &nv, &nc, &(A[matrix_index(m,n,i,i)]), &one,
  	   &(tau[i]), &(Q[matrix_index(m,k,i,i)]), &m, work) ;
  }
  
  /*make diagonals of R positive*/
  for ( i = 0 ; i < MIN(mr,nr) ; i ++ ) {
    tmp = -1.0 ;
    if ( R[matrix_index(mr,nr,i,i)] < 0 ) {
      blaswrap_dscal(nr, tmp, &(R[matrix_index(mr,nr,i,0)]), mr) ;
      blaswrap_dscal(m, tmp, &(Q[matrix_index(m,n,0,i)]), one) ;
    }
  }

  /*number of columns in R12*/
  n12 = nr - k ;
  
  /*linear solve with R11 and R12 packed in A upon exit from dgeqp3*/
  nab = nr - k ; mab = k ;
  AB = (gdouble *)g_malloc0(nab*mab*sizeof(gdouble)) ;
  for ( i = 0 ; i < mab ; i ++ ) {
    for ( j = 0 ; j < nab ; j ++ ) {
      AB[matrix_index(mab,nab,i,j)] = R[matrix_index(mr,nr,i,k+j)] ;
    }
  }

  /*initialize AB = A^{-1}B*/
  dtrtrs_("U", "N", "N", &k, &nab, R, &mr, AB, &mab, &info) ;

  gamma = (gdouble *)g_malloc0(n*sizeof(gdouble)) ;
  if ( k != m ) {
    for ( j = 0 ; j < n-k ; j ++ ) {
      i = n-k ;
      gamma[j] = blaswrap_dnrm2(i,&(A[matrix_index(m,n,k,k+j)]), one) ;
      /* g_assert(!isnan(gamma[j])) ; */
    }
  } else {
    memset(gamma, 0, (n-k)*sizeof(gdouble)) ;
  }

  /*copy R11 into work for inversion*/
  for ( i = 0 ; i < k ; i ++ ) {
    for ( j = 0 ; j < k ; j ++ ) {
      work[matrix_index(k,k,i,j)] = R[matrix_index(mr,nr,i,j)] ;
    }
  }

  dtrtri_("U", "N", &k, work, &k, &info) ;

  /* fprintf(stderr, "%lg\n", work[matrix_index(k,k,k-1,k-1)]) ; exit(0) ; */
  
  for ( i = 0 ; i < k ; i ++ ) {
    omega[i] =
      1.0/blaswrap_dnrm2(k, &(work[matrix_index(k,k,i,0)]), k) ;
  }
  
  iter = 0 ;
  while ( 1 ) {
    while ( 1 ) {
      iter ++ ;
      gint jj, Rm ;

      /* Rm = MIN(mr,nr) ; */
      Rm = MIN(m,n) ;
      for ( j = 0 ; j < n12 ; j ++ ) {
	for ( i = 0 ; i < k ; i ++ ) {
	  tmp = gamma[j]/omega[i] ;
	  tmp *= tmp ;
	  tmp += AB[matrix_index(mab,nab,i,j)]*AB[matrix_index(mab,nab,i,j)] ;
	  if ( tmp > f*f ) { break ; }
	}
	if ( tmp > f*f ) { break ; }
      }

      /*no entry > f*f, so RRQR is strong */
      if ( tmp < f*f ) { break ; }

      if ( j > 0 ) {
	/*swap columns k and k+j*/
	tmp = gamma[0] ; gamma[0] = gamma[j] ; gamma[j] = tmp ;
	itmp = pvt[k] ; pvt[k] = pvt[k+j] ; pvt[k+j] = itmp ;
	column_swap(AB, mab, nab, 0, j) ;
	column_swap(R, mr, nr, k+0, k+j) ;
      }

      if ( i < k-1 ) {
	array_shuffle_int(pvt, k-1, i) ;
	array_shuffle_double(omega, k-1, i) ;
	column_shuffle(R, mr, k-1, i, mr) ;
	row_shuffle(AB, mab, nab, i, k-1, k) ;
	
	/*Givens rotations*/
	for ( ii = i ; ii < k-1 ; ii ++ ) {
	  ta = R[matrix_index(mr,nr,ii  ,ii)] ;
	  tb = R[matrix_index(mr,nr,ii+1,ii)] ;
	  drotg_(&ta, &tb, &Gc, &Gs) ;
	  if ( Gc*R[matrix_index(mr,nr,ii  ,ii)] +
	       Gs*R[matrix_index(mr,nr,ii+1,ii)] < 0.0 ) {
	    Gc = -Gc ; Gs = -Gs ;
	  }
	  drot_(&nr,
		&(R[matrix_index(mr,nr,ii  ,0)]), &mr,
		&(R[matrix_index(mr,nr,ii+1,0)]), &mr,
		&Gc, &Gs) ;
	  drot_(&m,
	  	&(Q[matrix_index(m,k,0,ii  )]), &one,
	  	&(Q[matrix_index(m,k,0,ii+1)]), &one,
	  	&Gc, &Gs) ;
	}
	if ( R[matrix_index(mr,nr,k-1,k-1)] < 0 ) {
	  for ( ii = 0 ; ii < nr ; ii ++ )
	    R[matrix_index(mr,nr,k-1,ii)] *= -1 ;
	  for ( ii = 0 ; ii < m ; ii ++ )
	    Q[matrix_index(m,k,ii,k-1)] *= -1 ;
	}
      }
      /* fprintf(stderr, "%lg\n", R[matrix_index(mr,nr,k-1,k-1)]) ; exit(0) ; */
      
      if ( k < Rm ) {
	/* fprintf(stderr, "k=%d, Rm=%d, mr=%d\n", k, Rm, mr) ; */
	for ( ii = k + 1 ; ii < Rm ; ii ++ ) {
	  ta = R[matrix_index(mr,nr,k ,k)] ;
	  tb = R[matrix_index(mr,nr,ii,k)] ;
	  drotg_(&ta, &tb, &Gc, &Gs) ;
	  
	  if ( Gc*R[matrix_index(mr,nr, k,k)] +
	       Gs*R[matrix_index(mr,nr,ii,k)] < 0.0 ) {
	    Gc = -Gc ; Gs = -Gs ;
	  }
	  drot_(&nr,
		&(R[matrix_index(mr,nr, k,0)]), &mr,
		&(R[matrix_index(mr,nr,ii,0)]), &mr,
		&Gc, &Gs) ;
	  drot_(&m,
	  	&(Q[matrix_index(m,k,0,k )]), &one,
	  	&(Q[matrix_index(m,k,0,ii)]), &one,
	  	&Gc, &Gs) ;
	}
      }
      /* fprintf(stderr, "%lg\n", R[matrix_index(mr,nr,k-1,k-1)]) ; exit(0) ; */

      itmp = pvt[k-1] ; pvt[k-1] = pvt[k] ; pvt[k] = itmp ;
      gdouble ga, mu, nu, rho, ga_bar ;
      gdouble *b1, *b2, *u1, *u2, *c1T, *c2T, *c1Tbar, *c2Tbar ;

      b1 = &(work[0]) ; b2 = &(b1[n]) ; u1 = &(b2[n]) ; u2 = &(u1[n]) ;
      c1T = &(u2[n]) ; c2T = &(c1T[n]) ;
      c1Tbar = &(c2T[n]) ; c2Tbar = &(c1Tbar[n]) ;
      
      ga = R[matrix_index(mr, nr, k-1, k-1)] ;
      mu = R[matrix_index(mr, nr, k-1, k  )]/ga ;
      if ( k < Rm ) {
	nu = R[matrix_index(mr, nr, k, k)]/ga ;
      } else {
	nu = 0.0 ;
      }
      rho = sqrt(mu*mu + nu*nu) ;
      /* fprintf(stderr, "%lg %lg\n", mu, nu) ; exit(0) ; */
      
      ga_bar = ga*rho ;

      /* fprintf(stderr, "%lg %lg\n", ga, rho) ; exit(0) ; */
      
      ii = k - 1 ;
      /* fprintf(stderr, "%d %d\n", iter, ii) ; */
      blaswrap_dcopy(ii, &(R[matrix_index(mr, nr, 0, k-1)]), one, b1, one) ;
      blaswrap_dcopy(ii, &(R[matrix_index(mr, nr, 0, k  )]), one, b2, one) ;
      blaswrap_dcopy(ii, b2, one, &(R[matrix_index(mr, nr, 0, k-1)]), one) ;
      blaswrap_dcopy(ii, b1, one, &(R[matrix_index(mr, nr, 0, k  )]), one) ;

      /* fprintf(stderr, "%lg\n", R[matrix_index(mr,nr,k-1,k-1)]) ; exit(0) ; */
      ii = nr - k - 1 ;
      blaswrap_dcopy(ii, &(R[matrix_index(mr, nr, k-1, k+1)]), mr, c1T, one) ;
      if ( k+1 > Rm ) {
	memset(c2T, 0, ii*sizeof(gdouble)) ;	
      } else {
	blaswrap_dcopy(ii, &(R[matrix_index(mr, nr, k, k+1)]), mr, c2T, one) ;
      }

      for ( ii = 0 ; ii < n - k - 1 ; ii ++ ) {
      	c1Tbar[ii] = (mu*c1T[ii] + nu*c2T[ii])/rho ;
      	c2Tbar[ii] = (nu*c1T[ii] - mu*c2T[ii])/rho ;
      }

      /* fprintf(stderr, "%lg\n", R[matrix_index(mr,nr,k-1,k-1)]) ; exit(0) ; */

      R[matrix_index(mr,nr,k-1,k-1)] = ga_bar ;
      R[matrix_index(mr,nr,k-1,k  )] = ga*mu/rho ;
      R[matrix_index(mr,nr,k  ,k  )] = ga*nu/rho ;

      ii = n - k -1 ;
      blaswrap_dcopy(ii, c1Tbar, one,
		     &(R[matrix_index(mr, nr, k-1, k+1)]), mr) ;
      blaswrap_dcopy(ii, c2Tbar, one,
		     &(R[matrix_index(mr, nr, k  , k+1)]), mr) ;

      itmp = 1 ; jj = k - 1 ;
      dtrtrs_("U", "N", "N", &jj, &one, R, &mr, b1, &n, &info) ;
      for ( ii = 0 ; ii < k - 1 ; ii ++ ) {
	u1[ii] = AB[matrix_index(mab,nab,ii,0)] ;
	AB[matrix_index(mab,nab,ii,0)] = (nu*nu*b1[ii] - mu*u1[ii])/rho/rho ;
      }
      AB[matrix_index(mab,nab,k-1,0)] = mu/rho/rho ;
      for ( ii = 1 ; ii < nab ; ii ++ ) {
	AB[matrix_index(mab,nab,k-1,ii)] = c1Tbar[ii-1]/ga_bar ;
      }
      for ( ii = 0 ; ii < k-1 ; ii ++ ) {
	for ( jj = 1 ; jj < nab ; jj ++ ) {
	  AB[matrix_index(mab,nab,ii,jj)] +=
	    (nu*b1[ii]*c2Tbar[jj-1] - u1[ii]*c1Tbar[jj-1])/ga_bar ;
	}
      }

      gamma[0] = ga*nu/rho ;
      for ( ii = 1 ; ii < n12 ; ii ++ ) {
	/* g_assert(!isnan(gamma[ii])) ; */
	gamma[ii] = sqrt(gamma[ii]*gamma[ii] +
			 c2Tbar[ii-1]*c2Tbar[ii-1] -
			 c2T[ii-1]*c2T[ii-1]) ;
	/* if ( isnan(gamma[ii]) ) */
	/*   g_error("NaN: ii==%d; c2T=%lg; c2Tbar=%lg; nu=%lg; mu=%lg; rho=%lg; " */
	/* 	  "c1T=%lg;", */
	/* 	  ii, c2T[ii-1], c2Tbar[ii-1], nu, mu, rho, c1T[ii-1]) ; */
      }

      ii = k - 1 ; 
      blaswrap_daxpy(ii,mu,b1,one,u1,one) ;
		     
      omega[k-1] = ga_bar ;
      for ( ii = 0 ; ii < k -1 ; ii ++ ) {
	omega[ii] = 1.0/sqrt(1.0/(omega[ii]*omega[ii]) +
			     u1[ii]*u1[ii]/ga_bar/ga_bar -
			     b1[ii]*b1[ii]/ga/ga) ;
      }
      if ( k < Rm ) {
	/*this is a Q update*/
	gdouble Q1, Q2 ;
	for ( ii = 0 ; ii < m ; ii ++ ) {
	  Q1 = Q[matrix_index(m,k,ii,k-1)] ;
	  Q2 = Q[matrix_index(m,k,ii,k  )] ;
	  Q[matrix_index(m,k,ii,k-1)] = (Q1*mu + Q2*nu)/rho ;
	  Q[matrix_index(m,k,ii,k  )] = (Q1*nu - Q2*mu)/rho ;
	}
      }

      /* fprintf(stderr, "end loop\n")  ; */
      /* return 0 ; */
    }

    tmp = G_MAXDOUBLE ;
    for ( ii = 0 ; ii < k - 1 ; ii ++ ) {
      if ( omega[ii] < tmp ) {
	tmp = omega[ii] ; i = ii ;
      }
    }
  
    if ( tmp > tol ) {
      /* return 0 ; */
      /* fprintf(stderr, "%d %lg\n", i, tmp) ; */
      break ;
    }

    g_assert_not_reached() ;
  
    if ( i < k - 1 ) {
      array_shuffle_int(pvt, k-1, i) ;
      column_shuffle(R, mr, k-1, i, mr) ;
      for ( ii = i ; ii < k-1 ; ii ++ ) {
	ta = R[matrix_index(mr,nr,ii  ,ii)] ;
	tb = R[matrix_index(mr,nr,ii+1,ii)] ;
	drotg_(&ta, &tb, &Gc, &Gs) ;
	if ( Gc*R[matrix_index(mr,nr,ii  ,ii)] +
	     Gs*R[matrix_index(nr,nr,ii+1,ii)] < 0.0 ) {
	  Gc = -Gc ; Gs = -Gs ;
	}
	drot_(&n,
	      &(R[matrix_index(mr,nr,ii  ,0)]), &mr,
	      &(R[matrix_index(mr,nr,ii+1,0)]), &mr,
	      &Gc, &Gs) ;
      }
      if ( R[matrix_index(mr,nr,k-1,k-1)] < 0 ) {
	for ( ii = 0 ; ii < n ; ii ++ )
	  R[matrix_index(mr,nr,k-1,ii)] *= -1 ;
	/*needs to be applied to Q as well*/
      }
      k -- ;    
    }
  
    nab = n - k ; mab = k ;
    for ( i = 0 ; i < mab ; i ++ ) {
      for ( j = 0 ; j < nab ; j ++ ) {
	AB[matrix_index(mab,nab,i,j)] = R[matrix_index(mr,nr,i,k+j)] ;
      }
    }
  
    dtrtrs_("U", "N", "N", &k, &nab, R, &mr, AB, &mab, &info) ;

    for ( j = 0 ; j < MIN(m,n)-k ; j ++ ) {
      for ( i = k ; i < k+j+1 ; i ++ ) {
	gamma[j] +=
	  R[matrix_index(mr,nr,i,k+j)]*R[matrix_index(mr,nr,i,k+j)] ;
      }
      gamma[j] = sqrt(gamma[j]) ;
    }
    for ( i = 0 ; i < k ; i ++ ) {
      for ( j = 0 ; j < k ; j ++ ) {
	work[matrix_index(k,k,i,j)] = R[matrix_index(mr,nr,i,j)] ;
      }
    }

    dtrtri_("U", "N", &k, work, &k, &info) ;
    for ( i = 0 ; i < k ; i ++ ) {
      omega[i] =
	1.0/blaswrap_dnrm2((k), &(work[matrix_index(k,k,i,0)]),
			   (k)) ;
    }
    fprintf(stderr, "loop end\n") ;
  }

  *rank = k ;
  *ldr = mr ;
  
  return 0 ;
}