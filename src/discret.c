/* Copyright (C) 2007, 2020 by  Michael Carley */

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
  Adaptive discretization of arbitrary functions
*/

#include <stdio.h>
#include <math.h>

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
gint dlarfb_(gchar *side, gchar *trans, gchar *direct, gchar *storev,
	     gint *m, gint *n, gint *k, gdouble *v, gint *ldv,
	     gdouble *t, gint *ldt, gdouble *c, gint *ldc,
	     gdouble *work, gint *ldwork) ;
gint dlarft_(gchar *direct, gchar *storev, gint *n, gint *k,
	     gdouble *v, gint *ldv, gdouble *tau, gdouble *t, gint *ldt) ;

/*FORTRAN matrix indexing*/
#define matrix_index(_m,_n,_i,_j) ((_i) + (_j)*(_m))

static gint compare_double(const void *a, const void *b)

{
  if ( *((gdouble *)a) < *((gdouble *)b)) return -1 ;
  if ( *((gdouble *)a) > *((gdouble *)b)) return  1 ;
  
  return 0 ;
}

static gint legendre_fill(gint n, gdouble x, gdouble *P)

{
  gint i ;

  P[0] = 1.0 ; P[1] = x ;

  for ( i = 1 ; i < n ; i ++ ) {
    P[i+1] = (x*P[i]*(2*i+1) - P[i-1]*i)/(i+1) ;
  }
  
  return 0 ;
}

gint rrqr(gdouble *A, gint m, gint n, gdouble *tau, gint *jpvt,
	  gdouble *work, gint lwork)

{
  gint info ;
  /**
   * perform rrqr on A using LAPACK function A is supplied in FORTRAN
   * ordering, so vectors lie on rows in C ordering
   */

  if ( lwork < 3*n+1 )
    g_error("%s: workspace too small (%d < %d)",
	    __FUNCTION__, lwork, 3*n+1) ;
  dgeqp3_(&m, &n, A, &m, jpvt, tau, work, &lwork, &info) ;

  g_assert(info == 0) ;
  
  return 0 ;
}

gint rrqr_rank(gdouble *R, gint m, gint n, gdouble ee)

{
  gint rank ;

  for ( rank = 1 ; rank < MIN(m,n) ; rank ++ ) {
    /* if ( fabs(R[matrix_index(m,n,rank,rank)]) < */
    /* 	 fabs(R[matrix_index(m,n,   0,   0)])*ee) */
    if ( fabs(R[matrix_index(m,n,rank,rank)]) < ee )
      return rank ;
  }
  
  return rank ;
}

gint rrqr_qr(gdouble *A, gint m, gint n, gdouble *tau, gint rank,
	     gdouble *Q, gdouble *R11, gdouble *work, gint lwork)

{
  gint i, j, ldwork ;
  gdouble T[8192] ;
  
  for ( j = 0 ; j < rank ; j ++ ) {
    for ( i = 0 ; i < MIN(m,rank) ; i ++ ) {
      Q[matrix_index(m, rank, i, j)] = 0.0 ;
    }
    Q[matrix_index(m, rank, j, j)] = 1.0 ;
    for ( i = 0 ; i <= j ; i ++ ) {
      R11[matrix_index(rank,rank,i,j)] = A[matrix_index(m,n,i,j)] ;
    }
    for ( i = j+1 ; i < rank ; i ++ ) {
      R11[matrix_index(rank,rank,i,j)] = 0.0 ;
    }
  }

  /*this is correct for "L"*/
  ldwork = rank ;
  dlarft_("F", "C", &rank, &rank, A, &m, tau, T, &rank) ;
  dlarfb_("L", "N", "F", "C", &m, &rank, &rank, A, &m, T, &rank, Q, &m, work,
	  &ldwork) ;
  
  /*R11 is rank x rank; Q is m x rank*/
  /*make diagonals of R11 positive: multiply columns of R11 by sign of
    diagonal entry, multiply rows of Q*/
  for ( i = 0 ; i < rank ; i ++ ) {
    if ( R11[matrix_index(rank,rank,i,i)] < 0.0 ) {
      for ( j = 0 ; j < rank ; j ++ ) {
	R11[matrix_index(rank,rank,i,j)] *= -1 ;
      }
      for ( j = 0 ; j < m ; j ++ ) {
	Q[matrix_index(m,rank,j,i)] *= -1 ;
      }
    }
  }
  
  return 0 ;
}

gint gqr_discretize_orthogonalize(gdouble *A, gint m, gint n,
				  gdouble tol,
				  gint *rank, gint rankmax,
				  gdouble *Q, gdouble *R11,
				  gint *pvt, gint *ldr,
				  gdouble *work, gint lwork)

{
  gdouble *tau ;

  if ( lwork < (3*n+1) + rankmax )
    g_error("%s: workspace too small (%d < %d)",
	    __FUNCTION__, lwork, 3*n+1 + rankmax) ;


  gqr_srrqr(A, m, n, 1.0, tol, Q, R11, pvt, rank, ldr, work, lwork) ;

  return 0 ;
  
  tau = &(work[3*n+1]) ;
  
  rrqr(A, m, n, tau, pvt, work, 3*n+1) ;
  *rank = rrqr_rank(A, m, n, tol) ;
  if ( *rank > rankmax ) *rank = rankmax ;    
  rrqr_qr(A, m, n, tau, *rank, Q, R11, work, lwork) ;
  
  return 0 ;
}

gint gqr_discretize_adaptive(gqr_adapt_func_t func, gint idx, gpointer data,
			     gqr_rule_t *rule, gdouble *Q,
			     gdouble *ival, gint *ni, gint nimax, 
			     gdouble tol)

/*
  Given function func f_i(x) refine the division intervals until the
  required tolerance is achieved
*/

{
  gint nn, i, j, nq, one = 1 ;
  gdouble f[128], cfft[128], x, xb, dx, norm ;
  gdouble al = 1.0, bt = 0.0 ;
  
  nq = gqr_rule_length(rule) ;

  nn = 0 ;
  for ( i = 0 ; (i < *ni) && (*ni+nn < nimax) ; i ++ ) {
    /*check the interval*/
    dx = 0.5*(ival[i+1] - ival[i]) ;
    xb = 0.5*(ival[i+1] + ival[i]) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      x = xb + dx*gqr_rule_abscissa(rule, j) ;
      f[j] = func(x, idx, data) ;
    }
    blaswrap_dgemv(FALSE, nq, nq, al, Q, nq, f, one, bt, cfft, one) ;
    norm = 0.0 ;
    for ( j = nq/2 ; j < nq ; j ++ ) norm += cfft[j]*cfft[j] ;
    if ( norm > tol*tol ) {
      /*split the interval*/
      ival[*ni+1+nn] = xb ;
      nn ++ ;
    }
  }

  *ni += nn ;

  qsort(ival, *ni+1, sizeof(gdouble), compare_double) ;

  /*check for duplicates*/
  for ( i = 0 ; i < *ni ; i ++ ) {
    g_assert(ival[i] < ival[i+1]) ;
  }
  
  return nn ;
}

gint gqr_discretize_interp(gqr_rule_t *rule, gdouble *Q)

{
  gint i, j, nq ;
  gdouble x, w ;
  
  nq = gqr_rule_length(rule) ;

  for ( j = 0 ; j < nq ; j ++ ) {
    w = gqr_rule_weight(rule, j) ;
    x = gqr_rule_abscissa(rule, j) ;
    Q[0*nq+j] = 1.0*w*0.5 ;
    Q[1*nq+j] = x*  w*1.5 ;
  }

  /*recursion incorporating the normalization to give orthonormal
    basis*/
  for ( i = 1 ; i < nq-1 ; i ++ ) {
    for ( j = 0 ; j < nq ; j ++ ) {
      x = gqr_rule_abscissa(rule, j) ;
      Q[(i+1)*nq+j] = (x*Q[i*nq+j] - Q[(i-1)*nq+j]*i/(2*i-1))*(2*i+3)/(i+1) ;
    }
  }
  
  return 0 ;
}

static gint point_interval(gdouble *ival, gint ni, gdouble x)

{
  gint i0, i1, i ;

  i0 = 0 ; i1 = ni ;

  while ( i0 <= i1 ) {
    i = (i0+i1)/2 ;
    if ( ival[i] < x ) {
      i0 = i+1 ;
    } else {
      if ( ival[i] > x ) {
	i1 = i-1 ;
      } else {
	return i ;
      }
    }
  }
  
  return i1 ;
}

gint gqr_discretize_eval(gqr_rule_t *rule, gdouble *Q,
			 gdouble *ival, gint ni,
			 gdouble *u, gint nu, gint ustr,
			 gdouble x, gdouble *f)

{
  gint i, nq, one = 1 ;
  gdouble cfft[64], Pn[128], t ;
  gdouble al = 1.0, bt = 0.0 ;
  
  /*evaluation point outside domain*/
  if ( x < ival[0] || x > ival[ni] ) return -1 ;

  nq = gqr_rule_length(rule) ;
  i = point_interval(ival, ni+1, x) ;

  blaswrap_dgemv(FALSE, nq, nq, al, Q, nq, &(u[i*nq*ustr]), ustr, bt,
		 cfft, one) ;

  t = (2.0*x - ival[i+1] - ival[i])/(ival[i+1] - ival[i]) ;
  legendre_fill(nq, t, Pn) ;

  f[0] = 0.0 ;
  for ( i = 0 ; i < nq ; i ++ ) f[0] += cfft[i]*Pn[i] ;
  
  return 0 ;
}

gint gqr_discretize_fill(gqr_adapt_func_t func, gint idx, gpointer data,
			 gqr_rule_t *rule, gdouble *ival, gint ni,
			 gdouble *u, gint ustr)

{
  gint i, j, nq ;
  gdouble x, xb, dx ;
  
  nq = gqr_rule_length(rule) ;
  
  for ( i = 0 ; i < ni ; i ++ ) {
    dx = 0.5*(ival[i+1] - ival[i]) ;
    xb = 0.5*(ival[i+1] + ival[i]) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      x = xb + dx*gqr_rule_abscissa(rule, j) ;
      u[(i*nq+j)*ustr] = func(x, idx, data) ;
    }
  }

  return 0 ;
}

gint gqr_discretize_quadrature(gqr_rule_t *rule,
			       gdouble *ival, gint ni, gint i, gint j,
			       gdouble *x, gdouble *w)

{
  gdouble dx, xb ;

  /*for negative j, treat i as absolute index*/
  if ( j < 0 ) {
    j = i % gqr_rule_length(rule) ;
    i = i/gqr_rule_length(rule) ;
  }
  
  g_assert(i < ni && i >= 0) ;
  g_assert(j < gqr_rule_length(rule) && j >= 0) ;

  dx = 0.5*(ival[i+1] - ival[i]) ;
  xb = 0.5*(ival[i+1] + ival[i]) ;

  *x = dx*gqr_rule_abscissa(rule, j) + xb ;
  *w = dx*gqr_rule_weight(rule, j) ;

  return 0 ;
}
