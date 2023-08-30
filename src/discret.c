/* Copyright (C) 2007, 2020, 2023 by  Michael Carley */

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <blaswrap.h>

#include "gqr.h"
#include "gqr-private.h"

#ifdef HAVE_LIBRRQR
extern void dgeqpx_(gint *job, gint *m, gint *n, gint *k,
		    gdouble *A, gint *lda, gdouble *C, gint *ldc,
		    gint *jpvt, gdouble *ircond, gdouble *orcond,
		    gint *rank, gdouble *svlues, gdouble *work,
		    gint *lwork, gint *info) ;
#endif /*HAVE_LIBRRQR*/

/*FORTRAN matrix indexing*/
#define matrix_index(_m,_n,_i,_j) ((_i) + (_j)*(_m))

static gint compare_double(const void *a, const void *b)

{
  if ( *((gdouble *)a) < *((gdouble *)b)) return -1 ;
  if ( *((gdouble *)a) > *((gdouble *)b)) return  1 ;
  
  return 0 ;
}

static void legendre_fill(gint n, gdouble x, gdouble *P)

{
  gint i ;

  P[0] = 1.0 ; P[1] = x ;

  for ( i = 1 ; i < n ; i ++ ) {
    P[i+1] = (x*P[i]*(2*i+1) - P[i-1]*i)/(i+1) ;
  }
  
  return ;
}

gint gqr_discretize_orthogonalize(gdouble *A, gint m, gint n,
				  gdouble tol,
				  gint *rank, gint rankmax,
				  gdouble *Q, gdouble *R11,
				  gint *pvt, gint *ldr,
				  gdouble *work, gint lwork)

{
  gint k, lda, ldc, ma, na, info, job, i ;
  gdouble ircond, orcond, svlues[4] ;

  if ( lwork < (3*n+1) + rankmax )
    g_error("%s: workspace too small (%d < %d)",
	    __FUNCTION__, lwork, 3*n+1 + rankmax) ;
#ifndef HAVE_LIBRRQR
  gqr_srrqr(A, m, n, 1.0, tol, Q, R11, pvt, rank, ldr, work, lwork) ;

  return 0 ;
#else /*HAVE_LIBRRQR*/

  ma = m ; na = n ;
  lda = MAX(1, ma) ;
  ldc = ma ;
  memcpy(R11, A, ma*na*sizeof(gdouble)) ;

  job = 3 ;
  lda = ma ;
  k = ma ;
  ldc = k ;
  
  *ldr = ldc ;
  
  info = 0 ; ircond = tol ;
  memset(Q, 0, ldc*ldc*sizeof(gdouble)) ;
  for ( i = 0 ; i < ldc ; i ++ ) Q[i*ldc+i] = 1.0 ;
  dgeqpx_(&job, &ma, &na, &k, R11, &lda, Q, &ldc, pvt, &ircond, &orcond,
	  rank, svlues, work, &lwork, &info) ;  

  for ( i = 0 ; i < na ; i ++ ) pvt[i] -= 1 ;
#endif /*HAVE_LIBRRQR*/
  
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
