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
  Adaptive discretization of arbitrary functions
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

gdouble grule_bgr_func_scattering_r(gdouble t, gint i, gqr_parameter_t *p)

/*
 * quadratures for radial functions in Bremer and Gimbutas, On the
 * numerical evaluation of the singular integrals of scattering
 * theory.
 * 
 */
  
{
  gdouble d, f, x0, x1 ;

  g_assert(gqr_parameter_nf(p) > 3) ;
  d = gqr_parameter_double(p, 3) ;

  f = (1.0-d)*t + d ;
  
  if ( i == 0 ) return 1.0/f ;

  if ( i > 0 ) return pow(f, i-1) ;

  /*calculate integral of f over specified range*/
  x0  = gqr_parameter_double(p, 0) ;
  x1  = gqr_parameter_double(p, 1) ;

  i = -(i+1) ;

  if ( i == 0 ) {
    return (log((1.0-d)*x1+d) - log((1.0-d)*x0+d))/(1.0-d) ;
  }
  
  i -= 1 ;
  return (pow((1.0-d)*x1+d,i+1) - pow((1.0-d)*x0+d,i+1))/(i+1)/(1-d) ;
  
  return 0.0 ;
}

/*
 * implementation of Bremer, Gimbutas, Rokhlin, `A nonlinear
 * optimization procedure for generalized Gaussian quadratures', 
 * https://dx.doi.org/10.1137/080737046
 *
 */

gint grule_bgr(gdouble *x, gdouble *w, gqr_parameter_t *p)

/*
 * parameter p contains
 *
 * pointer 0: adapt_func for basis functions
 *
 * double  0: a, start of integration range
 * double  1: b,  of integration range
 * double  2: tol, discretization tolerance
 *
 * int     0: number of basis functions to include
 * int     1: nq, number of points in quadrature rule for discretization
 * int     2: nimax, maximum number of discretization intervals
 * int     3: rankmax, maximum length of output quadrature
 * 
 */
  
{
  gdouble x0, x1, tol, *ival, *Q, *R11, *u, *A, *work, *r, *z, xi, wi ;
  gint i, j, k, nf, nq, ni, nimax, rank, rankmax, *pvt, lwork, nn, nfunc ;
  gqr_rule_t *rule ;
  gqr_adapt_func_t func ;
  
  if ( gqr_parameter_np(p) < 1 ) 
    g_error("%s: at least one pointer must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_nf(p) < 3 ) 
    g_error("%s: at least three doubles must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_ni(p) < 4 ) 
    g_error("%s: at least four ints must be set in parameters",
	    __FUNCTION__) ;

  func = gqr_parameter_pointer(p, 0) ;
  
  x0  = gqr_parameter_double(p, 0) ;
  x1  = gqr_parameter_double(p, 1) ;
  tol = gqr_parameter_double(p, 2) ;

  nf      = gqr_parameter_int(p, 0) ;
  nq      = gqr_parameter_int(p, 1) ;
  nimax   = gqr_parameter_int(p, 2) ;
  rankmax = gqr_parameter_int(p, 3) ;

  rule = gqr_rule_alloc(nq) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, nq, NULL) ;
  
  ival = (gdouble *)g_malloc(nimax*sizeof(gdouble)) ;
  pvt = (gint *)g_malloc0(nimax*nq*sizeof(gint)) ;

  lwork = 3*nimax*nq+1 + rankmax ;
  work = (gdouble *)g_malloc(lwork*sizeof(gdouble)) ;
  
  ival[0] = x0 ; ival[1] = x1 ;
  
  Q   = (gdouble *)g_malloc(nq*nq*sizeof(gdouble)) ;
  R11 = (gdouble *)g_malloc(rankmax*rankmax*sizeof(gdouble)) ;
  r   = (gdouble *)g_malloc0(nf*sizeof(gdouble)) ;
  z   = (gdouble *)g_malloc0(rankmax*sizeof(gdouble)) ;

  gqr_discretize_interp(rule, Q) ;

  ni = 1 ;

  /*split interval for discretization*/
  for ( i = 0 ; i < nf ; i ++ ) {
    nn = 1 ;
    while ( nn != 0 ) {
      nn = gqr_discretize_adaptive(func, i, p, rule,
				   Q, ival, &ni, nimax, tol) ;
    }
  }

  /*number of points in discretized functions*/
  nfunc = nq*ni ;
  
  /*use FORTRAN indexing from here*/
  /*fill function array*/
  u = (gdouble *)g_malloc(nfunc*rankmax*sizeof(gdouble)) ;
  A = (gdouble *)g_malloc(nfunc*     nf*sizeof(gdouble)) ;
  /*fill columns of matrix A with functions*/
  for ( i = 0 ; i < nf ; i ++ ) {
    gqr_discretize_fill(func, i, p, rule, ival, ni, &(A[i*nfunc]), 1) ;
  }

  /*need to weight columns of A by quadrature weights in long rule*/
  for ( i = 0 ; i < ni ; i ++ ) {
    gdouble dx ;
    dx = 0.5*(ival[i+1] - ival[i]) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      for ( k = 0 ; k < nf ; k ++ ) {
	A[k*nfunc+i*nq+j] *= sqrt(gqr_rule_weight(rule, j)*dx) ;
      }
    }
  }

  gqr_discretize_orthogonalize(A, nfunc, nf, tol, &rank, rankmax,
			       u, R11, pvt, work, lwork) ;

  /*generate a quadrature rule*/
  for ( i = 0 ; i < ni ; i ++ ) {
    for ( j = 0 ; j < nq ; j ++ ) {
      gqr_discretize_quadrature(rule, ival, ni, i, j, &xi, &wi) ;
      for ( k = 0 ; k < nf ; k ++ ) {
	r[k] += u[k*nfunc+i*nq+j]*sqrt(wi) ;
      }
    }
  }

  /*there should be some way to do the RRQR on the transpose without
    explicitly transposing, but I haven't found it yet*/
  for ( i = 0 ; i < rank ; i ++ ) {
    for ( j = 0 ; j < nfunc ; j ++ ) {
      A[j*rank+i] = u[i*nfunc+j] ;
    }
  }

  memset(pvt, 0, nfunc*sizeof(gint)) ;
  gqr_discretize_orthogonalize(A, rank, nfunc, 1e-15, &rank, rankmax,
  			       u, R11, pvt, work, lwork) ;

  for ( i = 0 ; i < rank ; i ++ ) {
    /* z[i] = 0.0 ; */
    for ( j = 0 ; j < rank ; j ++ ) {
      /*this is the transpose multiplication in FORTRAN indexing*/
      z[i] += u[i*rank+j]*r[j] ;
    }
  }

  i = rank-1 ;
  z[rank-1] = z[rank-1]/R11[i*rank+i] ;
  for ( i = rank-2 ; i >= 0 ; i -- ) {
    for ( j = i+1 ; j < rank ; j ++ ) {
      z[i] -= R11[j*rank+i]*z[j] ;
    }
    z[i] /= R11[i*rank+i] ;
  }

  for ( i = 0 ; i < rank ; i ++ ) {
    /*subtract one to switch from FORTRAN indexing*/
    gqr_discretize_quadrature(rule, ival, ni, pvt[i]-1, -1, &xi, &wi) ;
    x[i] = xi ; w[i] = sqrt(wi)*z[i] ;
  }
  
  g_free(Q) ;
  g_free(R11) ;
  g_free(u) ;
  g_free(r) ;
  g_free(z) ;
  g_free(A) ;
  g_free(ival) ;
  g_free(pvt) ;
  
  return rank ;
}

gint gqr_rule_bgr_check(gqr_rule_t *rule, gqr_parameter_t *p,
			gint *imax, gdouble *emax,
			FILE *output)

{
  gdouble x0, x1, dx, xb, x, f, Iq[256], Ia[256] ;
  gint i, j, nf ;
  gqr_adapt_func_t func ;
  
  if ( gqr_parameter_np(p) < 1 ) 
    g_error("%s: at least one pointer must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_nf(p) < 3 ) 
    g_error("%s: at least three doubles must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_ni(p) < 4 ) 
    g_error("%s: at least four ints must be set in parameters",
	    __FUNCTION__) ;

  func = gqr_parameter_pointer(p, 0) ;
  
  x0 = gqr_parameter_double(p, 0) ;
  x1 = gqr_parameter_double(p, 1) ;

  nf = gqr_parameter_int(p, 0) ;

  g_assert(nf < 256) ;
  
  xb = gqr_rule_scale(rule, rule->a, rule->b, &xb, &dx) ;
  *imax = 0 ; *emax = 0.0 ;
  for ( i = 0 ; i < nf ; i ++ ) {
    Iq[i] = 0.0 ;
    for ( j = 0 ; j < gqr_rule_length(rule) ; j ++ ) {
      x = xb + dx*gqr_rule_abscissa(rule, j) ;
      f = func(x, i, p) ;
      Iq[i] += f*gqr_rule_weight(rule, j)*dx ;
    }

    Ia[i] = func(x, -i-1, p) ;
    if ( fabs(Ia[i] - Iq[i]) > *emax ) {
      *emax = fabs(Ia[i] - Iq[i]) ;
      *imax = i ;
    }
  }

  if ( output != NULL ) {
    for ( i = 0 ; i < nf ; i ++ ) 
      fprintf(output, "%d %lg %lg %lg\n",
	      i, Iq[i], Ia[i], fabs(Ia[i]-Iq[i])) ;
  }
  
  return 0 ;
}
