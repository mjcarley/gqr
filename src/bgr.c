/* Copyright (C) 2020, 2023 by Michael Carley */

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

#include "gqr.h"
#include "gqr-private.h"

#define matrix_index(_m,_n,_i,_j) ((_i) + (_j)*(_m))

/**
 * @defgroup bgr Generalized quadrature rules
 * @{
 * 
 */

static gint compare_double2(const void *a, const void *b)

{
  if ( *((gdouble *)a) < *((gdouble *)b)) return -1 ;
  if ( *((gdouble *)a) > *((gdouble *)b)) return  1 ;
  
  return 0 ;
}

/* static gint print_matrix(gdouble *A, gint m, gint n, FILE *out) */

/* { */
/*   gint i, j ; */

/*   for ( i = 0 ; i < m ; i ++ ) { */
/*     for ( j = 0 ; j < n ; j ++ ) */
/*       fprintf(out, "%1.16e ", A[matrix_index(m,n,i,j)]) ; */
/*     fprintf(out, "\n") ; */
/*   } */
  
/*   return 0 ; */
/* } */

gdouble grule_bgr_func_scattering_r(gdouble t, gint i, gqr_parameter_t *p)

/*
 * quadratures for radial functions in Bremer and Gimbutas, On the
 * numerical evaluation of the singular integrals of scattering
 * theory, http://dx.doi.org/10.1016/j.jcp.2013.05.048
 * 
 */
  
{
  gdouble d, f, x0, x1 ;

  g_assert(gqr_parameter_double_number(p) > 3) ;
  d = gqr_parameter_double(p, 3) ;

  g_assert(d > 0.0) ;
  
  f = (1.0-d)*t + d ;
  
  if ( i == 0 ) return 1.0/f ;

  if ( i > 0 ) {
    f = pow(f, i-1) ;
    return f ;
  }

  /*calculate integral of f over specified range*/
  x0  = gqr_parameter_double(p, 0) ;
  x1  = gqr_parameter_double(p, 1) ;

  i = -(i+1) ;

  if ( i == 0 ) {
    return (log((1.0-d)*x1+d) - log((1.0-d)*x0+d))/(1.0-d) ;
  }
  
  i -= 1 ;
  return (pow((1.0-d)*x1+d,i+1) -
	  pow((1.0-d)*x0+d,i+1))/(gdouble)(i+1)/(1.0-d) ;
  
  return 0.0 ;
}

gdouble grule_bgr_func_scattering_range_r(gdouble t, gint idx,
					  gqr_parameter_t *p)

/*
 * quadratures for radial functions in Bremer and Gimbutas, On the
 * numerical evaluation of the singular integrals of scattering
 * theory, http://dx.doi.org/10.1016/j.jcp.2013.05.048
 */
  
{
  gdouble d0, d1, d, f, x0, x1 ;
  gint nd, i, j ;
  gqr_rule_t *rule ;
  
  g_assert(gqr_parameter_int_number(p) > 4) ;
  nd = gqr_parameter_int(p, 4) ;

  if ( p->rules[0] == NULL ) p->rules[0] = gqr_rule_alloc(4*nd) ;
  rule = p->rules[0] ;

  if ( gqr_rule_length(rule) != nd ) {
    gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, nd, NULL) ;
  }

  g_assert(gqr_parameter_double_number(p) > 4) ;
  d0 = gqr_parameter_double(p, 3) ;
  d1 = gqr_parameter_double(p, 4) ;

  g_assert(d0 >= 0.0) ;
  g_assert(d1 > d0) ;

  if ( idx >= 0 ) {
    /*first nd functions are the singular ones*/

    if ( idx < nd ) {
      j = idx ;
      d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, j) ;
      f = (1.0-d)*t + d ;
      f = 1.0/f ;

      return f ;
    }

    /*powers of (1-d)*t + d are spanned by monomials so we only need
      to include one of them*/
    j = idx - nd ;
    d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, 0) ;
    f = (1.0-d)*t + d ;
    f = pow(f, j) ;

    return f ;
  }

  /*calculate integral of f over specified range*/
  x0  = gqr_parameter_double(p, 0) ;
  x1  = gqr_parameter_double(p, 1) ;

  idx = -(idx+1) ;

  if ( idx < nd ) {
    j = idx ;
    d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, j) ;
    return (log((1.0-d)*x1+d) - log((1.0-d)*x0+d))/(1.0-d) ;
  }

  j = idx - nd ;
  d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, 0) ;
  return (pow((1.0-d)*x1+d,j+1) -
	  pow((1.0-d)*x0+d,j+1))/(gdouble)(j+1)/(1.0-d) ;
    
  if ( idx >= 0 ) {
    /*select value of d in range*/
    i = idx/nd ;
    j = idx % nd ;
  
    d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, j) ;
  
    f = (1.0-d)*t + d ;
  
    if ( i == 0 ) return 1.0/f ;

    if ( i > 0 ) {
      f = pow(f, i-1) ;
      return f ;
    }
  }

  /*calculate integral of f over specified range*/
  x0  = gqr_parameter_double(p, 0) ;
  x1  = gqr_parameter_double(p, 1) ;

  idx = -(idx+1) ;

  i = idx/nd ;
  j = idx % nd ;
  
  d = 0.5*(d1 + d0) + 0.5*(d1 - d0)*gqr_rule_abscissa(rule, j) ;
  
  if ( i == 0 ) {
    return (log((1.0-d)*x1+d) - log((1.0-d)*x0+d))/(1.0-d) ;
  }
  
  i -= 1 ;
  return (pow((1.0-d)*x1+d,i+1) -
	  pow((1.0-d)*x0+d,i+1))/(gdouble)(i+1)/(1.0-d) ;
}

gdouble grule_bgr_func_scattering_range_th(gdouble t, gint idx,
					   gqr_parameter_t *p)

/*
 * quadratures for angular functions in Bremer and Gimbutas, On the
 * numerical evaluation of the singular integrals of scattering
 * theory, http://dx.doi.org/10.1016/j.jcp.2013.05.048
 * at range of values of r_0, \theta_0
 * 
 * indexing of functions is:
 * idx = 0, nj*(nj+3)*3/2 M_j(\theta_{0}u)\cos(i\theta_{0}u), i=0,3*(j+1)+2
 * idx -= nj*(nj+3)*3/2
 * idx = 0, nj*(3*nj+7)/2 M_j(\theta_{0}u)\sin(i\theta_{0}u), i=1,3*(j+1)+2
 * idx -= nj*(3*nj+7)/2
 * idx = 0, nk \cos(k\theta_{0}u), k=0, nk-1
 * idx -= nk
 * idx = 0, nk-1 \sin(k\theta_{0}u), k=1, nk-1
 * 
 */
  
{
  gdouble r0, th0, th, f, M, thmin, thmax, rmin, rmax ;
  gint nj, nk, i, j, k, off, nr, nt, nrt, tr, tt ;
  gqr_rule_t *rrule, *trule ;
  
  g_assert(gqr_parameter_int_number(p) > 7) ;
  nj = gqr_parameter_int(p, 4) ;
  nk = gqr_parameter_int(p, 5) ;
  nr = gqr_parameter_int(p, 6) ;
  nt = gqr_parameter_int(p, 7) ;

  g_assert(gqr_parameter_double_number(p) > 6) ;
  rmin = gqr_parameter_double(p, 3) ;
  rmax = gqr_parameter_double(p, 4) ;
  thmin = gqr_parameter_double(p, 5) ;
  thmax = gqr_parameter_double(p, 6) ;

  g_assert(rmax > rmin) ; 
  g_assert(thmax > thmin) ; 
  
  if ( p->rules[0] == NULL ) p->rules[0] = gqr_rule_alloc(4*nr) ;
  if ( p->rules[1] == NULL ) p->rules[1] = gqr_rule_alloc(4*nt) ;

  rrule = p->rules[0] ;
  trule = p->rules[1] ;

  if ( gqr_rule_length(rrule) != nr ) {
    gqr_rule_select(rrule, GQR_GAUSS_LEGENDRE, nr, NULL) ;
  }

  if ( gqr_rule_length(trule) != nt ) {
    gqr_rule_select(trule, GQR_GAUSS_LEGENDRE, nt, NULL) ;
  }

  /*total number of entries per (r_0,\theta_0) pair*/
  nrt = nj*(nj+3)*3/2 + nj*(3*nj+7)/2 + nk + nk - 1 ;

  i = idx/nrt ;
  tr = i/nt ;
  tt = i % nt ;

  r0 = 0.5*(rmax + rmin) +
    0.5*(rmax - rmin)*gqr_rule_abscissa(rrule, tr) ;
  th0 = 0.5*(thmax + thmin) +
    0.5*(thmax - thmin)*gqr_rule_abscissa(trule, tt) ;
  
  g_assert(r0 > 0.0) ;
  g_assert(r0 < 1.0) ;

  /*set idx in the range of values for a single (r_0,\theta_0) pair*/
  idx = idx % nrt ;
  
  off = 0 ; 

  th = th0*t ;
  M = r0*sin(th0)/(r0*sin(th0-th) + sin(th)) ;

  if ( idx < nj*(nj+3)*3/2 ) {
    /*M_j \sin*/
    for ( j = 0 ; j < nj ; j ++ ) {
      if ( j*(j+3)*3/2 > idx ) break ;
    }

    j -- ;

    i = idx - j*(j+3)*3/2 ;

    if ( j == 0 ) M = log(M) ;
    else M = pow(M, j) ;

    f = M*cos(i*th) ;
    g_assert(!isnan(f)) ;
    return f ;
  }

  off += nj*(nj+3)*3/2 ;
  idx -= nj*(nj+3)*3/2 ;

  if ( idx < nj*(3*nj+7)/2 ) {
    /*M_j \sin*/
    for ( j = 0 ; j < nj ; j ++ ) {
      if ( j*(3*j+7)/2 > idx ) break ;
    }

    j -- ;

    i = idx - j*(3*j+7)/2 ;
    i ++ ; /*don't use the \sin i\theta term for i==0*/

    if ( j == 0 ) M = log(M) ;
    else M = pow(M, j) ;

    f = M*cos(i*th) ;
    g_assert(!isnan(f)) ;

    return f ;
  }

  off += nj*(3*nj+7)/2 ;
  idx -= nj*(3*nj+7)/2 ;

  if ( idx < nk ) {
    /*cos k\theta*/

    k = idx ;

    f = cos(k*th) ;
    g_assert(!isnan(f)) ;

    return f ;
  }

  off += nk ;
  idx -= nk ;

  if ( idx < nk ) {
    /*sin k\theta*/

    k = idx + 1 ;
    f = sin(k*th) ;

    return f ;
  }

  g_error("%s: function index (%d) out of range (nj=%d, nk=%d)",
	  __FUNCTION__, idx+off, nj, nk) ;

  return 0.0 ;
}


gdouble grule_bgr_func_scattering_th(gdouble t, gint idx, gqr_parameter_t *p)

/*
 * quadratures for angular functions in Bremer and Gimbutas, On the
 * numerical evaluation of the singular integrals of scattering
 * theory, http://dx.doi.org/10.1016/j.jcp.2013.05.048
 * 
 * indexing of functions is:
 * idx = 0, nj*(nj+3)*3/2 M_j(\theta_{0}u)\cos(i\theta_{0}u), i=0,3*(j+1)+2
 * idx -= nj*(nj+3)*3/2
 * idx = 0, nj*(3*nj+7)/2 M_j(\theta_{0}u)\sin(i\theta_{0}u), i=1,3*(j+1)+2
 * 
 */
  
{
  gdouble r0, th0, th, f, M ;
  gint nj, nk, i, j, k, off ;
  
  g_assert(gqr_parameter_int_number(p) > 5) ;
  nj = gqr_parameter_int(p, 4) ;
  nk = gqr_parameter_int(p, 5) ;

  g_assert(gqr_parameter_double_number(p) > 4) ;
  r0 = gqr_parameter_double(p, 3) ;
  th0 = gqr_parameter_double(p, 4) ;

  g_assert(r0 > 0.0) ;
  g_assert(r0 < 1.0) ;

  off = 0 ; 

  th = th0*t ;
  M = r0*sin(th0)/(r0*sin(th0-th) + sin(th)) ;

  if ( idx < nj*(nj+3)*3/2 ) {
    /*M_j \sin*/
    for ( j = 0 ; j < nj ; j ++ ) {
      if ( j*(j+3)*3/2 > idx ) break ;
    }

    j -- ;

    i = idx - j*(j+3)*3/2 ;

    if ( j == 0 ) M = log(M) ;
    else M = pow(M, j) ;

    f = M*cos(i*th) ;
    return f ;
  }

  off += nj*(nj+3)*3/2 ;
  idx -= nj*(nj+3)*3/2 ;

  if ( idx < nj*(3*nj+7)/2 ) {
    /*M_j \sin*/
    for ( j = 0 ; j < nj ; j ++ ) {
      if ( j*(3*j+7)/2 > idx ) break ;
    }

    j -- ;

    i = idx - j*(3*j+7)/2 ;
    i ++ ; /*don't use the \sin i\theta term for i==0*/

    if ( j == 0 ) M = log(M) ;
    else M = pow(M, j) ;

    f = M*cos(i*th) ;
    return f ;
  }

  off += nj*(3*nj+7)/2 ;
  idx -= nj*(3*nj+7)/2 ;

  if ( idx < nk ) {
    /*cos k\theta*/

    k = idx ;

    f = cos(k*th) ;
    return f ;
  }

  off += nk ;
  idx -= nk ;

  if ( idx < nk ) {
    /*sin k\theta*/

    k = idx + 1 ;
    f = cos(k*th) ;

    return f ;
  }

  g_error("%s: function index (%d) out of range (nj=%d, nk=%d)",
	  __FUNCTION__, idx+off, nj, nk) ;

  return 0.0 ;
}


/** 
 * Generate a generalized quadrature rule using the method of Bremer,
 * Gimbutas, and Rokhlin, 2010, `A nonlinear optimization procedure
 * for generalized Gaussian quadratures',
 * https://dx.doi.org/10.1137/080737046
 * 
 * @param x on output contains nodes of quadrature rule;
 * @param w on output contains weights of quadrature rule;
 * @param p parameters of basis functions of quadrature, including a pointer
 * to the ::gqr_adapt_func_t which evaluates the functions.
 * 
 * @return number of nodes in rule.
 */

gint grule_bgr(gdouble *x, gdouble *w, gqr_parameter_t *p)

/*
 * parameter p contains
 *
 * pointer 0: adapt_func for basis functions
 *
 * double  0: a, start of integration range
 * double  1: b, end of integration range
 * double  2: tol, discretization tolerance
 * double  3 onwards, parameters to be supplied to adapt_func
 *
 * int     0: number of basis functions to include
 * int     1: nq, number of points in quadrature rule for discretization
 * int     2: nimax, maximum number of discretization intervals
 * int     3: rankmax, maximum length of output quadrature
 * int     r onwards, parameters to be supplied to adapt_func
 * 
 */
  
{
  gdouble x0, x1, tol, *ival, *Qi, *R11, *u, *A, *work, *r, *z, xi, wi ;
  gint i, j, k, nf, nq, ni, nimax, rank, rankmax, *pvt, lwork, nn, nfunc ;
  gint ldr, nu ;
  gqr_rule_t *rule ;
  gqr_adapt_func_t func ;
  gdouble test0 ;
  
  if ( gqr_parameter_pointer_number(p) < 1 ) 
    g_error("%s: at least one pointer must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_double_number(p) < 3 ) 
    g_error("%s: at least three doubles must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_int_number(p) < 4 ) 
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

  nimax += 4 ;
  ival = (gdouble *)g_malloc(nimax*sizeof(gdouble)) ;
  pvt = (gint *)g_malloc0(nimax*nq*sizeof(gint)) ;

  nimax -= 4 ;
  
  ival[0] = x0 ; ival[1] = x1 ;
  
  Qi = (gdouble *)g_malloc(nq*nq*sizeof(gdouble)) ;
  gqr_discretize_interp(rule, Qi) ;

  ni = 1 ;

  /*split interval for discretization*/
  for ( i = 0 ; i < nf ; i ++ ) {
    nn = 1 ;
    while ( nn != 0 ) {
      nn = gqr_discretize_adaptive(func, i, p, rule,
				   Qi, ival, &ni, nimax, tol) ;
    }
  }

  /*number of points in discretized functions*/
  nfunc = nq*ni ;
  
  lwork = 3*nimax*nq+1 + 8*rankmax*nfunc ;
  work = (gdouble *)g_malloc(lwork*sizeof(gdouble)) ;

  R11 = (gdouble *)g_malloc0(nfunc*(MAX(rankmax,nf)+1)*sizeof(gdouble)) ;
  r   = (gdouble *)g_malloc0(nf*sizeof(gdouble)) ;
  z   = (gdouble *)g_malloc0(rankmax*sizeof(gdouble)) ;

  /*use FORTRAN indexing from here*/
  /*fill function array*/
  nu = MAX(nfunc, nf) ;
  fprintf(stderr, "nfunc = %d; nf = %d\n", nfunc, nf) ;
  u = (gdouble *)g_malloc0(nu*nu*sizeof(gdouble)) ;
  g_assert(u != NULL) ;
  /* u = (gdouble *)g_malloc0(nfunc*nf*sizeof(gdouble)) ; */
  A = (gdouble *)g_malloc0(nfunc*(MAX(rankmax,nf)+1)*sizeof(gdouble)) ;
  /*fill columns of matrix A with functions*/
  for ( i = 0 ; i < nf ; i ++ ) {
    gqr_discretize_fill(func, i, p, rule, ival, ni, &(A[i*nfunc]), 1) ;
  }

  /*need to weight columns of A by quadrature weights in long rule*/
  test0 = 0.0 ;
  for ( i = 0 ; i < ni ; i ++ ) {
    gdouble dx ;
    dx = 0.5*(ival[i+1] - ival[i]) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      for ( k = 0 ; k < nf ; k ++ ) {
	/*A_jk = \phi_{k}(x_j)sqrt(w_j)*/
	A[matrix_index(nfunc,nf,i*nq+j,k)] *=
	  sqrt(gqr_rule_weight(rule, j)*dx) ;
	test0 += gqr_rule_weight(rule, j)*dx ;
      }
    }
  }

  gqr_discretize_orthogonalize(A, nfunc, nf, tol, &rank, rankmax,
  			       u, R11, pvt, &ldr, work, lwork) ;

  /*generate a quadrature rule*/
  for ( i = 0 ; i < ni ; i ++ ) {
    for ( j = 0 ; j < nq ; j ++ ) {
      gqr_discretize_quadrature(rule, ival, ni, i, j, &xi, &wi) ;
      for ( k = 0 ; k < nf ; k ++ ) {
	r[k] += u[matrix_index(nfunc,nf,i*nq+j,k)]*sqrt(wi) ;
      }
    }
  }

  /*there should be some way to do the RRQR on the transpose without
    explicitly transposing, but I haven't found it yet*/
  for ( i = 0 ; i < rank ; i ++ ) {
    for ( j = 0 ; j < nfunc ; j ++ ) {
      A[matrix_index(rank,nfunc,i,j)] =	u[matrix_index(nfunc,rank,j,i)] ;
    }
  }

  fprintf(stderr, "rank == %d (max = %d)\n", rank, rankmax) ;
  g_assert(rank <= rankmax) ;

  memset(pvt, 0, nfunc*sizeof(gint)) ;
  memset(u, 0, nfunc*nf*sizeof(gdouble)) ;
  memset(R11, 0, nfunc*nf*sizeof(gdouble)) ;
  memset(work, 0, lwork*sizeof(gdouble)) ;
  gqr_discretize_orthogonalize(A, rank, nfunc, tol, &rank, rankmax,
  			       u, R11, pvt, &ldr, work, lwork) ;

  /* fprintf(stderr, "rank == %d\n", rank) ; */
  
  for ( i = 0 ; i < rank ; i ++ ) {
    z[i] = 0.0 ;
    for ( j = 0 ; j < rank ; j ++ ) {
      /*this is the transpose multiplication in FORTRAN indexing*/
      z[i] += u[i*rank+j]*r[j] ;
    }
  }

  i = rank-1 ;
  z[rank-1] = z[rank-1]/R11[i*ldr+i] ;
  for ( i = rank-2 ; i >= 0 ; i -- ) {
    for ( j = i+1 ; j < rank ; j ++ ) {
      z[i] -= R11[j*ldr+i]*z[j] ;
    }
    z[i] /= R11[i*ldr+i] ;
  }

  for ( i = 0 ; i < rank ; i ++ ) {
    gqr_discretize_quadrature(rule, ival, ni, pvt[i], -1, &xi, &wi) ;
    work[2*i+0] = xi ; work[2*i+1] = sqrt(wi)*z[i] ;
  }

  qsort(work, rank, 2*sizeof(gdouble), compare_double2) ;
  for ( i = 0 ; i < rank ; i ++ ) {
    x[i] = work[2*i+0] ; w[i] = work[2*i+1] ;
  }
  
  return rank ;
}

/** 
 * Estimate error in generalized quadrature rules by integrating basis
 * functions and checking against evaluation using quadrature rules.
 * 
 * @param rule quadrature rule to check;
 * @param p parameters used to generate basis functions;
 * @param imax on exit, index of basis function with greatest error in 
 * integral; 
 * @param emax on exit, greatest error in integral evaluation;
 * @param analytic reserved for future use;
 * @param output if not NULL, file stream to write diagnostic information
 * and error estimates.
 * 
 * @return 0 on success.
 */

gint gqr_rule_bgr_check(gqr_rule_t *rule, gqr_parameter_t *p,
			gint *imax, gdouble *emax, gboolean analytic,
			FILE *output)

{
  gdouble x0, x1, dx, xb, x, f, *Iq, *Ia, tol, *Qi, *ival ;
  gdouble xi, wi, W ;
  gint i, j, k, nf, nn, ni, nimax, nq ;
  gqr_rule_t *r ;
  gqr_adapt_func_t func ;
  
  if ( gqr_parameter_pointer_number(p) < 1 ) 
    g_error("%s: at least one pointer must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_double_number(p) < 3 ) 
    g_error("%s: at least three doubles must be set in parameters",
	    __FUNCTION__) ;
  if ( gqr_parameter_int_number(p) < 4 ) 
    g_error("%s: at least four ints must be set in parameters",
	    __FUNCTION__) ;

  func = gqr_parameter_pointer(p, 0) ;
  
  x0 = gqr_parameter_double(p, 0) ;
  x1 = gqr_parameter_double(p, 1) ;
  tol = gqr_parameter_double(p, 2) ;

  nf = gqr_parameter_int(p, 0) ;

  Iq = (gdouble *)g_malloc(nf*sizeof(gdouble)) ;
  Ia = (gdouble *)g_malloc(nf*sizeof(gdouble)) ;
  
  gqr_rule_scale(rule, x0, x1, &xb, &dx, &W) ;
  *imax = 0 ; *emax = 0.0 ;
  for ( i = 0 ; i < nf ; i ++ ) {
    Iq[i] = 0.0 ;
    for ( j = 0 ; j < gqr_rule_length(rule) ; j ++ ) {
      x = xb + dx*gqr_rule_abscissa(rule, j) ;
      f = func(x, i, p) ;
      Iq[i] += f*gqr_rule_weight(rule, j)*W ;
    }
  }

  if ( analytic ) {
    x = 0.0 ;
    for ( i = 0 ; i < nf ; i ++ ) {
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

  /*integrate the discretized version of the functions*/
  nq      = gqr_parameter_int(p, 1) ;
  nimax   = gqr_parameter_int(p, 2) ;
  r = gqr_rule_alloc(nq) ;
  gqr_rule_select(r, GQR_GAUSS_LEGENDRE, nq, NULL) ;

  ival = (gdouble *)g_malloc((nimax+1)*sizeof(gdouble)) ;

  ival[0] = x0 ; ival[1] = x1 ;
  Qi = (gdouble *)g_malloc(nq*nq*sizeof(gdouble)) ;
  gqr_discretize_interp(r, Qi) ;

  ni = 1 ;

  /*split interval for discretization*/
  for ( i = 0 ; i < nf ; i ++ ) {
    nn = 1 ;
    while ( nn != 0 && ni < nimax ) {
      nn = gqr_discretize_adaptive(func, i, p, r,
				   Qi, ival, &ni, nimax, tol) ;
    }
  }

  /*integrals of discretized functions*/
  memset(Ia, 0, nf*sizeof(gdouble)) ;
  for ( j = 0 ; j < ni ; j ++ ) {
    for ( k = 0 ; k < nq ; k ++ ) {
      gqr_discretize_quadrature(r, ival, ni, j, k, &xi, &wi) ;
      for ( i = 0 ; i < nf ; i ++ ) {
	Ia[i] += func(xi, i, p)*wi ;
      }
    }
  }

  *emax = 0.0 ; *imax = 0 ;
  for ( i = 0 ; i < nf ; i ++ ) {
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


/**
 * @}
 * 
 */
