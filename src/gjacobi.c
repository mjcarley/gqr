/* Copyright (C) 2021 by Michael Carley */

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

#include <stdio.h>
#include <math.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"
#include "gqr-private.h"

static inline void grule_diff_jpoly(gint n, gdouble x,
				    gdouble al, gdouble bt,
				    gdouble *f, gdouble *df)

{
  gdouble Pnp1, Pn, Pnm1 ;
  gdouble dPnp1, dPn, dPnm1 ;
  gdouble A, B, C ;
  gint i ;

  if ( n == 0 ) {
    *f = *df = 0.0 ;
    return ;
  }

  i = 0 ;
  A = 0.5*(2.0*i+al+bt+1.0)*(2.0*i+al+bt+2.0)/(1.0+i)/(1.0+i+al+bt) ;
  /* B = 0.5*(al - bt)*(2.0*i+al+bt+1.0)/(1.0+i)/(1.0+al+bt+i) ; */
  B = 0.5*(bt - al)*(2.0*i+al+bt+1.0)/(1.0+i)/(1.0+al+bt+i) ;
  if ( n == 1 ) {
    *f = A*x + B ; *df = A ;
    return ;
  }
  Pnm1 = 1.0 ; Pn = A*x + B ;  
  dPnm1 = 0.0 ; dPn = A ;
    
  for ( i = 1 ; i < n ; i ++ ) {
    A = 0.5*(2.0*i+al+bt+1.0)*(2.0*i+al+bt+2.0)/(1.0+i)/(1.0+i+al+bt) ;
    /* B = 0.5*(al*al - bt*bt)*(2.0*i+al+bt+1.0)/(1.0+i)/(1.0+al+bt+i)/ */
    /*   (al+bt+2.0*i) ; */
    B = 0.5*(bt*bt - al*al)*(2.0*i+al+bt+1.0)/(1.0+i)/(1.0+al+bt+i)/
      (al+bt+2.0*i) ;
    C = (al+i)*(bt+i)*(2.0*i+al+bt+2.0)/(1.0+i)/(1.0+al+bt+i)/(2.0*i+al+bt) ;

    Pnp1 = (A*x + B)*Pn - C*Pnm1 ;
    dPnp1 = A*Pn + (A*x + B)*dPn - C*dPnm1 ;

    Pnm1 = Pn ; Pn = Pnp1 ;
    dPnm1 = dPn ; dPn = dPnp1 ;
  }

  *f = Pn ; *df = dPn ;
  
  return ;
}

static void jacobi_root(gint n, gdouble al, gdouble bt,
			gdouble *x, gdouble *df)

{
  gdouble tol, f ;
  gint i ;
  
  tol = 1e-14 ;
  grule_diff_jpoly(n, *x, al, bt, &f, df) ;

  for ( i = 0 ; i < 16 ; i ++ ) {
    *x -= f/(*df) ;
    grule_diff_jpoly(n, *x, al, bt, &f, df) ;
    if ( fabs(f) < tol ) return ;
  }

  g_assert_not_reached() ;
  
  return ;
}

/*from Hale and Townsend, 2013*/

static gdouble root_interior(gint n, gdouble al, gdouble bt, gint k)

{
  gdouble rho, phi, x ;

  rho = 0.5*(al+bt+1.0) + n ;
  phi = M_PI*(k + 0.5*al - 0.25)/rho ;
  
  x = (0.25 - al*al)/tan(0.5*phi) - (0.25 - bt*bt)*tan(0.5*phi) ;
  x = cos(phi + 0.25*x/rho/rho) ;
  
  return x ;
}

gint grule_jacobi(gint n, gdouble *x, gdouble *w, gqr_parameter_t *pm)

/*
 * coefficients of differential equation from DLMF
 * https://dlmf.nist.gov/18.8 
 * p = 1-x^2 ; 
 * q = \beta-\alpha -(\alpha+\beta+2)x; 
 * r = n(n+\alpha+\beta+1)
 */
    
{
  gdouble al, bt, A ;
  /* du, x0 ; */
  /* gdouble p[] = {1, 0, -1}, q[] = {0, -2, 0}, r[] = {0, 0, 0} ; */
  gint i ;
  
  if ( gqr_parameter_nf(pm) < 2 )
    g_error("%s: Jacobi quadrature rule requires two parameters",
	    __FUNCTION__) ;
  al = gqr_parameter_double(pm, 0) ;
  bt = gqr_parameter_double(pm, 1) ;

  /* fprintf(stderr, "%lg %lg\n", al, bt) ; */
  
  /* q[0] = bt - al ; */
  /* q[1] = -(al+bt+2.0) ; */

  /* r[0] = (al+bt+n+1.0)*(gdouble)n ; */

  /* x0 = 0 ; */
  /* jacobi_root(n, al, bt, &x0, &du) ; */
  /* gqr_function_roots(p, q, r, x0, du, n/2, x, w) ; */

  /*
   * bit crude but works: use asymptotic estimate of root (Hale and
   * Townsend) followed by Newton polishing
   *
   * definition of weights from
   * https://mathworld.wolfram.com/Jacobi-GaussQuadrature.html
   */
  A = lgamma(n+al+1) + lgamma(n+bt+1) - lgamma(n+al+bt+1) - lgamma(n+1)
    + (al+bt+1)*M_LN2 ;
  A = exp(A) ;
  for ( i = 0 ; i < n ; i ++ ) {
    x[i] = root_interior(n, al, bt, n-i) ;
    jacobi_root(n, al, bt, &(x[i]), &(w[i])) ;

    w[i] = A/(1-x[i]*x[i])/w[i]/w[i] ;    
  }
  
  return 0 ;
}
