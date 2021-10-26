/* Copyright (C) 2007, 2013 by  Michael Carley */

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

/** 
 * Compute a quadrature rule for the integration of functions with
 * singularities up to second order using the algorithm of Carley, M.,
 * SIAM Journal on Scientific Computing, 29(3):1207-1216, 2007,
 * http://dx.doi.org/10.1137/060666093, see also, Kolm, P. and
 * Rokhlin, V., `Numerical quadratures for singular and hypersingular
 * integrals', Computers and Mathematics with Applications,
 * 41:327--352, 2001.
 * 
 * @param n number of points in quadrature rule
 * @param m maximum order of polynomial to integrate
 * @param y position of singularity
 * @param x abscissae of quadrature rule
 * @param w weights of quadrature rule
 * 
 * @return 0 on success
 */

gint grule_kolm_rokhlin_new(gint n, gint m, gdouble y,
			    gdouble *x, gdouble *w)
			 
{
  gdouble *fp2 ;
  gdouble *P, *Q, *Py, *dPy ;
  gdouble Pn, Pnm1, Pnp1 ;
  gdouble fp, pv ;
  gdouble t ;
  gqr_rule_t *g, *h ;
  gint i, j, k ;
  gdouble *psi, *M ;

  g = gqr_rule_alloc(n) ;
  gqr_rule_select(g, GQR_GAUSS_LEGENDRE, n, NULL) ;
  h = gqr_rule_alloc(m) ;
  gqr_rule_select(h, GQR_GAUSS_LEGENDRE, m, NULL) ;
  k = gqr_rule_length(g) ;
  Py = (gdouble *)g_malloc((m+2)*sizeof(gdouble)) ;
  dPy = (gdouble *)g_malloc((m+1)*sizeof(gdouble)) ;
  P = (gdouble *)g_malloc((n+1)*k*sizeof(gdouble)) ;
  Q = (gdouble *)g_malloc((n+2)*sizeof(gdouble)) ;
  M = (gdouble *)g_malloc(MAX(4*m,n)*sizeof(gdouble)) ;
  psi = (gdouble *)g_malloc(4*m*n*sizeof(gdouble)) ;
  fp2 = (gdouble *)g_malloc(4*n*sizeof(gdouble)) ;

  for ( i = 0 ; i < 4*m*n ; i ++ ) psi[i] = 0.0 ;
  i = 0 ; 
  for ( j = 0 ; j < k ; j ++ ) {
    P[i*k + j] = 1.0 ;
    P[(i+1)*k + j] = gqr_rule_abscissa(g,j) ;
  }
  if ( ( y != -1.0 ) && ( y != 1.0 ) ) {
    Q[0] = 0.5*log(fabs((y+1)/(y-1))) ;
    Q[1] = y*Q[0] - 1.0 ;
  } else {
    Q[0] = M_LN2 ; Q[1] = M_LN2 - 0.5 ;
  }

  for ( i = 1 ; i < n ; i ++ ) {
    for ( j = 0 ; j < k ; j ++ ) {
      P[(i+1)*k + j] = 
	((gdouble)(2*i+1)*gqr_rule_abscissa(g,j)*P[i*k+j] - 
	 (gdouble)i*P[(i-1)*k+j])/(gdouble)(i+1) ;
    }
    Q[i+1] = ((gdouble)(2*i+1)*y*Q[i] - (gdouble)i*Q[i-1])/(gdouble)(i+1) ;
  }
  i = n ;
  Q[i+1] = ((gdouble)(2*i+1)*y*Q[i] - (gdouble)i*Q[i-1])/(gdouble)(i+1) ;

  Py[0] = 1.0 ; Py[1] = y ;
  for ( i = 1 ; i <= m ; i ++ ) {
    Py[i+1] = ((gdouble)(2*i+1)*y*Py[i] - (gdouble)i*Py[i-1])/(gdouble)(i+1) ;
  }

  dPy[0] = 0.0 ; dPy[1] = 1.0 ;
  if ( (y == 1.0) || (y == -1.0) ) 
    for ( i = 1 ; i < m ; i ++ )
      dPy[i+1] = y*gsl_pow_int(y,i+1)*0.5*(gdouble)((i+1)*(i+2)) ;
  else
    for ( i = 1 ; i < m ; i ++ )
      dPy[i+1] = (i+1)*(y*Py[i+1] - Py[i])/(y*y-1.0) ; 

  if ( (y == 1.0) || (y == -1.0) ) fp2[0] = 2.0*M_LN2 - 2.0 ;
  else fp2[0] = 2*Q[1] + log(fabs(y*y-1.0)) ;

  if ( (y == 1.0) || (y == -1.0) ) 
    for ( i = 1 ; i < n ; i ++ ) fp2[i] = gqr_finite_part_Pn_log(y,i) ;
  else
    for ( i = 1 ; i < n ; i ++ ) fp2[i] = 2.0/(2.0*i+1)*(Q[i+1]-Q[i-1]) ;

  pv = gqr_finite_part(-1, 1, y, 1.0) ;
  fp = gqr_finite_part(-1, 1, y, 2.0) ;
  for ( i = 0 ; i < m ; i ++ ) {
    M[i] = 0.0 ; M[1*m+i] = fp2[i] ;
    M[2*m+i] = Py[i]*pv ;
    M[3*m+i] = Py[i]*fp - dPy[i]*pv ;
  }

  for ( i = 0 ; i < m ; i ++ ) {
    t = gqr_rule_abscissa(h, i) ;
    Pnm1 = 1.0 ; Pn = t ;
    for ( j = 1 ; j < m ; j ++ ) {
      M[2*m+j] += gqr_rule_weight(h,i)*(Pn-Py[j])/(y-t) ;
      M[3*m+j] += 
	gqr_rule_weight(h,i)*(Pn-Py[j]+dPy[j]*(y-t))/(y-t)/(y-t) ;
      Pnp1 = ((gdouble)(2*j+1)*gqr_rule_abscissa(h, i)*Pn -
	      (gdouble)j*Pnm1)/(gdouble)(j+1) ;
      Pnm1 = Pn ; Pn = Pnp1 ;
    }
  }
  M[0*m+0] = 2.0 ;

  for ( i = 0 ; i < n ; i ++ ) {
    t = gqr_rule_abscissa(g,i) ;
    for ( j = 0 ; j < m ; j ++ ) {
      psi[i*4*m+j+0*m] = P[j*k+i] ;
      psi[i*4*m+j+1*m] = psi[i*4*m+j+0*m]*log(fabs(y-t)) ;
      psi[i*4*m+j+2*m] = psi[i*4*m+j+0*m]/(y-t) ;
      psi[i*4*m+j+3*m] = psi[i*4*m+j+2*m]/(y-t) ;
    }
  }
  gqr_lsqr_min_norm(psi, 4*m, n, M, MAX(4*m,n)) ;

  for ( i = 0 ; i < n ; i ++ ) {
    x[i] = gqr_rule_abscissa(g,i) ; w[i] = M[i] ;
  }

  gqr_rule_free(g) ; gqr_rule_free(h) ;
  g_free(P) ; g_free(Q) ; g_free(Py) ; g_free(dPy) ;
  g_free(M) ; g_free(fp2) ; g_free(psi) ;

  return 0 ;
}
