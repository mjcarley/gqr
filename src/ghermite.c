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



#include <stdio.h>
#include <glib.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "gqr.h"
#include "gqr-private.h"

#define SQRT_SQRT_1_PI 0.75112554446494250726118480088189

static inline gdouble grule_hpoly(gint n, gdouble x)

{
  gdouble Hnp1, Hn, Hnm1 ;
  gint i ;

  if ( n == 0 ) return SQRT_SQRT_1_PI*exp(-0.5*x*x) ;

  Hnm1 = 0.0 ; Hn = SQRT_SQRT_1_PI*exp(-0.5*x*x) ;
  for ( i = 0 ; i < n ; i ++ ) {
    Hnp1 = sqrt(2.0/(i+1))*x*Hn - sqrt((gdouble)i/(gdouble)(i+1))*Hnm1 ;
    Hnm1 = Hn ; Hn = Hnp1 ;
  }

  return Hn ;
}

static inline gdouble grule_diff_hpoly(gint n, gdouble x)

{
  gdouble Hnp1, Hn, Hnm1 ;
  gdouble dHnp1, dHn, dHnm1 ;
  gint i ;

  if ( n == 0 ) return 0.0 ;

  Hnm1 = 0.0 ; Hn = SQRT_SQRT_1_PI*exp(-0.5*x*x) ;
  dHnm1 = 0.0 ; dHn = 0.0 ;
  for ( i = 0 ; i < n ; i ++ ) {
    Hnp1 = sqrt(2.0/(gdouble)(i+1))*x*Hn - 
      sqrt((gdouble)i/(gdouble)(i+1))*Hnm1 ;
    dHnp1 = sqrt(2.0/(gdouble)(i+1))*(x*dHn + Hn) -
      sqrt((gdouble)i/(gdouble)(i+1))*dHnm1 ;
    Hnm1 = Hn ; Hn = Hnp1 ;
    dHnm1 = dHn ; dHn = dHnp1 ;
  }

  return dHn ;
}

gint grule_hermite(gint n, gdouble *x, gdouble *w)

{
  gdouble p[] = {1, 0, 0}, q[] = {0, 0, 0}, r[] = {0, 0, -1} ;
  gdouble x0, du ;
  gint i ;

  r[0] = 2*n + 1 ; 
  if ( 2*(n/2) != n ) {
    du = grule_diff_hpoly(n, 0.0) ;
    gqr_function_roots(p, q, r, 0, du, n/2+1, &(x[n/2]), &(w[n/2])) ;
    for ( i = n ; i > n/2 ; i -- ) {
      w[n-i-1] = w[i] = 2.0*exp(-x[i]*x[i])/(w[i]*w[i]) ;
      x[n-i-1] = -x[i] ;
    }
    w[n/2] = 2.0/(w[n/2]*w[n/2]) ;
    return 0 ;
  }

  gqr_function_nextroot(p, q, r, 0.0, grule_hpoly(n, 0.0), &x0, &du) ;
  gqr_function_roots(p, q, r, x0, du, n/2, &(x[n/2]), &(w[n/2])) ;

  for ( i = n-1 ; i >= n/2 ; i -- ) {
    w[n-i-1] = w[i] = 2.0*exp(-x[i]*x[i])/(w[i]*w[i]) ;
    x[n-i-1] = -x[i] ; 
  }

  return 0 ;
}


