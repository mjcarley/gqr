/* Copyright (C) 2023 by  Michael Carley */

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
#include <glib.h>
#include <math.h>

#include "gqr.h"
#include "gqr-private.h"

static inline gdouble grule_Lpoly(gint n, gdouble al, gdouble x)

{
  gdouble Lnp1, Ln, Lnm1, A, B, C ;
  gint i ;

  if ( x == 0 ) return tgamma(al+1+n)/tgamma(al+1)/tgamma(n+1) ;
  
  Lnm1 = 1.0 ; i = 0 ;
  Lnm1 = exp(-0.5*x) ;
  if ( n == 0 ) return Lnm1 ;

  A = -1.0/(i+1) ; B = (2.0*i+al+1)/(i+1) ; C = (al + i)/(i+1) ;
  Ln = A*x + B ;
  Ln *= exp(-0.5*x) ;
  if ( n == 1 ) return Ln ;
  for ( i = 1 ; i < n ; i ++ ) {
    A = -1.0/(i+1) ; B = (2.0*i+al+1)/(i+1) ; C = (al + i)/(i+1) ;
    Lnp1 = (A*x + B)*Ln - C*Lnm1 ;
    Lnm1 = Ln ; Ln = Lnp1 ;
  }

  return Ln ;
}

static inline gdouble grule_diff_Lpoly(gint n, gdouble al, gdouble x)

{
  if ( n == 0 ) return 0.0 ;
  
  return -grule_Lpoly(n-1, al+1, x) ;
}

gint grule_laguerre(gint n, gdouble *x, gdouble *w, gqr_parameter_t *pm)

{
  gdouble p[] = {0, 1, 0}, q[] = {0, 0, 0}, r[] = {0, 0, 0} ;
  gdouble x0, du ;
  gdouble f0, f1, fmid, x1, xmid, al, ee, nu, f, df, jroot ;
  gint i, j, m ;

  if ( gqr_parameter_nf(pm) < 1 ) al = 0.0 ;
  else al = gqr_parameter_double(pm, 0) ;

  /* p[0] = 0 ; p[1] = 0 ; p[2] = 1 ; */
  /* q[0] = 0 ; q[1] = al + 1 ; q[2] = 0 ; */
  /* r[0] = 0 ; r[1] = 0.5*(2*n + al + 1) ; r[2] = -0.25 ; */
  p[2] = 0 ; p[0] = 0 ; p[1] = 1 ;
  q[2] = 0 ; q[0] = al + 1 ; q[1] = 0 ;
  r[2] = 0 ; r[0] = 0.5*(2*n + al + 1) ; r[1] = -0.25 ;

  ee = 0 ; f0 = grule_Lpoly(n, al, ee) ;
  /* ee = 0 ;  f0 = grule_Lpoly(n, al, ee) ; */
  gqr_function_nextroot(p, q, r, ee, f0, &x0, &du) ;
  gqr_function_roots(p, q, r, x0, du, n, x, w) ;

  for ( i = 0 ; i < n ; i ++ ) {
    f = grule_Lpoly(n, al, x[i]) ;
    w[i] = f ;
  }
  
  return 0 ;

  /*from DLMF https://dlmf.nist.gov/18.16#iv*/
  nu = 4.0*n + 2.0*al + 2 ;
  for ( i = 0 ; i < n ; i ++ ) {
    m = i + 1 ;
    x1 = (4.0*m+2*al+2)*(2.0*m + al + 1 +
			 sqrt((2*m+al+1)*(2*m+al+1) + 0.25 - al*al))/nu ;
    f1 = grule_Lpoly(n, al, x1) ; 
    jroot = bessel_zero_mcmahon(al, m) ;
    x0 = jroot*jroot/nu ;
    f0 = grule_Lpoly(n, al, x0) ; 

    /* if ( (f1 > 0 && f0 > 0) || (f1 < 0 && f1 < 0) ) { */
    /*   x[i] = 0.5*(x0+x1) ;  */
    /* } else { */
    {
      for ( j = 0 ; (j < 4) ; j ++ ) {
	xmid = 0.5*(x0 + x1) ;
	fmid = grule_Lpoly(n, al, xmid) ;
	if ( (fmid < 0 && f0 < 0) || (fmid > 0 && f0 > 0) ) {
	  x0 = xmid ; f0 = fmid ;
	} else {
	  x1 = xmid ; f1 = fmid ;
	}
      }
      x[i] = xmid ;
    }
    /* x[i] = 0.5*(x0 + x1) ; */

    for ( j = 0 ; j < 8 ; j ++ ) {
      f  = grule_Lpoly(n, al, x[i]) ;
      df = grule_diff_Lpoly(n, al, x[i]) ;
      x[i] = (df*x[i] - f)/df ;
    }
    /* f = grule_Lpoly(n+1, al, x[i]) ; */
    /* w[i] = tgamma(n+al+1)*x[i]/tgamma(n+1)/(n+1)/(n+1)/f/f ; */
    f = grule_Lpoly(n, al, x[i]) ;
    w[i] = f ;
  }

  return 0 ;
}


