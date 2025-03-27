/* Copyright (C) 2025 by  Michael Carley */

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
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include <ctype.h>

#include <glib.h>

#include "gqr.h"
#include "gqr-private.h"

/*
 * Paget, D.F. A quadrature rule for finite-part integrals. BIT 21,
 * 212–220 (1981). https://doi.org/10.1007/BF01933166
 */

gint grule_paget(gint n, gdouble *x, gdouble *w, gqr_parameter_t *p)

{
  gint a, k, i ;
  gdouble b[1024] = {0}, c[1024] = {0}, d[1024] = {0}, E2, E1, E ;
  
  if ( p == NULL || gqr_parameter_int_number(p) == 0 ) {
    a = 2 ;
  } else {
    a = gqr_parameter_int(p,0) ;
  }

  grule_legendre(n, x, w) ;

  if ( a > 1 ) b[0] = 1.0/(1-a) ;

  for ( k = 1 ; k <= a-2 ; k ++ ) {
    b[k] = -(k+a-1)/(k-a+1)*b[k-1] ;
  }

  c[a-1] = (2*a-2)*b[a-2] ;

  for ( k = 1 ; k <= a-1 ; k ++ ) d[a-1] -= 1.0/k/(k+a-1) ;
  d[a-1] *= a-1 ;
  b[a-1] = c[a-1]*d[a-1] ;
  
  for ( k = a ; k < n+2 ; k ++ ) {
    c[k] = -(gdouble)(k+a-1)/(k-a+1)*c[k-1] ;
    d[k] = d[k-1] + 1.0/(k+a-1) + 1.0/(k-a+1) ;
    b[k] = c[k]*d[k] ;
  }

  /* for ( i = 0 ; i <= n ; i ++ ) { */
  /*   fprintf(stderr, "%e %e %e\n", b[i], c[i], d[i]) ; */
  /* } */
  
  for ( i = 0 ; i < n ; i ++ ) {
    E2 = E1 = 0 ;
    for ( k = n-1 ; k >= 0 ; k -- ) {
      E = (k+0.5)*(b[k] + 2.0*x[i]/(k+1)*E1) - (gdouble)(k+1)/(k+2)*E2 ;
      E2 = E1 ; E1 = E ;      
      /* fprintf(stderr, "%lg\n", E) ; */
    }
    w[i] *= E ;
    x[i] = 0.5*(1.0+x[i]) ;
  }
  
  return 0 ;
}
