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
#include <math.h>

#include <glib.h>

#include <config.h>

#include "gqr.h"
#include "gqr-private.h"

/** 
 * Compute a quadrature rule of Gauss-Chebyshev type.
 * 
 * @param n number of points in quadrature rule
 * @param x abscissae of quadrature rule
 * @param w weights of quadrature rule
 * 
 * @return 0 on success
 */
gint grule_chebyshev_1(gint n, gdouble *x, gdouble *w)

{
  gint i ;
  gdouble ww ;
  
  ww = M_PI/(gdouble)n ;
  for ( i = 1 ; i <= n ; i ++ ) {
    x[i-1] = cos((gdouble)(2*i-1)*M_PI/(gdouble)(2*n)) ;
    w[i-1] = ww ;
  }

  return 0 ;
}

gint grule_chebyshev_2(gint n, gdouble *x, gdouble *w)

{
  gint i ;
  
  for ( i = 1 ; i <= n ; i ++ ) {
    x[i-1] = cos((gdouble)i*M_PI/(gdouble)(n+1)) ;
    w[i-1] = (1-x[i-1]*x[i-1])/(gdouble)(n+1)*M_PI;
  }

  return 0 ;
}
