/* Copyright (C) 2007 by  Michael Carley */

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

#include <math.h>

#include <glib.h>

#include "gqr.h"
#include "gqr-private.h"

/*
  Functions to implement the method of Korsunsky, A. M.,
  `Gauss-Chebyshev formulae for strongly singular integrals',
  Quarterly of Applied Mathematics, LVI(3):461--472, 1998.
*/

gint grule_korsunsky(gint n, gint k, gdouble *xk, gdouble *w)

{
  gint i ;
  gdouble ti, sgn ;

  g_assert(k <= n+1) ; g_assert(k >= 1) ;

  *xk = cos(0.5*(gdouble)(2*k-1)*M_PI/(gdouble)(n+1)) ;
  sgn = -1.0 ;
  for ( (sgn = -1.0), (i = 1) ; i <= k ; (i ++), (sgn = -sgn)) ;
  for ( i = 1 ; i <= n ; (i ++), (sgn = -sgn)) {
    ti = cos((gdouble)i*M_PI/(gdouble)(n+1)) ;
    w[i-1] = sgn*(1-ti*ti)/(ti-(*xk))/sqrt(1-(*xk)*(*xk))*M_PI ;
  }  

  return 0 ;
}
