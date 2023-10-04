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
#include <math.h>

#include <glib.h>

#include "gqr.h"
#include "gqr-private.h"

gdouble bessel_zero_mcmahon(gdouble nu, gint m)

/*DLMF: https://dlmf.nist.gov/10.21*/
  
{
  gdouble a, mu, j, a8 ;

  a = (m + 0.5*nu - 0.25)*M_PI ; a8 = 8.0*a ;
  mu = 4.0*nu*nu ;

  j = a ;
  j -= (mu-1)/a8 ;               a8 *= 64*a*a ;
  j -=  4*(mu-1)*(7*mu-31)/3/a8 ; a8 *= 64*a*a ;
  j -= 32*(mu-1)*(3779 + mu*(-982 + mu*83))/15/a8 ; a8 *= 64*a*a ;
  j -= 64*(mu-1)*(-6277237 + mu*(1585743 + mu*(-153855 + mu*6949)))/105/a8 ;
  
  return j ;
}
