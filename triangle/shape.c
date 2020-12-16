/* This file is part of GQR, a library for Gaussian quadratures
 *
 * Copyright (C) 2020 Michael Carley
 *
 * GQR is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  GQR is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GQR.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <glib.h>

#define vector_cross(GQR_C,GQR_A,GQR_B)					\
  ((GQR_C)[0] = (GQR_A)[1]*(GQR_B)[2] - (GQR_A)[2]*(GQR_B)[1],		\
   (GQR_C)[1] = (GQR_A)[2]*(GQR_B)[0] - (GQR_A)[0]*(GQR_B)[2],		\
   (GQR_C)[2] = (GQR_A)[0]*(GQR_B)[1] - (GQR_A)[1]*(GQR_B)[0])

#define vector_scalar(GQR_A,GQR_B)				\
  (((GQR_A)[0])*((GQR_B)[0])+					\
   ((GQR_A)[1])*((GQR_B)[1])+					\
   ((GQR_A)[2])*((GQR_B)[2]))

#define vector_length(GQR_A)					\
  (sqrt(((GQR_A)[0])*((GQR_A)[0])+				\
	((GQR_A)[1])*((GQR_A)[1]) +				\
	((GQR_A)[2])*((GQR_A)[2])))


gint element_shape_3d(gint ne, gdouble s, gdouble t,
		      gdouble *L, gdouble *dLds, gdouble *dLdt)

{
  switch ( ne ) {
  default: g_assert_not_reached() ; break ;
  case 3:
    L[0] = 1.0 - s - t ; 
    L[1] =       s     ; 
    L[2] =           t ;
    if ( dLds == NULL ) break ;
    dLds[0] = -1.0 ; dLds[1] =  1.0 ; dLds[2] =  0.0 ;
    dLdt[0] = -1.0 ; dLdt[1] =  0.0 ; dLdt[2] =  1.0 ;
    break ;
  case 6:
    L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
    L[1] = 2.0*s*(s-0.5) ;
    L[2] = 2.0*t*(t-0.5) ;
    L[3] = 4.0*s*(1.0-s-t) ;
    L[4] = 4.0*s*t ;
    L[5] = 4.0*t*(1.0-s-t) ;

    if ( dLds == NULL ) break ;

    dLds[0] = -3.0 + 4.0*s + 4.0*t ;
    dLds[1] = -1.0 + 4.0*s ;
    dLds[2] =  0.0 ;
    dLds[3] =  4.0 - 8.0*s - 4.0*t ;
    dLds[4] =  4.0*t ;
    dLds[5] = -4.0*t ;

    dLdt[0] = -3.0 + 4.0*s + 4.0*t ;
    dLdt[1] =  0.0 ;
    dLdt[2] = -1.0 + 4.0*t ;
    dLdt[3] = -4.0*s ;
    dLdt[4] =  4.0*s ;
    dLdt[5] =  4.0 - 4.0*s - 8.0*t ;
    break ;
  }

  return 0 ;
}

gint element_point_interp_3d(gdouble *xe, gint xstr, gint ne,
			     gdouble *L, gdouble *dLds, gdouble *dLdt,
			     gdouble *y, gdouble *n, gdouble *J)

{
  gint i ;
  gdouble dyds[3]={0.0}, dydt[3]={0.0} ;
  
  y[0] = y[1] = y[2] = n[0] = n[1] = n[2] = 0.0 ;
  for ( i = 0 ; i < ne ; i ++ ) {
    y[0] += xe[i*xstr+0]*L[i] ; y[1] += xe[i*xstr+1]*L[i] ;
    y[2] += xe[i*xstr+2]*L[i] ; 
    dyds[0] += xe[i*xstr+0]*dLds[i] ; dyds[1] += xe[i*xstr+1]*dLds[i] ;
    dyds[2] += xe[i*xstr+2]*dLds[i] ; 
    dydt[0] += xe[i*xstr+0]*dLdt[i] ; dydt[1] += xe[i*xstr+1]*dLdt[i] ;
    dydt[2] += xe[i*xstr+2]*dLdt[i] ; 
  }

  vector_cross(n, dyds, dydt) ;
  
  *J = vector_length(n) ;

  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ;

  return 0 ;
}

gint element_point_3d(gdouble *xe, gint xstr, gint ne,
		      gdouble s, gdouble t,
		      gdouble *y, gdouble *n, gdouble *J)

{
  gdouble L[32]={0.0}, dLds[32]={0.0}, dLdt[32]={0.0} ;
  
  element_shape_3d(ne, s, t, L, dLds, dLdt) ;
  element_point_interp_3d(xe, xstr, ne, L, dLds, dLdt, y, n, J) ;

  return 0 ;
}
