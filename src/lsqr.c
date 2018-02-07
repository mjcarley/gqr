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



#include <stdio.h>
#include <math.h>

#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include "gqr.h"
#include "gqr-private.h"

extern void dgelsd_(int *, int *, int *,
		    double *, int *, double *, int *, double *,
		    double *, int *, double *, int *, int *,
		    int *) ;

#define log2(_x) (log((_x))/M_LN2)

/* gdouble log2(gdouble x)  */

/* { */
/*   return (log(x)/M_LN2) ; */
/* } */

#ifdef HAVE_LAPACK
gint gqr_lsqr_min_norm(gdouble *A, gint m, gint n, gdouble *b, gint nb)

{
  gdouble *s, rcond, *w, tmp ;
  gint rank, nw, *iw, i, tti ;
  gint nlvl, smlsiz = 64 ;
  gint minmn ;
  gint lda ;
  gint c_1 = 1 ;

  s = (gdouble *)g_malloc(2*MAX(m,n)*sizeof(gdouble)) ;
  lda = m ;
  nlvl = MAX(0, (gint)ceil(log2((gdouble)MIN(m,n)/(gdouble)(smlsiz+1)))) ;

  nw = -1 ; tti = 0 ;
  dgelsd_((int *)&m, (int *)&n, (int *)&c_1,
  	  A, (int *)&lda, b, (int *)&nb, s,
  	  &rcond, (int *)&rank, &tmp, (int *)&nw, (int *)&tti,
  	  (int *)&i) ;
  nw = tmp ;
  minmn = MIN(n,m) ;
  w = (gdouble *)g_malloc(nw*sizeof(gdouble)) ;

  i = 3*minmn*nlvl + 11*minmn ;
  iw = (gint *)g_malloc(2*i*sizeof(glong)) ; 
  rcond = GQR_LAPACK_COND_TOL ;
  dgelsd_((int *)&m, (int *)&n, (int *)&c_1,
  	  A, (int *)&lda, b, (int *)&nb, s,
  	  &rcond, (int *)&rank, w, (int *)&nw, (int *)iw,
  	  (int *)&i) ;

  g_free(s) ; g_free(w) ; g_free(iw) ;

  return 0 ;
}
#else /*HAVE_LAPACK*/
gint gqr_lsqr_min_norm(gdouble *A, gint m, gint n, gdouble *b, gint nb)

{
  g_assert_not_reached() ;

  return 0 ;
}
#endif /*HAVE_LAPACK*/
