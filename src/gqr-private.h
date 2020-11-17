/* Copyright (C) 2008 by  Michael Carley */

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

#ifndef GQR_PRIVATE_H_INCLUDED
#define GQR_PRIVATE_H_INCLUDED

#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#if HAVE_LAPACK
#ifndef GQR_LAPACK_COND_TOL
#define GQR_LAPACK_COND_TOL 1e-15
#endif
#endif /*HAVE_LAPACK*/

#ifndef g_debug
#define g_debug(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, format)
#endif
#ifndef g_message
#define g_message(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, format)
#endif
#ifndef g_warning
#define g_warning(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_WARNING, format)
#endif
#ifndef g_error
#define g_error(format...) g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, format)
#endif 
#define GQR_LOGGING_DATA_WIDTH     4
#define GQR_LOGGING_DATA_FID       0
#define GQR_LOGGING_DATA_PREFIX    1
#define GQR_LOGGING_DATA_LEVEL     2
#define GQR_LOGGING_DATA_EXIT_FUNC 3

gdouble log2(gdouble x) ;
gint grule_korsunsky(gint n, gint k, gdouble *xk, gdouble *w) ;
gint gqr_lsqr_min_norm(gdouble *A, gint m, gint n, gdouble *b, gint nb) ;
gint gqr_legendre_singular(gint N, gint x, gint m, gdouble *I) ;
gint gqr_quad_log(gint N, gdouble x, gdouble y, gdouble *I) ;
gint gqr_quad_1_r(gint N, gdouble x, gdouble y, gdouble *I) ;
gint gqr_quad_2_r(gint N, gdouble x, gdouble y, gdouble *I) ;
gint gqr_quad_3_r(gint N, gdouble x, gdouble y, gdouble *I) ;

gint grule_kolm_rokhlin_new(gint n, gint m, gdouble y,
			     gdouble *x, gdouble *w) ;
gint grule_near_singular(gint n, gint m, gdouble x, gdouble y,
			  gdouble *xk, gdouble *wt) ;
gint grule_legendre(gint n, gdouble *x, gdouble *w) ;
gint grule_chebyshev_1(gint n, gdouble *x, gdouble *w) ;
gint grule_chebyshev_2(gint n, gdouble *x, gdouble *w) ;
gint grule_logarithmic_smith(gint n, gdouble *x, gdouble *w) ;
gint grule_hermite(gint n, gdouble *x, gdouble *w) ;

gint grule_multi_singular(gint n, gint m, gdouble x,
			  gint ns, gint *s, 
			  gdouble *xk, gdouble *wt) ;
gint grule_multi_nsingular(gint n, gint m, gdouble x, gdouble y,
			   gint ns, gint *s, 
			   gdouble *xk, gdouble *wt) ;
gint grule_bgr(gdouble *x, gdouble *w, gqr_parameter_t *p) ;

gint rrqr(gdouble *A, gint m, gint n, gdouble *tau, gint *jpvt,
	  gdouble *work, gint lwork) ;
gint rrqr_rank(gdouble *R, gint m, gint n, gdouble ee) ;
gint rrqr_qr(gdouble *A, gint m, gint n, gdouble *tau, gint rank,
	     gdouble *Q, gdouble *R11, gdouble *work, gint lwork) ;
gint gqr_srrqr(gdouble *A, gint m, gint n, gdouble f, gdouble tol,
	       gdouble *Q, gdouble *R, gint *pvt, gint *rank, gint *ldr,
	       gdouble *work, gint lwork) ;

gdouble grule_bgr_func_scattering_r(gdouble t, gint i, gqr_parameter_t *p) ;

#endif /*GQR_PRIVATE_H_INCLUDED*/
