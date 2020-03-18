/* Copyright (C) 2007, 2010 by  Michael Carley */

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

#ifndef GQR_H_INCLUDED
#define GQR_H_INCLUDED

#include <stdio.h>

#include <glib.h>

#if HAVE_LIBMPFR
#ifndef GQR_MPFR_BITS
#define GQR_MPFR_BITS 128
#define GQR_MPFR_ROUND_MODE GMP_RNDN
#endif
#endif /*HAVE_LIBMPFR*/

enum {
  GQR_ERROR = -1,
  GQR_SUCCESS = 0,
  GQR_NULL_PARAMETER,
  GQR_PARAMETER_OUT_OF_RANGE,
  GQR_INVALID_STRING
} ;

/**
 * @typedef gqr_func
 * @ingroup gqr
 *
 * A type for evaluating integrands and their derivatives. 
 * @param t evaluation point;
 * @param i the order of derivative to be evaluated;
 * @param data user data to be passed to the function.
 *
 * @return value of the \f$i\f$th derivative of the function at
 * \f$t\f$.
 */

  typedef gdouble (* gqr_func)(gdouble t, gint i, gpointer data) ;

/**
 * @typedef gqr_t
 * @ingroup gqr
 *
 * A type defining the quadrature rules available. A quadrature rule
 * is defined as a basic type possibly combined with a
 * singularity. For example GQR_GAUSS_LEGENDRE |
 * GQR_GAUSS_HYPERSINGULAR specifies a Gauss-Legendre rule which
 * handles singularities up to and including second order. Note that
 * not all of these rules are implemented yet. GQR will return an
 * error message if you try to use an unimplemented rule.
 * 
 */

typedef enum {
  GQR_GAUSS_LEGENDRE = 1,	/**< Gauss-Legendre*/
  GQR_GAUSS_CHEBYSHEV_1 = 2,	/**< Gauss-Chebyshev of the first kind */
  GQR_GAUSS_CHEBYSHEV_2 = 3,	/**< Gauss-Chebyshev of the second kind */
  GQR_GAUSS_HERMITE = 4,	/**< Hermite */
  GQR_GAUSS_LAGUERRE = 5,	/**< Laguerre */
  GQR_GAUSS_JACOBI = 6,		/**< Jacobi */
  GQR_GAUSS_LOGARITHMIC = 1 << 8, /**< Logarithmic singularities */
  GQR_GAUSS_SINGULAR = 1 << 9,	/**< First order and logarithmic singularities */
  GQR_GAUSS_HYPERSINGULAR = 1 << 10, /**<  Singularities up to second order*/
  GQR_GAUSS_MULTISINGULAR = 1 << 11 /**< Multiple singularities which must be specified in the rule initialization */
} gqr_t ;

#define GQR_GAUSS_REGULAR 0

#define GQR_SINGULARITY_MASK \
  (GQR_GAUSS_LOGARITHMIC | GQR_GAUSS_SINGULAR | GQR_GAUSS_HYPERSINGULAR | \
   GQR_GAUSS_MULTISINGULAR)

#define GQR_RULE_MASK 255

/**
 * @typedef gqr_parameter_t
 * @ingroup gqr
 * 
 * A data structure containing parameters for quadrature rules which
 * require them. 
 */

typedef  struct _gqr_parameter_t gqr_parameter_t ;

struct _gqr_parameter_t {
  gqr_t type ;
  gint i[8], ni, nf ;
  gdouble f[8] ;
} ;

/**
 * @typedef gqr_rule_t
 * @ingroup gqr
 * 
 * Data structure containing the nodes and weights of a quadrature
 * rule. The rule is allocated using ::gqr_rule_alloc and initialized
 * using ::gqr_rule_select. To transform the range of integration use
 * ::gqr_rule_scale.
 */

typedef  struct _gqr_rule_t gqr_rule_t ;

struct _gqr_rule_t {
  gqr_t type ;
  gint n, nmax ;
  gdouble a, b ;
  gdouble *x, *w ;
} ;

#define gqr_rule_length(g) (g->n)
#define gqr_rule_abscissa(g,_i) (g->x[(_i)])
#define gqr_rule_weight(g,_i) (g->w[(_i)])
#define gqr_rule_type(g) (g->type)
#define gqr_parameter_int(p,_i) (p->i[(_i)])
#define gqr_parameter_double(p,_i) (p->f[(_i)])
#define gqr_parameter_clear(p) ((p)->ni=(p)->nf=0)
#define gqr_parameter_clear_int(p) ((p)->ni=0)
#define gqr_parameter_clear_float(p) ((p)->nf=0)
#define gqr_parameter_set_double(p,_f) ((p)->f[((p)->nf)++] = (_f))
#define gqr_parameter_set_int(p,_i) ((p)->i[((p)->ni)++] = (_i))
#define gqr_parameter_ni(p) (((p)->ni))
#define gqr_parameter_nf(p) (((p)->nf))
#define gqr_rule_base_type(t) (t & GQR_RULE_MASK)
#define gqr_rule_singularity_type(t) (t & GQR_SINGULARITY_MASK)

gint gqr_rule_free(gqr_rule_t *g) ;
gqr_rule_t *gqr_rule_alloc(gint n) ;
gqr_rule_t *gqr_rule_realloc(gqr_rule_t *g, gint n) ;
gint gqr_rule_write(gqr_rule_t *g, FILE *f) ;
gint gqr_rule_select(gqr_rule_t *g, gqr_t type, gint n, 
		     gqr_parameter_t *p) ;
gint gqr_rule_scale(gqr_rule_t *g, gdouble a, gdouble b,
		    gdouble *xbar, gdouble *dx) ;
gqr_t gqr_rule_from_name(gchar *str, gqr_parameter_t *p) ;

gdouble gqr_finite_part(gdouble a, gdouble b, gdouble y, gdouble g) ;
gdouble gqr_finite_part_integral(gqr_func f, gpointer data,
				 gdouble y, gdouble gm,
				 gdouble a, gdouble b, gqr_rule_t *g) ;
gdouble gqr_finite_part_1mt_n(gint n, gint m) ;
gdouble gqr_legendre_dPdx(gdouble x, gint n, gint m) ;
gdouble gqr_finite_part_Pn1(gint n, gint m) ;
gdouble gqr_finite_part_Pn_log(gdouble x, gint n) ;
gchar *gqr_rule_name(gqr_t type) ;
gchar *gqr_rule_name_base(gqr_t type) ;
gchar *gqr_rule_name_singularity(gqr_t type) ;

gdouble gqr_finite_part_tn_log(gdouble x, gint n) ;
gint gqr_legendre_integrals(gint N, gdouble *I, gdouble *P) ;

gint gqr_function_nextroot(gdouble *P, gdouble *Q, gdouble *R,
			   gdouble x0, gdouble u0,
			   gdouble *x, gdouble *du) ;
gint gqr_function_roots(gdouble *P, gdouble *Q, gdouble *R,
			gdouble x0, gdouble du0, gint N,
			gdouble *x, gdouble *du) ;
gint gqr_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func) ;

#endif /*GQR_H_INCLUDED*/
