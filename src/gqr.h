/* Copyright (C) 2007, 2010, 2020, 2023 by  Michael Carley */

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

#ifdef DOXYGEN

/**
 * @file   gqr.h
 * @author  <michael@michael.paraffinalia.co.uk>
 * @date   Wed Mar  6 16:41:09 2024
 * 
 * @brief  
 * 
 * 
 */
#endif /*DOXYGEN*/

enum {
  GQR_ERROR = -1,
  GQR_SUCCESS = 0,
  GQR_NULL_PARAMETER,
  GQR_PARAMETER_OUT_OF_RANGE,
  GQR_INVALID_STRING
} ;

#ifdef DOXYGEN
/**
 * @typedef gqr_t
 * 
 * A type defining the quadrature rules available. A quadrature rule
 * is defined as a basic type possibly combined with a
 * singularity. For example ::GQR_GAUSS_LEGENDRE |
 * ::GQR_GAUSS_HYPERSINGULAR specifies a Gauss-Legendre rule which
 * handles singularities up to and including second order. Note that
 * not all of these rules are implemented or going to be
 * implemented. GQR will return an error message if you try to use an
 * unimplemented rule.
 */
#endif /*DOXYGEN*/
typedef enum
  {
    GQR_GAUSS_LEGENDRE = 1,	/**< Gauss-Legendre*/
    GQR_GAUSS_CHEBYSHEV_1 = 2,	/**< Gauss-Chebyshev of the first kind */
    GQR_GAUSS_CHEBYSHEV_2 = 3,	/**< Gauss-Chebyshev of the second kind */
    GQR_GAUSS_HERMITE = 4,	/**< Hermite */
    GQR_GAUSS_LAGUERRE = 5,	/**< Laguerre */
    GQR_GAUSS_JACOBI = 6,	/**< Jacobi with weight function \f$(1+x)^\alpha (1-x)^\beta\f$, with \f$(\alpha,\beta)\f$ supplied as the first two floating point parameters */
    GQR_GAUSS_LOGARITHMIC = 1 << 8, /**< logarithmic singularities */
    GQR_GAUSS_SINGULAR = 1 << 9,	/**< first order and logarithmic singularities */
    GQR_GAUSS_HYPERSINGULAR = 1 << 10, /**< singularities up to second order*/
    GQR_GAUSS_MULTISINGULAR = 1 << 11, /**< multiple singularities which must be specified in the rule initialization */
    GQR_GAUSS_GENERALIZED = 1 << 12 /**< generalized Gaussian quadratures using the algorithm of Bremer, Gimbutas and Rokhlin */
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

#define GQR_PARAMETER_POINTER_NUMBER_MAX 8
#define GQR_PARAMETER_INTEGER_NUMBER_MAX 8
#define GQR_PARAMETER_FLOAT_NUMBER_MAX   8

typedef  struct _gqr_parameter_t gqr_parameter_t ;

struct _gqr_parameter_t {
  gqr_t type ;
  gint i[GQR_PARAMETER_INTEGER_NUMBER_MAX], ni, nf, np ;
  gdouble f[GQR_PARAMETER_FLOAT_NUMBER_MAX] ;
  gpointer p[GQR_PARAMETER_POINTER_NUMBER_MAX], rules[2] ;
} ;

#define gqr_parameter_int(_p,_i) ((_p)->i[(_i)])
#define gqr_parameter_double(_p,_i) ((_p)->f[(_i)])
#define gqr_parameter_pointer(_p,_i) ((_p)->p[(_i)])
#define gqr_parameter_clear(_p)			\
  do {						\
  (_p)->ni=(_p)->nf=(_p)->np=0 ;		\
  (_p)->rules[0] = (_p)->rules[1] = NULL ;	\
  } while ( 0 )
#define gqr_parameter_clear_int(_p) ((_p)->ni=0)
#define gqr_parameter_clear_double(_p) ((_p)->nf=0)
#define gqr_parameter_set_double(_p,_f) ((_p)->f[((_p)->nf)++] = (_f))
#define gqr_parameter_set_pointer(_p,_f) ((_p)->p[((_p)->np)++] = (_f))
#define gqr_parameter_set_int(_p,_i) ((_p)->i[((_p)->ni)++] = (_i))
#define gqr_parameter_ni(_p) (((_p)->ni))
#define gqr_parameter_nf(_p) (((_p)->nf))
#define gqr_parameter_np(_p) (((_p)->np))

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
  gdouble a, b, *x, *w ;
  gqr_parameter_t data ;
} ;

#define gqr_rule_length(_g)           ((_g)->n)
#define gqr_rule_abscissa(_g,_i)      ((_g)->x[(_i)])
#define gqr_rule_weight(_g,_i)        ((_g)->w[(_i)])
#define gqr_rule_type(_g)             ((_g)->type)
#define gqr_rule_base_type(_t)        ((_t) & GQR_RULE_MASK)
#define gqr_rule_singularity_type(_t) ((_t) & GQR_SINGULARITY_MASK)
#define gqr_rule_data(_g)             (&((_g)->data))

/**
 * @typedef gqr_func_t
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

  typedef gdouble (* gqr_func_t)(gdouble t, gint i, gpointer data) ;

/**
 * @typedef gqr_adapt_func_t
 * @ingroup gqr
 *
 * A type for evaluating functions to be used in adaptive
 * discretizations. Optionally, if \a i is negative, the function call
 * should be interpreted as a request for the integral of the function
 * between limits defined in the first two floating point entries in
 * the parameter list \a p, with the function index taken as -(i+1).
 * 
 * @param t evaluation point;
 * @param i the index of the function to be evaluated;
 * @param p parameter list to be passed to the function.
 *
 * @return value of the function \f$f_i(t)\f$, or the integral of
 * \f$f_{-i-1}(t)\f$ if \f$i<0\f$.
 */

typedef gdouble (* gqr_adapt_func_t)(gdouble t, gint i, gqr_parameter_t *p) ;

typedef gdouble (* gqr_test_integrand_t)(gdouble t, gqr_parameter_t *p,
					 gint j) ;
typedef gint (* gqr_test_integral_t)(gdouble a, gdouble b, gqr_parameter_t *p,
				     gdouble *I, gint n) ;

gint gqr_rule_free(gqr_rule_t *g) ;
gqr_rule_t *gqr_rule_alloc(gint n) ;
gqr_rule_t *gqr_rule_realloc(gqr_rule_t *g, gint n) ;
gint gqr_rule_write(gqr_rule_t *g, FILE *f) ;
gint gqr_rule_select(gqr_rule_t *g, gqr_t type, gint n, 
		     gqr_parameter_t *p) ;
gint gqr_rule_scale(gqr_rule_t *g, gdouble a, gdouble b,
		    gdouble *xbar, gdouble *dx, gdouble *W) ;
gqr_t gqr_rule_from_name(gchar *str, gqr_parameter_t *p) ;

gdouble gqr_finite_part(gdouble a, gdouble b, gdouble y, gdouble g) ;
gdouble gqr_finite_part_integral(gqr_func_t f, gpointer data,
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
gint gqr_legendre_integrals(gint N, gdouble *i, gdouble *P) ;

gint gqr_function_nextroot(gdouble *P, gdouble *Q, gdouble *R,
			   gdouble x0, gdouble u0,
			   gdouble *x, gdouble *du) ;
gint gqr_function_roots(gdouble *P, gdouble *Q, gdouble *R,
			gdouble x0, gdouble du0, gint N,
			gdouble *x, gdouble *du) ;
gint gqr_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func) ;
gpointer gqr_pointer_parse(gchar *str) ;

gint gqr_discretize_interp(gqr_rule_t *rule, gdouble *Q) ;
gint gqr_discretize_adaptive(gqr_adapt_func_t func, gint idx, gpointer data,
			     gqr_rule_t *rule, gdouble *Q,
			     gdouble *ival, gint *ni, gint nimax, 
			     gdouble tol) ;
gint gqr_discretize_eval(gqr_rule_t *rule, gdouble *Q,
			 gdouble *ival, gint ni,
			 gdouble *u, gint nu, gint ustr,
			 gdouble x, gdouble *f) ;
gint gqr_discretize_fill(gqr_adapt_func_t func, gint idx, gpointer data,
			 gqr_rule_t *rule, gdouble *ival, gint ni,
			 gdouble *u, gint ustr) ;
gint gqr_discretize_quadrature(gqr_rule_t *rule,
			       gdouble *ival, gint ni, gint i, gint j,
			       gdouble *x, gdouble *w) ;
gint gqr_discretize_orthogonalize(gdouble *A, gint m, gint n,
				  gdouble tol,
				  gint *rank, gint rankmax,
				  gdouble *Q, gdouble *R11,
				  gint *pvt, gint *ldr,
				  gdouble *work, gint lwork) ;
gint gqr_rule_bgr_check(gqr_rule_t *rule, gqr_parameter_t *p,
			gint *imax, gdouble *emax, gboolean analytic,
			FILE *output) ;
gint gqr_pointers_list(FILE *output, gboolean verbose) ;

gboolean gqr_test_rule_base(gqr_rule_t *g, gdouble a, gdouble b,
			    gqr_test_integrand_t ifunc,
			    gqr_test_integral_t intfunc,
			    gint nfunc, gdouble tol,
			    gdouble *emax, gdouble *Imax, gint *imax) ;
gint gqr_test_integral_legendre(gdouble a, gdouble b, gqr_parameter_t *p,
				gdouble *I, gint n) ;
gint gqr_test_integral_jacobi(gdouble a, gdouble b, gqr_parameter_t *p,
			      gdouble *I, gint n) ;
gint gqr_test_integral_chebyshev_1(gdouble a, gdouble b, gqr_parameter_t *p,
				   gdouble *I, gint n) ;
gint gqr_test_integral_chebyshev_2(gdouble a, gdouble b, gqr_parameter_t *p,
				   gdouble *I, gint n) ;
gint gqr_test_integral_hermite(gdouble a, gdouble b, gqr_parameter_t *p,
			       gdouble *I, gint n) ;
gint gqr_test_integral_laguerre(gdouble a, gdouble b, gqr_parameter_t *p,
				gdouble *I, gint n) ;
gdouble gqr_test_integrand_monomial(gdouble t, gqr_parameter_t *p, gint j) ;
gboolean gqr_test_rule(gqr_rule_t *g, gdouble a, gdouble b,
		       gdouble tol, gdouble *emax, gdouble *Imax, 
		       gint *imax) ;
#endif /*GQR_H_INCLUDED*/
