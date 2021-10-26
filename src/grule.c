/* Copyright (C) 2007, 2008, 2010, 2021 by  Michael Carley */

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

/**
 * @defgroup gqr GQR quadrature rules
 * @{
 * 
 * The GQR library computes Gaussian and Gaussian-type quadrature
 * rules including some which handle singular integrands. A simple
 * example of its use to integrate a function is: 
 * <a href="example.c-example.html">example.c</a>.
 *
 * The functions described here are used to generate the quadrature
 * rules and do basic operations such as reading them from and writing
 * them to a file. 
 */

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

gchar *gqr_pointers[] =
  {
   "scattering_radial",
   "basis functions for the radial integrals in Bremer and Gimbutas,\n"
   "\"On the numerical evaluation of the singular integrals of scattering\n"
   "theory\", 2013, https://doi.org/10.1016/j.jcp.2013.05.048\n"
   "(equation 26)\n\n"
   "parameter list:\n\n"
   "  integer: [number of basis functions]\n"
   "           [base quadrature length]\n"
   "           [maximum number of discretization intervals]\n"
   "           [maximum rank for orthogonalization (rule length)]\n"
   "  float:   [start of quadrature interval]\n"
   "           [end of quadrature interval]\n"
   "           [tolerance for discretization]\n"
   "           [ratio of radii (see Bremer and Gimbutas)]\n",
   "scattering_radial_range",
   "basis functions for the radial integrals in Bremer and Gimbutas,\n"
   "\"On the numerical evaluation of the singular integrals of scattering\n"
   "theory\", 2013, https://doi.org/10.1016/j.jcp.2013.05.048\n"
   "(equation 26, with multiple values of d distributed between minimum\n"
   "and maximum)\n\n"
   "parameter list:\n\n"
   "  integer: [number of basis functions]\n"
   "           [base quadrature length]\n"
   "           [maximum number of discretization intervals]\n"
   "           [maximum rank for orthogonalization (rule length)]\n"
   "           [number of radius ratios d]\n"
   "  float:   [start of quadrature interval]\n"
   "           [end of quadrature interval]\n"
   "           [tolerance for discretization]\n"
   "           [minimum ratio of radii]\n"
   "           [maximum ratio of radii]\n",
   "scattering_angular_range",
   "basis functions for the angular integrals in Bremer and Gimbutas,\n"
   "\"On the numerical evaluation of the singular integrals of scattering\n"
   "theory\", 2013, https://doi.org/10.1016/j.jcp.2013.05.048\n"
   "(equation 27, with multiple values of r0 and theta_0)\n\n"
   "parameter list:\n\n"
   "  integer: [number of basis functions]\n"
   "           [base quadrature length]\n"
   "           [maximum number of discretization intervals]\n"
   "           [maximum rank for orthogonalization (rule length)]\n"
   "           [number of basis functions M(u)]\n"
   "           [number of plain trigonometric functions cos, sin]\n"
   "           [number of values of r_0]\n"
   "           [number of values of \theta_0]\n"
   "  float:   [start of quadrature interval]\n"
   "           [end of quadrature interval]\n"
   "           [tolerance for discretization]\n"
   "           [minium radius r_0 (see Bremer and Gimbutas)]\n"
   "           [maximum radius r_0 (see Bremer and Gimbutas)]\n"
   "           [minimum angular range theta_0 (see Bremer and Gimbutas)]\n"
   "           [maximum angular range theta_0 (see Bremer and Gimbutas)]\n",
   "scattering_angular",
   "basis functions for the angular integrals in Bremer and Gimbutas,\n"
   "\"On the numerical evaluation of the singular integrals of scattering\n"
   "theory\", 2013, https://doi.org/10.1016/j.jcp.2013.05.048\n"
   "(equation 27)\n\n"
   "parameter list:\n\n"
   "  integer: [number of basis functions]\n"
   "           [base quadrature length]\n"
   "           [maximum number of discretization intervals]\n"
   "           [maximum rank for orthogonalization (rule length)]\n"
   "           [number of basis functions M(u)]\n"
   "           [number of plain trigonometric functions cos, sin]\n"
   "  float:   [start of quadrature interval]\n"
   "           [end of quadrature interval]\n"
   "           [tolerance for discretization]\n"
   "           [radius r_0 (see Bremer and Gimbutas)]\n"
   "           [angular range theta_0 (see Bremer and Gimbutas)]\n",
   ""			 
} ;

/** 
 * Allocate a Gaussian quadrature rule.
 * 
 * @param n maximum length of rule (number of abscissae and weights).
 * 
 * @return newly allocated rule.
 */

gqr_rule_t *gqr_rule_alloc(gint n)

{
  gqr_rule_t *g ;

  g_return_val_if_fail(n > 0, NULL) ;

  g = (gqr_rule_t *)g_malloc0(sizeof(gqr_rule_t)) ;
  g->nmax = n ; g->n = 0 ;

  g->x = (gdouble *)g_malloc0(sizeof(gdouble)*(n+1)) ;
  g->w = (gdouble *)g_malloc0(sizeof(gdouble)*(n+1)) ;

  return g ;
}

/** 
 * Reallocate memory for a Gaussian quadrature rule to allow for a
 * change in rule length.
 * 
 * @param g ::gqr_rule_t to reallocate;
 * @param n new maximum length of rule (number of abscissae and weights).
 * 
 * @return newly allocated rule.
 */

gqr_rule_t *gqr_rule_realloc(gqr_rule_t *g, gint n)

{
  g_return_val_if_fail(n > 0, NULL) ;

  if ( g == NULL ) return gqr_rule_alloc(n) ;
  if ( g->nmax >= n ) return g ;

  g_free(g->x) ; g_free(g->w) ;
  g->nmax = n ; g->n = 0 ;

  g->x = (gdouble *)g_malloc0(sizeof(gdouble)*(n+1)) ;
  g->w = (gdouble *)g_malloc0(sizeof(gdouble)*(n+1)) ;
  
  /* g_free(g) ; */

  return g ;
  /* (gqr_rule_alloc(n)) ; */
}

/** 
 * Free all memory associated with a Gaussian quadrature rule.
 * 
 * @param g rule to be freed.
 * 
 * @return 0 on success.
 */

gint gqr_rule_free(gqr_rule_t *g)

{
  g_return_val_if_fail(g != NULL, GQR_NULL_PARAMETER) ;

  g_free(g->x) ; g_free(g->w) ; g_free(g) ;
  return 0 ;
}

/** 
 * Write the abscissae and weights of a Gaussian quadrature rule to a
 * file. 
 * 
 * @param g Gaussian quadrature rule;
 * @param f file stream to write to.
 * 
 * @return 0 on success.
 */

gint gqr_rule_write(gqr_rule_t *g, FILE *f)

{
  gint i ;

  g_return_val_if_fail(g != NULL, GQR_NULL_PARAMETER) ;
  g_return_val_if_fail(f != NULL, GQR_NULL_PARAMETER) ;

  for ( i = 0 ; i < gqr_rule_length(g) ; i ++ )
    fprintf(f, "%1.16e %1.16e\n", 
	    gqr_rule_abscissa(g,i),
	    gqr_rule_weight(g,i)) ;
	    

  return 0 ;
}

/** 
 * Fill a Gaussian quadrature rule with abscissae and weights of a
 * chosen type.
 * 
 * @param g rule to fill;
 * @param type a gqr_t for the rule;
 * @param n number of points in rule;
 * @param p gqr_parameter_t of parameters for rule, or NULL for rules
 * with no parameter (e.g. standard Legendre quadrature).
 * 
 * @return 0 on success.
 */

gint gqr_rule_select(gqr_rule_t *g, gqr_t type, gint n, 
		     gqr_parameter_t *p)

{
  g_return_val_if_fail(g != NULL, GQR_NULL_PARAMETER) ;
  
  if ( n > g->nmax ) 
    g_error("%s: number of abscissae (%d) must be less than maximum "
	    "for g (%d)", __FUNCTION__, n, g->nmax) ;
  switch (type) {
  default: 
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	  "quadrature rule %d (%s) not implemented ", 
	  type, gqr_rule_name(type)) ;
    break ;
  case GQR_GAUSS_LEGENDRE:
    grule_legendre(n, g->x, g->w) ; g->n = n ;
    g->a = -1 ; g->b = 1 ; g->type = type ;
    break ;
  case GQR_GAUSS_JACOBI:
    grule_jacobi(n, g->x, g->w, p) ; g->n = n ;
    g->a = -1 ; g->b = 1 ; g->type = type ;
    break ;
  case GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR:
    if ( p->ni <= 0 ) {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	    "%s: integration order M not set for hypersingular "
	    "quadrature rule", 
	    __FUNCTION__) ;
    }
    if ( p->nf <= 0 ) {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	    "%s: singularity position x not set for quadrature rule", 
	    __FUNCTION__) ;
    }
#if HAVE_LAPACK
    if ( p->nf == 1 ) 
      grule_kolm_rokhlin_new(n, 
			      gqr_parameter_int(p,0), 
			      gqr_parameter_double(p,0),
			      g->x, g->w) ;
    else 
      grule_near_singular(n, 
			      gqr_parameter_int(p,0), 
			      gqr_parameter_double(p,0),
			      gqr_parameter_double(p,1),
			      g->x, g->w) ;      
    g->a = -1 ; g->b = 1 ; g->n = n ; g->type = type ;
#else /*HAVE_LAPACK*/
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	  "quadrature rule %d (%s) not implemented (requires LAPACK)", 
	  type, gqr_rule_name(type)) ;
#endif /*HAVE_LAPACK*/
    break ;
  case GQR_GAUSS_LEGENDRE | GQR_GAUSS_MULTISINGULAR:
    g_assert(p->ni > 1) ; g_assert(p->nf > 0) ;
#if HAVE_LAPACK
    if ( p->nf == 1 ) 
      grule_multi_singular(n, 
			   gqr_parameter_int(p,0),
			   gqr_parameter_double(p,0),
			   gqr_parameter_ni(p)-1,
			   &(gqr_parameter_int(p,1)),
			   g->x, g->w) ; 
    else 
      grule_multi_nsingular(n, 
			    gqr_parameter_int(p,0),
			    gqr_parameter_double(p,0),
			    gqr_parameter_double(p,1),
			    gqr_parameter_ni(p)-1,
			    &(gqr_parameter_int(p,1)),
			    g->x, g->w) ;      
    g->a = -1 ; g->b = 1 ; g->n = n ; g->type = type ;
#else /*HAVE_LAPACK*/
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR, 
	  "quadrature rule %d (%s) not implemented (requires LAPACK)", 
	  type, gqr_rule_name(type)) ;
#endif /*HAVE_LAPACK*/
    break ;
  case GQR_GAUSS_LEGENDRE | GQR_GAUSS_LOGARITHMIC:
    if ( (p == NULL) || ((p->ni == 0) && (p->nf == 0)) ) {
      /*if no parameters are specified, default to Smith 
	quadrature rule*/
      grule_logarithmic_smith(n, g->x, g->w) ;
      g->a = 0.0 ; g->b = 1.0 ; g->n = n ; g->type = type ;
    } else {
      g_assert_not_reached() ;
    }
    break ;
  case GQR_GAUSS_CHEBYSHEV_1:
    grule_chebyshev_1(n, g->x, g->w) ; g->n = n ;
    g->a = -1 ; g->b = 1 ; g->type = type ;
    break ;
  case GQR_GAUSS_CHEBYSHEV_2:
    grule_chebyshev_2(n, g->x, g->w) ; g->n = n ;
    g->a = -1 ; g->b = 1 ; g->type = type ;
    break ;
  case GQR_GAUSS_HERMITE:
    grule_hermite(n, g->x, g->w) ; g->n = n ;
    g->a = -1 ; g->b = 1 ; g->type = type ;    
    break ;
  case GQR_GAUSS_GENERALIZED:
    g->n = grule_bgr(g->x, g->w, p) ;
    g_assert(g->n < g->nmax) ;
    g->a = gqr_parameter_double(p, 0) ;
    g->b = gqr_parameter_double(p, 1) ;
    g->type = type ;    
    break ;
  }

  return 0 ;
}

/** 
 * Set the scaling for a rule so that the integral is approximated by:
 * \f$I\approx\Delta x \sum_{j=0}^{n}w_{j}f(x_{j})\f$ where
 * \f$x=\bar{x}+\Delta x t_{j}\f$. 
 * 
 * @param g Gaussian quadrature rule;
 * @param a lower limit of integration;
 * @param b upper limit of integration;
 * @param xbar \f$\bar{x}\f$;
 * @param dx \f$\Delta x\f$.
 * 
 * @return 0 on success.
 */

gint gqr_rule_scale(gqr_rule_t *g, gdouble a, gdouble b,
		    gdouble *xbar, gdouble *dx)

{
  g_return_val_if_fail(g != NULL, GQR_NULL_PARAMETER) ;
  g_return_val_if_fail(xbar != NULL, GQR_NULL_PARAMETER) ;
  g_return_val_if_fail(dx != NULL, GQR_NULL_PARAMETER) ;

  *dx = (b-a)/(g->b-g->a) ; *xbar = a - (g->a)*(*dx) ;

  return 0 ;
}

gchar *gqr_rule_name_singularity(gqr_t type)

{
  gqr_t t ;

  /*strip the basic information*/
  t = type & GQR_SINGULARITY_MASK ;

  if ( t == 0 ) return "non-singular" ;
  if ( t == GQR_GAUSS_LOGARITHMIC ) return "logarithmic" ;
  if ( t == GQR_GAUSS_SINGULAR ) return "Cauchy singular" ;
  if ( t == GQR_GAUSS_HYPERSINGULAR ) return "hypersingular" ;
    
  return "unknown singularity" ;
}

gchar *gqr_rule_name_base(gqr_t type)

{
  gqr_t t ;

  /*strip the singularity information from the leading bits*/
  t = type & GQR_RULE_MASK ;

  if ( t == GQR_GAUSS_LEGENDRE ) 
    return "Gauss Legendre" ;
  if ( t == GQR_GAUSS_CHEBYSHEV_1 )
    return "Gauss Chebyshev of the first kind" ;
  if ( t == GQR_GAUSS_CHEBYSHEV_2 )
    return "Gauss Chebyshev of the second kind" ;
  if ( t == GQR_GAUSS_HERMITE )
    return "Gauss Hermite" ;
  if ( t == GQR_GAUSS_LAGUERRE )
    return "Gauss Laguerre" ;
  if ( t == GQR_GAUSS_JACOBI )
    return "Gauss Jacobi" ;

  return "unknown" ;
}

/** 
 * Generate a string giving the name of the quadrature rule (Legendre,
 * Hermite, etc.) and any singularities built into it.
 * 
 * @param type a ::gqr_t describing a quadrature rule.
 * 
 * @return a string giving the name of rule \a type.
 */

gchar *gqr_rule_name(gqr_t type)

{
  gchar *s ;

  s = g_strdup_printf("%s %s",
		      gqr_rule_name_singularity(type),
		      gqr_rule_name_base(type)) ;

  return s ;
}

static gqr_t string_to_gqr(gchar *s)

{
  if ( !strcmp(s, "GQR_GAUSS_LEGENDRE") ) return GQR_GAUSS_LEGENDRE ;
  if ( !strcmp(s, "GQR_GAUSS_CHEBYSHEV_1") ) return GQR_GAUSS_CHEBYSHEV_1 ;
  if ( !strcmp(s, "GQR_GAUSS_CHEBYSHEV_2") ) return GQR_GAUSS_CHEBYSHEV_2 ;
  if ( !strcmp(s, "GQR_GAUSS_HERMITE") ) return GQR_GAUSS_HERMITE ;
  if ( !strcmp(s, "GQR_GAUSS_LAGUERRE") ) return GQR_GAUSS_LAGUERRE ;
  if ( !strcmp(s, "GQR_GAUSS_JACOBI") ) return GQR_GAUSS_JACOBI ;
  if ( !strcmp(s, "GQR_GAUSS_LOGARITHMIC") ) return GQR_GAUSS_LOGARITHMIC ;
  if ( !strcmp(s, "GQR_GAUSS_SINGULAR") ) return GQR_GAUSS_SINGULAR ;
  if ( !strcmp(s, "GQR_GAUSS_HYPERSINGULAR") ) return GQR_GAUSS_HYPERSINGULAR ;
  if ( !strcmp(s, "GQR_GAUSS_MULTISINGULAR") ) return GQR_GAUSS_MULTISINGULAR ;

  g_error("%s: unrecognized quadrature type %s", __FUNCTION__, s) ;

  return 0 ;
}

static gboolean string_is_int(gchar *str)

{
  gint i ;

  for ( i = 0 ; i < strlen(str) ; i ++ ) {
    if ( !isdigit(str[i]) && str[i] != '+' && str[i] != '-' ) 
      return FALSE ;
  }

  return TRUE ;
}

/** 
 * Parse a string description of a GQR rule, returning a ::gqr_t
 * description, intended for use in parsing command line options.
 * 
 * @param str a string containing a description of a GQR rule,
 * e.g. "GQR_GAUSS_LEGENDRE | GQR_GAUSS_LOGARITHMIC";
 * @param p NULL or a ::gqr_parameter_t which will be filled with 
 * the data relevant to the selected rule type.
 * 
 * @return a ::gqr_t description of the rule in \a str, or an error code. 
 */

gqr_t gqr_rule_from_name(gchar *str, gqr_parameter_t *p)

{
  gqr_t rule, mod ;
  gchar **tokens ;
  gint i ;

  g_return_val_if_fail(str != NULL, GQR_NULL_PARAMETER) ;

  if ( p != NULL ) gqr_parameter_clear(p) ;

/*   g_strdelimit(str, " ,", ' ') ; */
  tokens = g_strsplit(str, " ", 0) ;
  if ( tokens[0] == NULL ) return GQR_INVALID_STRING ;

  rule = 0 ; i = 0 ;
/*   rule = string_to_gqr(tokens[0]) ; */
  while ( tokens[i] != NULL ) {
    if ( string_is_int(tokens[i]) ) {
      if ( p != NULL ) 
	gqr_parameter_set_int(p, atoi(tokens[i])) ;
      else
	g_error("%s: cannot set int values with p == NULL", 
		__FUNCTION__) ;
    } else {
      if ( i == 0 ) rule = string_to_gqr(tokens[i]) ;
      else {
	mod = string_to_gqr(tokens[i+1]) ;
	if ( !strcmp(tokens[i], "|") ) rule = rule | mod ;
	else {
	  if ( !strcmp(tokens[i], "&") ) rule = rule & mod ;
	  else {
	    g_error("%s: cannot apply operator %s in string %s",
		    __FUNCTION__, tokens[i], str) ;
	  }
	}
	i ++ ;
      }
    }
    i ++ ;
  }

  if ( p != NULL ) p->type = rule ;

  return rule ;
}

gpointer gqr_pointer_parse(gchar *str)

{
  if ( !strcmp(str, "scattering_radial") )
    return grule_bgr_func_scattering_r ;
  
  if ( !strcmp(str, "scattering_radial_range") )
    return grule_bgr_func_scattering_range_r ;
  
  if ( !strcmp(str, "scattering_angular") )
    return grule_bgr_func_scattering_th ;

  if ( !strcmp(str, "scattering_angular_range") )
    return grule_bgr_func_scattering_range_th ;

  return NULL ;
}

gint gqr_pointers_list(FILE *output, gboolean verbose)

{
  gint i ;

  i = 0 ; 
  while ( strlen(gqr_pointers[2*i]) != 0 ) {
    fprintf(output, "%s\n", gqr_pointers[2*i+0]) ;
    if ( verbose ) {
      fprintf(output, "\n%s\n", gqr_pointers[2*i+1]) ;
    }
    i ++ ; 
  }
  
  return 0 ;
}

/**
 * @}
 * 
 */
