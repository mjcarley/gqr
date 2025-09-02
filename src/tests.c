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
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include <glib.h>

#include "gqr.h"
#include "gqr-private.h"

/**
 * @defgroup tests Testing GQR quadrature rules
 * @{
 * 
 * GQR has a facility for checking most of its standard rules against
 * analytical integration
 * 
 */

static gdouble beta_func(gdouble x, gdouble y)

/*
 * Gradshteyn and Ryzhik 8.384.1 
 */
  
{
  gdouble B ;

  B = tgamma(x)*tgamma(y)/tgamma(x+y) ;
  
  return B ;
}

static void quad_monomial_weighted(gdouble a, gdouble b,
				   gdouble al, gdouble bt,
				   gdouble *I, gint n)

/*
 *
 * integrals of x^n*(b-x)^\alpha*(x-a)^\beta
 * 
 * Gradshteyn and Ryzhik 3.196.3 and a recursion
 *
 */
  
{
  gdouble mu, nu ;
  gint j ;
  
  mu = al + 1.0 ; nu = bt + 1.0 ;
  I[0] = beta_func(mu, nu)*pow(b-a, mu+nu-1) ;
  I[1] = (mu*b + nu*a)*I[0]/(mu + nu) ;

  for ( j = 1 ; j < n ; j ++ ) {
    I[j+1] = ((a+b)*j + b*mu + a*nu)*I[j] - a*b*j*I[j-1] ;
    I[j+1] /= mu + nu + j ;
  }
  
  return ;
}

static gdouble quad_eval(gqr_rule_t *g, gdouble a, gdouble b,
			 gqr_parameter_t *p,
			 gqr_test_integrand_t func,
			 gint j)

{
  gdouble tbar, dt, I, t, W ;
  gint i ;
  
  gqr_rule_scale(g, a, b, &tbar, &dt, &W) ;

  I = 0.0 ;
  for ( i = 0 ; i < gqr_rule_length(g) ; i ++ ) {
    t = gqr_rule_abscissa(g,i)*dt + tbar ;
    I += W*gqr_rule_weight(g,i)*func(t, p, j) ;
  }
  
  return I ;
}

gdouble gqr_test_integrand_monomial(gdouble t, gqr_parameter_t *p, gint j)

{
  gdouble f ;

  f = pow(t, j) ;
  
  return f ;
}


/** 
 * Test function for Gauss-Jacobi quadratures, evaluating
 * \f$\int_{a}^{b}x^i*(b-x)^\alpha*(x-a)^\beta\,\mathrm{d}x\f$,
 * \f$i=0,\ldots,n\f$. Singularity strengths \f$\alpha\f$ and
 * \f$\beta\f$ are the first two double entries in the parameter list \a p. 
 * 
 * @param a lower limit of integration (singularity location);
 * @param b upper limit of integration (singularity location);
 * @param p parameter list;
 * @param I on exit contains integrals of weighted monomials;
 * @param n number of integrals to evaluate.
 * 
 * @return 0 on success.
 */

gint gqr_test_integral_jacobi(gdouble a, gdouble b, gqr_parameter_t *p,
			      gdouble *Q, gint n)

{
  gdouble al, bt ;

  al = gqr_parameter_double(p, 0) ;
  bt = gqr_parameter_double(p, 1) ;

  if ( al <= -1 )
    g_error("%s: alpha (%lg) out of range for Gauss-Jacobi quadrature",
	    __FUNCTION__, al) ;
  if ( bt <= -1 )
    g_error("%s: beta (%lg) out of range for Gauss-Jacobi quadrature",
	    __FUNCTION__, al) ;

  quad_monomial_weighted(a, b, al, bt, Q, n) ;
  
  return 0 ;
}

gint gqr_test_integral_hermite(gdouble a, gdouble b, gqr_parameter_t *p,
			       gdouble *Q, gint n)

{
  gint i ;

  for ( i = 0 ; i <= n ; i += 2 ) {
    Q[i+0] = tgamma(0.5*(i+1)) ;
    Q[i+1] = 0 ;
  }
  
  return 0 ;
}

gint gqr_test_integral_laguerre(gdouble a, gdouble b, gqr_parameter_t *p,
			       gdouble *Q, gint n)

{
  gint i ;

  for ( i = 0 ; i < n ; i += 2 ) {
    Q[i+0] = tgamma(0.5*(i+1)) ;
    Q[i+1] = 0 ;
  }
  
  return 0 ;
}

gint gqr_test_integral_chebyshev_1(gdouble a, gdouble b, gqr_parameter_t *p,
				   gdouble *Q, gint n)

{
  gdouble al, bt ;

  al = -0.5 ; bt = -0.5 ;

  quad_monomial_weighted(a, b, al, bt, Q, n) ;
  
  return 0 ;
}

gint gqr_test_integral_chebyshev_2(gdouble a, gdouble b, gqr_parameter_t *p,
				   gdouble *Q, gint n)

{
  gdouble al, bt ;

  al = 0.5 ; bt = 0.5 ;

  quad_monomial_weighted(a, b, al, bt, Q, n) ;
  
  return 0 ;
}

/** 
 * Test function for Gauss-Legendre quadratures,
 * \f$\int_{a}^{b}x^{i}\,\mathrm{d}x\f$, \f$i=0,\ldots,n\f$.
 * 
 * @param a lower limit of integration (singularity location);
 * @param b upper limit of integration (singularity location);
 * @param p parameter list;
 * @param I on exit contains integrals of monomials;
 * @param n number of integrals to evaluate.
 * 
 * @return 0 on success.
 */

gint gqr_test_integral_legendre(gdouble a, gdouble b, gqr_parameter_t *p,
				gdouble *Q, gint n)

{
  gdouble Ib, Ia ;
  gint i ;

  Ib = b ; Ia = a ;
  for ( i = 0 ; i < n ; i ++ ) {
    Q[i] = (Ib - Ia) ;
    Ia *= a*(i+1)/(i+2) ;
    Ib *= b*(i+1)/(i+2) ;
  }
  
  return 0 ;
}

/** 
 * Check a GQR quadrature rule against analytic evaluation
 * 
 * @param g rule to be tested;
 * @param a lower limit of integration;
 * @param b upper limit of integration;
 * @param ifunc function to be used as integrand in numerical integration;
 * @param intfunc value of integral on interval, calculated analytically or
 * otherwise;
 * @param n number of integrals to evaluate;
 * @param tol tolerance for error check
 * @param emax maximum difference between numerical and analytical evaluation;
 * @param Imax maximum absolute value of integral;
 * @param imax index of integral with largest absolute error.
 * 
 * @return TRUE if *emax < tol*Imax, FALSE otherwise.
 */

gboolean gqr_test_rule_base(gqr_rule_t *g, gdouble a, gdouble b,
			    gqr_test_integrand_t ifunc,
			    gqr_test_integral_t intfunc,
			    gint n, gdouble tol,
			    gdouble *emax, gdouble *Imax, gint *imax)

{
  gdouble Ic, Ia[512] ;
  gint i ;
  gqr_parameter_t *p = gqr_rule_data(g) ;

  *imax = 0 ; *emax = 0.0 ; *Imax = 0.0 ;

  intfunc(a, b, p, Ia, n) ;

  for ( i = 0 ; i < n ; i ++ ) {
    Ic = quad_eval(g, a, b, p, ifunc, i) ;
    if ( fabs(Ic - Ia[i]) > (*emax) ) {
      *imax = i ; *emax = fabs(Ic - Ia[i]) ;
    }
    *Imax = MAX(Ia[i], (*Imax)) ;
  }

  if ( *emax > tol*(*Imax) ) return FALSE ;
  
  return TRUE ;
}

/** 
 * Test a GQR quadrature rule by calling ::gqr_test_rule_base with
 * appropriate integrand and integral functions, based on the type of
 * quadrature in the supplied rule
 * 
 * @param g rule to test;
 * @param a lower limit of integration;
 * @param b upper limit of integration;
 * @param tol tolerance for error check;
 * @param emax maximum difference between exact and computed integrals;
 * @param Imax maximum absolute value of an integral;
 * @param imax index of integral with largest error.
 * 
 * @return TRUE if emax < tol, FALSE otherwise.
 */

gboolean gqr_test_rule(gqr_rule_t *g, gdouble a, gdouble b,
		       gdouble tol, gdouble *emax, gdouble *Imax,
		       gint *imax)

{
  gqr_t baserule ;
  gint n ;
  gqr_test_integrand_t ifunc ;
  gqr_test_integral_t intfunc ;
  
  baserule = gqr_rule_type(g) & GQR_RULE_MASK ;

  ifunc = NULL ;
  intfunc = NULL ;
  
  switch ( baserule ) {
  default:
    g_error("%s: unrecognized rule %u", __FUNCTION__, baserule) ;
    break ;
  case GQR_GAUSS_LEGENDRE:
    if ( baserule == gqr_rule_type(g) ) 
      n = gqr_rule_length(g)*2 ;
    else
      n = gqr_parameter_int(&(g->data),0) ;
    ifunc = gqr_test_integrand_monomial ;
    intfunc = gqr_test_integral_legendre ;
    break ;
  case GQR_GAUSS_CHEBYSHEV_1:
    n = gqr_rule_length(g)*2 ;
    ifunc = gqr_test_integrand_monomial ;
    intfunc = gqr_test_integral_chebyshev_1 ;
    break ;
  case GQR_GAUSS_CHEBYSHEV_2:
    n = gqr_rule_length(g)*2 ;
    ifunc = gqr_test_integrand_monomial ;
    intfunc = gqr_test_integral_chebyshev_2 ;
    break ;
  case GQR_GAUSS_HERMITE:
    n = gqr_rule_length(g)*2 ;
    ifunc = gqr_test_integrand_monomial ;
    intfunc = gqr_test_integral_hermite ;
    break ;
  case GQR_GAUSS_LAGUERRE: break ;
  case GQR_GAUSS_JACOBI:
    n = gqr_rule_length(g)*2 ;
    ifunc = gqr_test_integrand_monomial ;
    intfunc = gqr_test_integral_jacobi ;
    break ;	
  }

  if ( ifunc == NULL )
    g_error("%s: tests for %s rule not yet implemented",
	    __FUNCTION__, gqr_rule_name_base(baserule)) ;
  
  return gqr_test_rule_base(g, a, b, ifunc, intfunc, n, tol, emax, Imax, imax) ;
}

/**
 * @}
 * 
 */
