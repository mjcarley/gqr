/**
 * @example example.c
 *
 * An example program which does an integration.
 * 
 */

#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include <gqr.h>

gdouble func(gdouble t, gdouble x, gdouble y)

{
  return (1.0/((x-t)*(x-t) + y*y)) ;
}

gint main()

{
  gqr_rule_t *g ;
  gqr_parameter_t p ;
  gdouble x, y ;
  gdouble a, b ;
  gdouble dx, xbar, W ;
  gdouble I, t ;
  gint i, N, M ;
  gqr_t rule ;

  M = 4 ; N = 16 ; 
  x = 0.3 ; y = 0.2 ;
  a = -1 ; b = 1 ;
  rule = GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR ;
  g = gqr_rule_alloc(N) ;

  gqr_parameter_clear(&p) ;
  gqr_parameter_set_int(&p, M) ;
  gqr_parameter_set_double(&p, x) ;
  gqr_parameter_set_double(&p, y) ;

  gqr_rule_select(g, rule, N, &p) ;

  gqr_rule_scale(g, a, b, &xbar, &dx, &W) ;

  I = 0.0 ;
  for ( i = 0 ; i < gqr_rule_length(g) ; i ++ ) {
    t = gqr_rule_abscissa(g,i)*dx + xbar ;
    I += W*gqr_rule_weight(g,i)*func(t, x, y) ;
  }

  fprintf(stdout, "I: %lg\n", I) ;

  return 0 ;
}
