#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"

gint _gqr_quad_2_r(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log_new(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log(gint N, gdouble x, gdouble y, gdouble *I) ;

gdouble _grule_hpoly(gint n, gdouble x)

{
  gdouble Hnp1, Hn, Hnm1 ;
  gint i ;

  if ( n == 0 ) return sqrt(sqrt(M_1_PI))*exp(-0.5*x*x) ;

  Hnm1 = 0.0 ; Hn = sqrt(sqrt(M_1_PI))*exp(-0.5*x*x) ;
  for ( i = 0 ; i < n ; i ++ ) {
    Hnp1 = sqrt(2.0/(i+1))*x*Hn - sqrt((gdouble)i/(gdouble)(i+1))*Hnm1 ;
    Hnm1 = Hn ; Hn = Hnp1 ;
  }

  return Hn ;
}

gdouble func_exp(gdouble x, gint n, gpointer data)

{
  return exp(x) ;
}

gdouble func_one(gdouble x, gint n, gpointer data)

{
  if ( n == 0 ) return 1.0 ;
  return 0.0 ;
}

gdouble func_Pn(gdouble x, gint q, gpointer data)

{
  gdouble P, t1, t2 ;
  gint n, k ; 

  n = 5 ;

  t1 = 1-x ; t2 = 1+x ;
  for ( (P = 0.0), (k = q) ; k <= n ; k ++ ) {
    P += (gsl_pow_int(-1,q)*gsl_pow_int(t1,k-q) + 
	  gsl_pow_int(-1,n)*gsl_pow_int(t2,k-q))*
      gsl_pow_int(-1,k)/gsl_sf_fact(k)*gsl_sf_fact(n+k)/gsl_sf_fact(n-k)/
      gsl_pow_int(2,k+1)/gsl_sf_fact(k-q) ;
  }

  return P ;
}

gdouble func_1mt(gdouble x, gint n, gpointer data)

{
  gdouble m, pm, y ;

  m = 4.0 ; pm = -1.0 ;
  y = 1.0+pm*x ;
  switch (n) {
  case 0: return pow(y,m) ; break ;
  case 1: return pm*m*pow(y,m-1) ; break ;
  case 2: return m*(m-1)*pow(y,m-2) ; break ;
  default: g_assert_not_reached() ; break ;
  }
  
  return 0.0 ;
}

gdouble func_power(gdouble x, gint m, gint *n)

{
  switch (m) {
  default: g_assert_not_reached() ; break ;
  case 0: return gsl_pow_int(x,(*n)) ; break ;
  case 1: if ( *n == 0) return 0.0 ; 
    return (*n)*gsl_pow_int(x,(*n)-1) ; 
    break ;
  case 2: if ( *n <= 1 ) return 0.0 ; 
    return (*n)*(*n-1)*gsl_pow_int(x,(*n)-2) ; 
    break ;
  }

  return 0.0 ;
}

gdouble func_krsin(gdouble x, gint m, gpointer data)

{
  switch (m) {
  default: g_assert_not_reached() ; break ;
  case 0: return (sin(2*x)+cos(3*x)) ; break ;
  case 1: return (2*cos(2*x)-3*sin(3*x)) ; break ;
  }

  return 0.0 ;  
}

gdouble func_log_test(gdouble t, gint m, gpointer *data)

{
  gdouble x = *(gdouble *)data[0] ;
  gint n = *(gint *)data[1] ;
  gdouble f ;

  f = 0.0 ;
  switch (m) {
  case 0:
    if ( x == t ) f = 0.0 ;
    else f = gsl_pow_int(t, n)*log(fabs(x-t)) ;
    break ;
  case 1:
    f = gsl_pow_int(t, n)/(x-t) ;
    break ;
  case 2:
    f = gsl_pow_int(t, n)/(x-t)/(x-t) ;
    break ;
  default: g_assert_not_reached() ; break ;
  }
  return f ;
}

gdouble log_test(gdouble x, gint q, gint n) 

{
  gint m, k ;
  gdouble I, t, xn, xn1 ;
  static gqr_rule_t *g = NULL ;

  m = n/2 ; I = 0.0 ;

  if ( g == NULL ) {
    g = gqr_rule_alloc(1024) ;
    gqr_rule_select(g, GQR_GAUSS_LEGENDRE, 1024, NULL) ;
  }

  switch (q) {
  default: g_assert_not_reached() ; break ;
  case 0:
    if ( x == -1.0 ) I = (1-gsl_pow_int(-1,n+1))*M_LN2 ;
    if ( x == 1.0 ) I = (1+gsl_pow_int(-1,n))*M_LN2 ;

    if ( (x != -1.0) && (x != 1.0) ) {
      xn = gsl_pow_int(x,n+1) ;
      I = (1-xn)*log(fabs(1-x)) + (gsl_pow_int(-1,n)+xn)*log(fabs(1+x)) ;
    }

    if ( (n) == 2*m ) { /*even n*/
      for ( (k = 0), (t = 2.0) ; k <= m ; k ++ ) {
	I -= t/(gdouble)(2*m-2*k+1) ;
	t *= x*x ;
      }
      I /= (gdouble)(2*m+1) ;
    } else {
      for ( (k = 0), (t = 2.0*x) ; k <= m ; k ++ ) {
	I -= t/(gdouble)(2*m-2*k+1) ;
	t *= x*x ;
      }
      I /= (gdouble)(2*m+2) ;
    }
    break ;
  case 1:
    if ( (x == 1.0) || (x == -1.0) ) {
      return gqr_finite_part_integral(func_power, &n, x, 1.0, -1, 1, g) ;
    }

    if ( (x != -1.0) && (x != 1.0) ) {
      xn = gsl_pow_int(x,n) ; xn1 = gsl_pow_int(x,n+1) ;
      I = -(1-xn1)/(1-x) - (n+1)*xn*log(fabs(1-x))
	+ (gsl_pow_int(-1,n)+xn1)/(1+x) + (n+1)*xn*log(fabs(1+x)) ;
    }

    if ( (n) == 2*m ) { /*even n*/
      for ( (k = 1), (t = 2.0*x) ; k <= m ; k ++ ) {
	I -= 2*k*t/(gdouble)(2*m-2*k+1) ;
	t *= x*x ;
      }
      I /= (gdouble)(2*m+1) ;
    } else {
      for ( (k = 0), (t = 2.0) ; k <= m ; k ++ ) {
	I -= (2*k+1)*t/(gdouble)(2*m-2*k+1) ;
	t *= x*x ;
      }
      I /= (gdouble)(2*m+2) ;
    }
    break ;
  case 2:
    return gqr_finite_part_integral(func_power, &n, x, 2.0, -1, 1, g) ;
    break ;
  }

  return I ;
}

gint main(gint argc, gchar **argv)

{
  gqr_rule_t *g, *gleg ;
  gqr_func func ;
  gqr_parameter_t p ;
  gqr_t type ;
  gpointer fdata[2] ;
  gdouble a, b, y, yy ;
  gdouble y0, y1, dy ;
  gdouble dx, xbar, x ;
  gdouble g0, g1, g2, glog ;
  gint i, j, ns ;
  gint ngp, ngps ;
  gint N, M ;
  gint n, nmax ;
  gchar ch ;
  gdouble *I ;
  gdouble P[3], Q[3], R[3] ;
  gdouble r[512], du[512] ;
  gdouble x0, du0 ;
  gint np ;

  N = 32 ;
  x = 0.4 ; y = 0.3 ;
  np = 64 ;

  N = 3 ;

  /* type = gqr_rule_from_name("GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR",&p) */

/*   P[2] = -1 ; P[1] =  0 ; P[0] = 1 ; */
/*   Q[2] =  0 ; Q[1] = -2 ; Q[0] = 0 ; */
/*   R[2] =  0 ; R[1] =  0 ; R[0] = np*(np+1) ; */

/*   P[2] = 0 ; P[1] = 0 ; P[0] = 1 ; */
/*   Q[2] = 0 ; Q[1] = 0 ; Q[0] = 0 ; */
/*   R[2] = -1 ; R[1] = 0 ; R[0] = 2*np+1 ; */

/*   x0 = 0.0 ; */
/*   du0 = _grule_diff_hpoly(np, x0) ; */

/* /\*   gqr_diff_legendre(np, x0, &du0) ; *\/ */
/*   gqr_function_nextroot(P, Q, R, 0.0, _grule_hpoly(np, 0.0), &x0, &du0) ; */
/*   gqr_function_roots(P, Q, R, x0, du0, N, r, du) ; */

/*   for ( np = 0 ; np < N ; np ++ )  */
/*     fprintf(stdout, "%1.32e %1.32e\n", r[np], du[np]) ; */

  gqr_parameter_clear(&p) ;

  g = gqr_rule_alloc(np) ;
  gqr_rule_select(g, GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR, np, &p) ;
  gqr_rule_write(g, stdout) ;

  return 0 ;

/*   I = (gdouble *)g_malloc((N+1)*sizeof(gdouble)) ; */
/*   _gqr_quad_log(N, x, y, I) ; */
/*   for ( i = 0 ; i <= N ; i ++ ) */
/*     fprintf(stdout, "%1.16e\n", I[i]) ; */
/*   _gqr_quad_log_new(N, x, y, I) ; */
/*   for ( i = 0 ; i <= N ; i ++ ) */
/*     fprintf(stdout, "%1.16e\n", I[i]) ; */
/*   return 0 ; */
  
  /*Chebyshev tests*/
/*   g = gqr_rule_alloc(32) ; */
/*   for ( i = 1 ; i < 5 ; i ++ ) { */
/*     gqr_rule_select(g, GQR_GAUSS_CHEBYSHEV_2, i, NULL) ; */
/*     gqr_rule_write(g, stdout) ; */
/*   } */
/*   a = -1.0 ; b = 1.0 ; */
/*   for ( i = 0 ; i < 16384 ; i ++ )  */
/*     gqr_rule_select(g, GQR_GAUSS_CHEBYSHEV_1, 16384, NULL) ; */
/*   gqr_rule_select(g, GQR_GAUSS_LEGENDRE, 8, NULL) ; */
/*   gqr_rule_write(g, stdout) ; */
/*   return 0 ; */
/*   gqr_rule_scale(g, a, b, &xbar, &dx) ; */
/*   for ( n = 0 ; n < 8 ; n ++ ) { */
/*     for ( (i = 0), (g0 = 0.0) ; i < gqr_rule_length(g) ; i ++ ) { */
/*       x = gqr_rule_abscissa(g,i)*dx + xbar ; */
/*       g0 += gqr_rule_weight(g,i)*dx*gsl_pow_int(x,n) ; */
/*     } */
/*     fprintf(stdout, "%d %1.16e\n", n, g0) ; */
/*   } */
/*   return 0 ; */

  /*logarithmic Green's function test*/
  func = func_log_test ;
  fdata[0] = &y ; fdata[1] = &n ;
  y0 = -2.0 ; y1 = 1.0 ; dy = 0.0125 ;
/*   dy = 4.0/257.0 ;  */
  ns = (gint)floor((y1 - y0)/dy)+1 ;
  ngp = 1024 ; ngps = 2048 ;
  a = -1.0 ; b = 1.0 ;
  g = gqr_rule_alloc(ngps) ; 
  gleg = gqr_rule_alloc(ngp) ;

  M = 4 ; N = 16 ; nmax = 4 ; y = 0.0 ;
  while ( (ch = getopt(argc, argv, "M:N:n:y:")) != EOF ) {
    switch(ch) {
    default: g_assert_not_reached() ; break ;
    case 'M': M = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'n': nmax = atoi(optarg) ; break ;
    case 'y': y = atof(optarg) ; break ;
    }
  }
  
  gqr_parameter_clear(&p) ;
  gqr_parameter_set_int(&p, M) ;
  gqr_parameter_set_double(&p, 0.4) ;
  gqr_parameter_set_double(&p, 0.3) ;
  gqr_rule_select(g, GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR, N, &p) ;
  gqr_rule_write(g, stdout) ;


  gqr_parameter_clear(&p) ;
  gqr_parameter_set_int(&p, M) ;
  gqr_parameter_set_int(&p, 0) ;
/*   gqr_parameter_set_int(&p, 1) ; */
  gqr_parameter_set_int(&p, 2) ;
  gqr_parameter_set_double(&p, 0.4) ;
  gqr_parameter_set_double(&p, 0.3) ;
  gqr_rule_select(g, GQR_GAUSS_LEGENDRE | GQR_GAUSS_MULTISINGULAR, N, &p) ;
  gqr_rule_write(g, stdout) ;
  

  return 0 ;

  for ( j = 0 ; j < ns ; j ++ ) {
    yy = y0 + j*dy ;
    gqr_parameter_clear(&p) ;
    gqr_parameter_set_int(&p, M) ;
    if ( gsl_fcmp(yy, 1.0, 2e-16) == 0) yy = 1.0 ;
    gqr_parameter_set_double(&p, yy) ;
    gqr_rule_select(g, GQR_GAUSS_LEGENDRE | GQR_GAUSS_HYPERSINGULAR, N, &p) ;
    gqr_rule_select(gleg, GQR_GAUSS_LEGENDRE, ngp, NULL) ;

    gqr_rule_scale(g, a, b, &xbar, &dx) ;
    y = xbar + dx*yy ;
    for ( n = 0 ; n < nmax ; n ++ ) {
      for ( ( i = 0 ), ( glog = g0 = g1 = g2 = 0.0 ) ; 
	    i < gqr_rule_length(g) ; i ++ ) {
	x = gqr_rule_abscissa(g,i)*dx + xbar ;
	glog += (func_log_test(x, 2, fdata)+
		 func_log_test(x, 1, fdata)+
		 func_log_test(x, 0, fdata))
	  *gqr_rule_weight(g,i)*dx ;
      }
      fprintf(stdout, "%d %g %g %g %g\n",
	      n, y, glog, 
	      log_test(y, 2, n)+log_test(y, 1, n)+log_test(y, 0, n),
	      log_test(y, 2, n)+log_test(y, 1, n)+log_test(y, 0, n)-glog) ;
    }
  }

  return 0 ;
}
