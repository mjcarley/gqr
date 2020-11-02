#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"

gint _gqr_quad_2_r(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log_new(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log(gint N, gdouble x, gdouble y, gdouble *I) ;

gint dgeqpx_(gint *job, gint *m, gint *n, gint *k,
	     gdouble *A, gint *lda, gdouble *C, gint *ldc,
	     gint *jpvt, gdouble *ircond, gdouble *orcond,
	     gint *rank, gdouble *svlues, gdouble *work,
	     gint *lwork, gint *info) ;

gchar *tests[] = {"discretization",
		  "orthogonalization",
		  ""} ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}


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

gint legendre_interp_matrix_test(gint n, gdouble x0, gdouble x1)

{
  gdouble Q[8192] ;
  gqr_rule_t *rule ;
  gint i, j ;
  
  rule = gqr_rule_alloc(n) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, n, NULL) ;

  gqr_discretize_interp(rule, Q) ;

  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ )
      fprintf(stdout, "%lg ", Q[i*n+j]) ;
    fprintf(stdout, "\n") ;
  }
  
  return 0 ;
}

gdouble discretization_test_func(gdouble t, gint i, gpointer data)

{
  gdouble f ;
  
  f = sin(t) ;
  
  return f ;
}

gint legendre_discretization_test(gint nq, gdouble x0, gdouble x1, gdouble tol)

{
  gdouble ival[8193], Q[8192], emax, *u, x, f[32], g[32] ;
  gqr_adapt_func_t func = discretization_test_func ;
  gqr_rule_t *rule ;
  gint ni, nimax, nn, nx, ustr, nu, i, idx ;
  gpointer data = NULL ;
  
  nimax = 1024 ; ustr = 3 ; nu = 1 ; idx = 0 ; nx = 1024 ;

  ival[0] = x0 ; ival[1] = x1 ;
  
  rule = gqr_rule_alloc(nq) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, nq, NULL) ;

  gqr_discretize_interp(rule, Q) ;

  ni = 1 ;

  nn = 1 ;
  while ( nn != 0 ) {
    nn = gqr_discretize_adaptive(func, idx, data, rule,
				 Q, ival, &ni, nimax, tol) ;
  }

  fprintf(stderr, "ni = %d\n", ni) ;

  /*fill function array*/
  u = (gdouble *)g_malloc(nq*(ni+1)*ustr*sizeof(gdouble)) ;
  gqr_discretize_fill(func, 0, data, rule, ival, ni, u, ustr) ;

  /*check evaluation*/
  emax = 0.0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    x = x0 + (x1 - x0)*i/nx ;
    gqr_discretize_eval(rule, Q, ival, ni, u, nu, ustr, x, g) ;
    f[0] = func(x, idx, data) ;
    emax = MAX(emax, fabs(f[0] - g[0])) ;
    fprintf(stdout, "%lg %lg %lg\n", x, f[0], g[0]) ;
  }

  fprintf(stderr, "max error: %lg\n", emax) ;
  
  return 0 ;
}

gdouble orthogonalization_test_func(gdouble t, gint i, gpointer data)

{
  if ( t == 0 ) return 1.0/t ;

  return pow(t, i-1) ;
}

gint orthogonalization_test(gint nq, gdouble x0, gdouble x1, gdouble tol)

{
  gdouble ival[8193], Q[8192], emax, *u, x, f[32], g[32], *A, *C, *work ;
  gqr_adapt_func_t func = orthogonalization_test_func ;
  gqr_rule_t *rule ;
  gint ni, nimax, nn, nx, ustr, nu, i, j, idx, nf, n, m, k, lwork ;
  gint jpvt[1024], rank, info, job, nfunc ;
  gdouble ircond, orcond, svlues[1024] ;
  gpointer data = NULL ;

  nf = 8 ; ustr = nf ;
  nimax = 1024 ; nu = 1 ; idx = 0 ; nx = 1024 ;

  ival[0] = x0 ; ival[1] = x1 ;
  
  rule = gqr_rule_alloc(nq) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, nq, NULL) ;

  gqr_discretize_interp(rule, Q) ;

  ni = 1 ;

  /*split interval for discretization*/
  for ( i = 0 ; i < nf ; i ++ ) {
    nn = 1 ;
    while ( nn != 0 ) {
      nn = gqr_discretize_adaptive(func, idx, data, rule,
				   Q, ival, &ni, nimax, tol) ;
    }
  }

  fprintf(stderr, "ni = %d\n", ni) ;

  /*number of points in discretized functions*/
  nfunc = nq*ni ;
  
  /*use FORTRAN indexing from here*/
  /*fill function array*/
  u = (gdouble *)g_malloc(nfunc*nf*sizeof(gdouble)) ;
  A = (gdouble *)g_malloc(nfunc*nf*sizeof(gdouble)) ;
  /*this needs checking*/
  C = (gdouble *)g_malloc0(nfunc*nf*sizeof(gdouble)) ;
  /*fill columns of matrix A with functions*/
  ustr = 1 ;
  for ( i = 0 ; i < nf ; i ++ ) {
    gqr_discretize_fill(func, i, data, rule, ival, ni, &(A[i*nfunc]), ustr) ;
  }

  /*need to weight columns of A by quadrature weights in long rule*/
  for ( i = 0 ; i < ni ; i ++ ) {
    gdouble dx ;
    dx = 0.5*(ival[i+1] - ival[i]) ;
    for ( j = 0 ; j < nq ; j ++ ) {
      for ( k = 0 ; k < nf ; k ++ ) {
	A[k*nfunc+i*nq+j] *= sqrt(gqr_rule_weight(rule, j)*dx) ;
      }
    }
  }
  
  /*perform RRQR*/
  job = 2 ;    /*no transpose operations (in FORTRAN indexing)*/
  n = nf ;     /*number of columns in A*/
  m = nfunc ;  /*number of rows in A*/
  k = nfunc ;     /*number of rows in C*/
  lwork = 2*m*n + 2*n + MAX(k, n) ;
  work = (gdouble *)g_malloc0(lwork*sizeof(gdouble)) ;
  dgeqpx_(&job, &m, &n, &k, A, &m, C, &k, jpvt, &ircond, &orcond,
	  &rank, svlues, work, &lwork, &info) ;

  /*check orthogonality*/

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gqr_rule_t *g, *gleg ;
  gqr_func_t func ;
  gqr_parameter_t p ;
  gqr_t type ;
  gpointer fdata[2] ;
  gdouble a, b, y, yy ;
  gdouble y0, y1, dy ;
  gdouble dx, xbar, x ;
  gdouble g0, g1, g2, glog, tol ;
  gint i, j, ns, ngp, ngps, N, M, n, nmax, nq, test, nx ;
  gchar ch, *progname ;
  gdouble *I ;
  gdouble P[3], Q[3], R[3], r[512], du[512], x0, x1, du0, s0, t0 ;
  gint np ;
  FILE *input ;
  
  N = 32 ;
  x = 0.4 ; y = 0.3 ;
  np = 64 ;

  N = 3 ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;

  a = 1.0 ; b = 3.0 ; tol = 1e-9 ; nq = 32 ;
  
  while ( (ch = getopt(argc, argv, "a:b:e:n:q:s:T:t:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'a': a = atof(optarg) ; break ;
    case 'b': b = atof(optarg) ; break ;
    case 'e': tol = atof(optarg) ; break ;
    case 'n': nx  = atoi(optarg) ; break ;
    case 'q': nq  = atoi(optarg) ; break ;
    case 's': s0  = atof(optarg) ; break ;
    case 'T': test = parse_test(optarg) ; break ;      
    case 't': t0  = atof(optarg) ; break ;
    }
  }

  if ( test == 0 ) {
    legendre_discretization_test(nq, a, b, tol) ;

    return 0 ;
  }
  
  if ( test == 1 ) {
    orthogonalization_test(nq, a, b, tol) ;
    
    return 0 ;
  }
  
  return 0 ;

  
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

  /* legendre_interp_matrix_test(8, 0.3, 0.9) ; */


  return 0 ;
  
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
