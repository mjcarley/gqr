#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"
#include "gqr-private.h"

gint _gqr_quad_2_r(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log_new(gint N, gdouble x, gdouble y, gdouble *I) ;
gint _gqr_quad_log(gint N, gdouble x, gdouble y, gdouble *I) ;

gchar *tests[] = {"discretization",
		  "orthogonalization",
		  ""} ;

gint rrqr(gdouble *A, gint m, gint n, gdouble *tau, gint *jpvt,
	  gdouble *work, gint lwork) ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}


static gdouble _grule_hpoly(gint n, gdouble x)

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

static gdouble func_exp(gdouble x, gint n, gpointer data)

{
  return exp(x) ;
}

/* static gdouble func_one(gdouble x, gint n, gpointer data) */

/* { */
/*   if ( n == 0 ) return 1.0 ; */
/*   return 0.0 ; */
/* } */

/* static gdouble func_Pn(gdouble x, gint q, gpointer data) */

/* { */
/*   gdouble P, t1, t2 ; */
/*   gint n, k ; */

/*   n = 5 ; */

/*   t1 = 1-x ; t2 = 1+x ; */
/*   for ( (P = 0.0), (k = q) ; k <= n ; k ++ ) { */
/*     P += (gsl_pow_int(-1,q)*gsl_pow_int(t1,k-q) + */
/* 	  gsl_pow_int(-1,n)*gsl_pow_int(t2,k-q))* */
/*       gsl_pow_int(-1,k)/gsl_sf_fact(k)*gsl_sf_fact(n+k)/gsl_sf_fact(n-k)/ */
/*       gsl_pow_int(2,k+1)/gsl_sf_fact(k-q) ; */
/*   } */

/*   return P ; */
/* } */

/* static gdouble func_1mt(gdouble x, gint n, gpointer data) */

/* { */
/*   gdouble m, pm, y ; */

/*   m = 4.0 ; pm = -1.0 ; */
/*   y = 1.0+pm*x ; */
/*   switch (n) { */
/*   case 0: return pow(y,m) ; break ; */
/*   case 1: return pm*m*pow(y,m-1) ; break ; */
/*   case 2: return m*(m-1)*pow(y,m-2) ; break ; */
/*   default: g_assert_not_reached() ; break ; */
/*   } */
  
/*   return 0.0 ; */
/* } */

/* static gdouble func_power(gdouble x, gint m, gint *n) */

/* { */
/*   switch (m) { */
/*   default: g_assert_not_reached() ; break ; */
/*   case 0: return gsl_pow_int(x,(*n)) ; break ; */
/*   case 1: if ( *n == 0) return 0.0 ; */
/*     return (*n)*gsl_pow_int(x,(*n)-1) ; */
/*     break ; */
/*   case 2: if ( *n <= 1 ) return 0.0 ; */
/*     return (*n)*(*n-1)*gsl_pow_int(x,(*n)-2) ; */
/*     break ; */
/*   } */

/*   return 0.0 ; */
/* } */

/* static gdouble func_krsin(gdouble x, gint m, gpointer data) */

/* { */
/*   switch (m) { */
/*   default: g_assert_not_reached() ; break ; */
/*   case 0: return (sin(2*x)+cos(3*x)) ; break ; */
/*   case 1: return (2*cos(2*x)-3*sin(3*x)) ; break ; */
/*   } */

/*   return 0.0 ; */
/* } */

/* static gdouble func_log_test(gdouble t, gint m, gpointer *data) */

/* { */
/*   gdouble x = *(gdouble *)data[0] ; */
/*   gint n = *(gint *)data[1] ; */
/*   gdouble f ; */

/*   f = 0.0 ; */
/*   switch (m) { */
/*   case 0: */
/*     if ( x == t ) f = 0.0 ; */
/*     else f = gsl_pow_int(t, n)*log(fabs(x-t)) ; */
/*     break ; */
/*   case 1: */
/*     f = gsl_pow_int(t, n)/(x-t) ; */
/*     break ; */
/*   case 2: */
/*     f = gsl_pow_int(t, n)/(x-t)/(x-t) ; */
/*     break ; */
/*   default: g_assert_not_reached() ; break ; */
/*   } */
/*   return f ; */
/* } */

/* static gdouble log_test(gdouble x, gint q, gint n) */

/* { */
/*   gint m, k ; */
/*   gdouble I, t, xn, xn1 ; */
/*   static gqr_rule_t *g = NULL ; */

/*   m = n/2 ; I = 0.0 ; */

/*   if ( g == NULL ) { */
/*     g = gqr_rule_alloc(1024) ; */
/*     gqr_rule_select(g, GQR_GAUSS_LEGENDRE, 1024, NULL) ; */
/*   } */

/*   switch (q) { */
/*   default: g_assert_not_reached() ; break ; */
/*   case 0: */
/*     if ( x == -1.0 ) I = (1-gsl_pow_int(-1,n+1))*M_LN2 ; */
/*     if ( x == 1.0 ) I = (1+gsl_pow_int(-1,n))*M_LN2 ; */

/*     if ( (x != -1.0) && (x != 1.0) ) { */
/*       xn = gsl_pow_int(x,n+1) ; */
/*       I = (1-xn)*log(fabs(1-x)) + (gsl_pow_int(-1,n)+xn)*log(fabs(1+x)) ; */
/*     } */

/*     if ( (n) == 2*m ) { /\*even n*\/ */
/*       for ( (k = 0), (t = 2.0) ; k <= m ; k ++ ) { */
/* 	I -= t/(gdouble)(2*m-2*k+1) ; */
/* 	t *= x*x ; */
/*       } */
/*       I /= (gdouble)(2*m+1) ; */
/*     } else { */
/*       for ( (k = 0), (t = 2.0*x) ; k <= m ; k ++ ) { */
/* 	I -= t/(gdouble)(2*m-2*k+1) ; */
/* 	t *= x*x ; */
/*       } */
/*       I /= (gdouble)(2*m+2) ; */
/*     } */
/*     break ; */
/*   case 1: */
/*     if ( (x == 1.0) || (x == -1.0) ) { */
/*       return gqr_finite_part_integral(func_power, &n, x, 1.0, -1, 1, g) ; */
/*     } */

/*     if ( (x != -1.0) && (x != 1.0) ) { */
/*       xn = gsl_pow_int(x,n) ; xn1 = gsl_pow_int(x,n+1) ; */
/*       I = -(1-xn1)/(1-x) - (n+1)*xn*log(fabs(1-x)) */
/* 	+ (gsl_pow_int(-1,n)+xn1)/(1+x) + (n+1)*xn*log(fabs(1+x)) ; */
/*     } */

/*     if ( (n) == 2*m ) { /\*even n*\/ */
/*       for ( (k = 1), (t = 2.0*x) ; k <= m ; k ++ ) { */
/* 	I -= 2*k*t/(gdouble)(2*m-2*k+1) ; */
/* 	t *= x*x ; */
/*       } */
/*       I /= (gdouble)(2*m+1) ; */
/*     } else { */
/*       for ( (k = 0), (t = 2.0) ; k <= m ; k ++ ) { */
/* 	I -= (2*k+1)*t/(gdouble)(2*m-2*k+1) ; */
/* 	t *= x*x ; */
/*       } */
/*       I /= (gdouble)(2*m+2) ; */
/*     } */
/*     break ; */
/*   case 2: */
/*     return gqr_finite_part_integral(func_power, &n, x, 2.0, -1, 1, g) ; */
/*     break ; */
/*   } */

/*   return I ; */
/* } */

static gdouble discretization_test_func(gdouble t, gint i, gqr_parameter_t *p)

{
  gdouble f ;
  
  f = sin(t) ;
  
  return f ;
}

static gint legendre_discretization_test(gint nq, gdouble x0, gdouble x1,
					 gdouble tol)

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

static gdouble orthogonalization_test_func(gdouble t, gint i, gqr_parameter_t *p)

{
  gdouble d, f ;

  d = 1e-7 ;

  f = (1.0-d)*t + d ;
  
  if ( i == 0 ) return 1.0/f ;

  return pow(f, i-1) ;
}

static gint orthogonalization_test(gint nq, gdouble x0, gdouble x1, gdouble tol)

{
  gdouble *ival, *Q, *R11, *u, x, *A, *work ;
  gdouble w, r[16384], z[16384] ;
  gqr_adapt_func_t func = orthogonalization_test_func ;
  gqr_rule_t *rule ;
  gint ni, nimax, nn, i, j, nf, k, lwork ;
  gint rank, rankmax, nfunc, *pvt, ldr ;
  gpointer data = NULL ;

  nf = 32 ;
  nimax = 16384*4 ; 

  ival = (gdouble *)g_malloc(nimax*sizeof(gdouble)) ;
  pvt = (gint *)g_malloc0(nimax*nq*sizeof(gint)) ;
  
  rankmax = nf ;

  fprintf(stderr, "orthogonalization test\n") ;
  fprintf(stderr, "======================\n") ;
  fprintf(stderr, "interval:   %lg, %lg\n", x0, x1) ;
  fprintf(stderr, "tolerance:  %lg\n",      tol) ;
  fprintf(stderr, "quadrature: %d\n",       nq) ;
  fprintf(stderr, "functions:  %d\n",       nf) ;
  
  lwork = 3*nimax*nq+1 + rankmax ;
  work = (gdouble *)g_malloc(lwork*sizeof(gdouble)) ;
  fprintf(stderr, "lwork:      %d\n",       lwork) ;
  
  ival[0] = x0 ; ival[1] = x1 ;
  
  rule = gqr_rule_alloc(nq) ;
  gqr_rule_select(rule, GQR_GAUSS_LEGENDRE, nq, NULL) ;

  Q   = (gdouble *)g_malloc(nq*nq*sizeof(gdouble)) ;
  R11 = (gdouble *)g_malloc(rankmax*rankmax*sizeof(gdouble)) ;

  gqr_discretize_interp(rule, Q) ;

  ni = 1 ;

  /*split interval for discretization*/
  for ( i = 0 ; i < nf ; i ++ ) {
    nn = 1 ;
    while ( nn != 0 ) {
      nn = gqr_discretize_adaptive(func, i, data, rule,
				   Q, ival, &ni, nimax, tol) ;
    }
  }

  fprintf(stderr, "intervals:  %d\n", ni) ;
  
  /*number of points in discretized functions*/
  nfunc = nq*ni ;
  
  /*use FORTRAN indexing from here*/
  /*fill function array*/
  u = (gdouble *)g_malloc(nfunc*rankmax*sizeof(gdouble)) ;
  A = (gdouble *)g_malloc(nfunc*     nf*sizeof(gdouble)) ;
  /*fill columns of matrix A with functions*/
  for ( i = 0 ; i < nf ; i ++ ) {
    gqr_discretize_fill(func, i, data, rule, ival, ni, &(A[i*nfunc]), 1) ;
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

  gqr_discretize_orthogonalize(A, nfunc, nf, tol, &rank, rankmax,
			       u, R11, pvt, &ldr, work, lwork) ;
  fprintf(stderr, "rank:       %d\n", rank) ;

  /*generate a quadrature rule*/
  memset(r, 0, nf*sizeof(gdouble)) ;
  for ( i = 0 ; i < ni ; i ++ ) {
    for ( j = 0 ; j < nq ; j ++ ) {
      gqr_discretize_quadrature(rule, ival, ni, i, j, &x, &w) ;
      for ( k = 0 ; k < nf ; k ++ ) {
	r[k] += u[k*nfunc+i*nq+j]*sqrt(w) ;
      }
    }
  }

  /*there should be some way to do the RRQR on the transpose without
    explicitly transposing, but I haven't found it yet*/
  for ( i = 0 ; i < rank ; i ++ ) {
    for ( j = 0 ; j < nfunc ; j ++ ) {
      A[j*rank+i] = u[i*nfunc+j] ;
    }
  }

  memset(pvt, 0, nfunc*sizeof(gint)) ;
  gqr_discretize_orthogonalize(A, rank, nfunc, 1e-15, &rank, rankmax,
  			       u, R11, pvt, &ldr, work, lwork) ;

  for ( i = 0 ; i < rank ; i ++ ) {
    z[i] = 0.0 ;
    for ( j = 0 ; j < rank ; j ++ ) {
      /*this is the transpose multiplication in FORTRAN indexing*/
      z[i] += u[i*rank+j]*r[j] ;
    }
  }

  i = rank-1 ;
  z[rank-1] = z[rank-1]/R11[i*rank+i] ;
  for ( i = rank-2 ; i >= 0 ; i -- ) {
    for ( j = i+1 ; j < rank ; j ++ ) {
      z[i] -= R11[j*rank+i]*z[j] ;
    }
    z[i] /= R11[i*rank+i] ;
  }

  for ( i = 0 ; i < rank ; i ++ ) {
    /*subtract one to switch from FORTRAN indexing*/
    gqr_discretize_quadrature(rule, ival, ni, pvt[i]-1, -1, &x, &w) ;
    fprintf(stdout, "%1.24e %1.24e\n", x, sqrt(w)*z[i]) ;
  }
  
  return 0 ;
#if 0  
  /*check orthogonality*/
  edmax = eoffmax = 0.0 ;
  for ( n = 0 ; n < rank ; n ++ ) {
    for ( m = 0 ; m < n ; m ++ ) {
      Iq = 0.0 ;

      for ( i = 0 ; i < ni ; i ++ ) {
	for ( j = 0 ; j < nq ; j ++ ) {
	  gqr_discretize_quadrature(rule, ival, ni, i, j, &x, &w) ;
	  Iq += u[n*nfunc+i*nq+j]*u[m*nfunc+i*nq+j] ;
	}
      }
      eoffmax = MAX(eoffmax, fabs(Iq)) ;
    }
    m = n ; 
    Iq = 0.0 ;
    
    for ( i = 0 ; i < ni ; i ++ ) {
      for ( j = 0 ; j < nq ; j ++ ) {
	gqr_discretize_quadrature(rule, ival, ni, i, j, &x, &w) ;
	Iq += u[n*nfunc+i*nq+j]*u[m*nfunc+i*nq+j] ;
      }
    }
    edmax = MAX(edmax, fabs(Iq-1.0)) ;

    for ( m = n+1 ; m < rank ; m ++ ) {
      Iq = 0.0 ;

      for ( i = 0 ; i < ni ; i ++ ) {
	for ( j = 0 ; j < nq ; j ++ ) {
	  gqr_discretize_quadrature(rule, ival, ni, i, j, &x, &w) ;
	  Iq += u[n*nfunc+i*nq+j]*u[m*nfunc+i*nq+j] ;
	}
      }
      eoffmax = MAX(eoffmax, fabs(Iq)) ;
    }
  }

  fprintf(stderr, "orthogonality error: %lg, %lg\n", edmax, eoffmax) ;
#endif
  
  return 0 ;
}

static gint adaptive_function_check(void)

{
  gint nj, nk, nimax, rmax, nf, nq ;
  gint i ;
  gqr_parameter_t p ;
  gqr_adapt_func_t func = grule_bgr_func_scattering_th ;
  
  nj = 4 ; nk = 5 ; nq = 16 ; rmax = 32 ; nimax = 8 ;

  nf = nj*(nj+3)*3/2 + nj*(3*nj+7)/2 + nk + nk-1 ;
  
  gqr_parameter_clear(&p) ;

  gqr_parameter_set_double(&p, 0.0) ;
  gqr_parameter_set_double(&p, 1.0) ;
  gqr_parameter_set_double(&p, 1e-9) ;
  gqr_parameter_set_double(&p, 0.7) ;
  gqr_parameter_set_double(&p, 1.5) ;

  gqr_parameter_set_int(&p, nf) ;
  gqr_parameter_set_int(&p, nq) ;
  gqr_parameter_set_int(&p, nimax) ;
  gqr_parameter_set_int(&p, rmax) ;
  gqr_parameter_set_int(&p, nj) ;
  gqr_parameter_set_int(&p, nk) ;

  for ( i = 0 ; i < nf ; i ++ ) {
    func(0.4, i, &p) ;
  }
  
  return 0 ;
}

static gdouble func_monomial(gdouble t,
			     gdouble *d, gint nd, gint *i, gint ni, gint j)

{
  gdouble f ;

  f = pow(t, j) ;
  
  return f ;
}

static gdouble beta_func(gdouble x, gdouble y)

/*
 * Gradshteyn and Ryzhik 8.384.1 
 */
  
{
  gdouble B ;

  B = tgamma(x)*tgamma(y)/tgamma(x+y) ;
  
  return B ;
}

static gint quad_monomial(gdouble a, gdouble b,
			  gdouble *d, gint nd, gint *i, gint ni,
			  gdouble *I, gint n)

{
  gint j ;

  for ( j = 0 ; j < n ; j ++ ) {
    I[j] = (pow(b, j+1) - pow(a, j+1))/(j+1) ;
  }
  
  return 0 ;
}

static gint quad_monomial_weighted(gdouble a, gdouble b,
				   gdouble *d, gint nd,
				   gint *i, gint ni,
				   gdouble *I, gint n)

/*
 *
 * integrals of x^n*(b-x)^\alpha*(x-a)^\beta
 * 
 * Gradshteyn and Ryzhik 3.196.3 and a recursion
 *
 */
  
{
  gdouble al, bt, mu, nu ;
  gint j ;
  
  al = d[0] ; bt = d[1] ;

  mu = al + 1.0 ; nu = bt + 1.0 ;
  I[0] = beta_func(mu, nu)*pow(b-a, mu+nu-1) ;
  I[1] = (mu*b + nu*a)*I[0]/(mu + nu) ;

  for ( j = 1 ; j < n ; j ++ ) {
    I[j+1] = ((a+b)*j + b*mu + a*nu)*I[j] - a*b*j*I[j-1] ;
    I[j+1] /= mu + nu + j ;
  }
  
  return 0 ;
}

static gint set_parameters(gqr_t baserule, gqr_t singularity,
			   gint N,
			   gdouble pdouble[], gint npd,
			   gint pint[], gint npi,
			   gqr_parameter_t *p,
			   gpointer *func, gpointer *quad, gint *nfunc,
			   gchar **errstr)

{
  switch ( baserule ) {
  default: g_assert_not_reached() ; break ;
  case GQR_GAUSS_LEGENDRE:
    *func = func_monomial ;
    *quad = quad_monomial ;
    *nfunc = 2*N-1 ;
    return 0 ;
  case GQR_GAUSS_JACOBI:
    if ( npd < 2 ) {
      *errstr = "Gauss-Jacobi rule requires at least two double parameters" ;
      return -1 ;
    }
    *func = func_monomial ;
    *quad = quad_monomial_weighted ;
    gqr_parameter_set_double(p, pdouble[0]) ;
    gqr_parameter_set_double(p, pdouble[1]) ;
    *nfunc = 2*N-1 ;
    return 0 ;
    break ;
  }

  return -1 ;
}

static gdouble quad_eval(gqr_rule_t *g, gdouble a, gdouble b,
			 gdouble pdouble[], gint npd,
			 gint pint[], gint npi,
			 gdouble (* func)(gdouble t, gdouble *d, gint nd,
					  gint *i, gint ni, gint j),
			 gint j)

{
  gdouble tbar, dt, I, t, W ;
  gint i ;
  
  gqr_rule_scale(g, a, b, &tbar, &dt, &W) ;

  I = 0.0 ;
  for ( i = 0 ; i < gqr_rule_length(g) ; i ++ ) {
    t = gqr_rule_abscissa(g,i)*dt + tbar ;
    I += W*gqr_rule_weight(g,i)*func(t, pdouble, npd, pint, npi, j) ;
  }
  
  return I ;
}

static gboolean quadrature_check(gqr_t baserule, gqr_t singularity,
				 gdouble a, gdouble b, gint N,
				 gdouble pdouble[], gint npd,
				 gint pint[], gint npi, gdouble tol,
				 gdouble *emax, gint *imax)

{
  gqr_rule_t *g ;
  gqr_parameter_t p ;
  gdouble (* func)(gdouble t, gdouble *d, gint nd, gint *i, gint ni, gint j) ;
  gint (* quad)(gdouble a, gdouble b,
		gdouble *d, gint nd, gint *i, gint ni, gdouble *I, gint n) ;
  gint nfunc, i ;
  gdouble Ic, Ia[128] ;
  gchar *errstr ;

  g = gqr_rule_alloc(N) ;
  gqr_parameter_clear(&p) ;

  if ( set_parameters(baserule, singularity, N, pdouble, npd, pint, npi,
		      &p, (gpointer *)(&func), (gpointer *)(&quad), &nfunc,
		      &errstr) != 0 ) {
    fprintf(stderr, "error: %s\n", errstr) ;
    *emax = 0 ; *imax = -1 ;
    return FALSE ;
  }

  gqr_rule_select(g, baserule | singularity, N, &p) ;

  *imax = 0 ; *emax = 0.0 ;
  quad(a, b, pdouble, npd, pint, npi, Ia, nfunc) ;
  for ( i = 0 ; i <= nfunc ; i ++ ) {
    Ic = quad_eval(g, a, b, pdouble, npd, pint, npi, func, i) ;
    if ( fabs(Ic - Ia[i]) > (*emax) ) {
      *imax = i ; *emax = fabs(Ic - Ia[i]) ;
    }
    /* fprintf(stderr, "%d %lg %lg\n", i, Ic, Ia[i]) ; */
  }

  if ( *emax <= tol ) return TRUE ;
  
  return FALSE ;
}

gint main(gint argc, gchar **argv)

{
  gdouble a, b, tol, pdouble[32], emax ;
  gint N, nq, test, pint[32], npi, npd, imax ;
  gchar ch, *progname ;
  gqr_t baserule, singularity ;
  gboolean pass ;
  
  N = 32 ;
  test = -1 ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  baserule = GQR_GAUSS_LEGENDRE ; singularity = GQR_GAUSS_REGULAR ;
  baserule = GQR_GAUSS_JACOBI ; singularity = GQR_GAUSS_REGULAR ;

  a = -1.0 ; b = 1.0 ; tol = 1e-9 ; nq = 32 ;
  npi = npd = 0 ;
  while ( (ch = getopt(argc, argv, "a:b:d:e:i:N:q:T:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'a': a = atof(optarg) ; break ;
    case 'b': b = atof(optarg) ; break ;
    case 'd': pdouble[npd] = atof(optarg) ; npd ++ ; break ;
    case 'e': tol = atof(optarg) ; break ;
    case 'i': pint[npi] = atoi(optarg) ; npi ++ ; break ;
    case 'N': N  = atoi(optarg) ; break ;
    case 'q': nq  = atoi(optarg) ; break ;
    case 'T': test = parse_test(optarg) ; break ;      
    }
  }

  pass = quadrature_check(baserule, singularity, a, b, N,
			  pdouble, npd, pint, npi, tol, &emax, &imax) ;

  if ( pass ) {
    fprintf(stderr, "%s: pass\n", progname) ;
  } else {
    fprintf(stderr, "%s: fail, maximum error %lg in integrand %d\n",
	    progname, emax, imax) ;    
  }
  
  if ( test == 0 ) {
    legendre_discretization_test(nq, a, b, tol) ;

    return 0 ;
  }
  
  if ( test == 1 ) {
    orthogonalization_test(nq, a, b, tol) ;
    
    return 0 ;
  adaptive_function_check() ;
  
  }

  return 0 ;
}
