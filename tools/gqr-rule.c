#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"

gint fill_array(gchar *arg, gint **s, gint *ns) ;

gint fill_array(gchar *arg, gint **s, gint *ns)

{
  gint i ;
  gchar **tokens ;

  tokens = g_strsplit(arg, ",", 0) ;

  for ( (*ns) = 0 ; tokens[(*ns)] != NULL ; (*ns) ++ ) ;

  *s = (gint *)g_malloc((*ns)*sizeof(gint)) ;

  for ( i = 0 ; i < *ns ; i ++ ) (*s)[i] = atoi(tokens[i]) ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gqr_rule_t *g ;
  gqr_parameter_t p ;
  gdouble x, y, emax, Imax, a, b, tol ;
  gint N, M, *s, i, ns, imax ;
  gchar ch ;
  gboolean print_help, error_check ;
  gqr_t rule, baserule, singularity ;
  gchar *progname ;
  
  M = 4 ; N = 16 ; x = G_MAXDOUBLE ; y = G_MAXDOUBLE ;
  baserule = GQR_GAUSS_LEGENDRE ; singularity = GQR_GAUSS_REGULAR ;
  s = NULL ;
  a = -1 ; b = 1 ; tol = 1e-12 ;
  print_help = FALSE ; error_check = FALSE ;
  gqr_logging_init(stderr, argv[0], G_LOG_LEVEL_WARNING, NULL) ;

  gqr_parameter_clear(&p) ;
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  while ( (ch = getopt(argc, argv, "ha:b:CEef:GHi:JLM:N:p:Ps:Tt:x:y:Z"))
	  != EOF ) {
    switch(ch) {
    case 'h':
    default: 
      print_help = TRUE ;
      break ;
    case 'a': a = atof(optarg) ; break ;
    case 'b': b = atof(optarg) ; break ;
    case 'C': baserule = GQR_GAUSS_CHEBYSHEV_1 ; break ;
    /* case 'E': baserule = GQR_GAUSS_LAGUERRE ; break ; */
    case 'e': error_check = TRUE ; break ;
    case 'f': gqr_parameter_set_double(&p, atof(optarg)) ; break ;
    case 'G': baserule = GQR_GAUSS_GENERALIZED ; break ;
    case 'H': baserule = GQR_GAUSS_HERMITE ; break ;
    case 'i': gqr_parameter_set_int(&p, atoi(optarg)) ; break ;
    case 'J': baserule = GQR_GAUSS_JACOBI ; break ;
    case 'L': baserule = GQR_GAUSS_LEGENDRE ; break ;
    case 'M': M = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'T': baserule = GQR_GAUSS_CHEBYSHEV_2 ; break ;
    case 't': tol = atof(optarg) ; break ;
    case 'P': gqr_pointers_list(stderr, TRUE) ; return 0 ; break ;
    case 'p':
      gqr_parameter_set_pointer(&p, gqr_pointer_parse(optarg)) ;
      break ;
    case 's': fill_array(optarg, &s, &ns) ; break ;
    case 'x': x = atof(optarg) ; break ;
    case 'y': y = atof(optarg) ; break ;
    case 'Z': singularity = GQR_GAUSS_PAGET ; break ;
    }
  }

  if ( print_help || argc == 1) {
    fprintf(stderr, 
	    "%s: compute quadrature rules\n\n"
	    "Usage: %s <options>\n\n"
	    "Options:\n"
	    "        -h print this message and exit\n"
	    "        -a lower limit of integration in tests (%lg)\n"
	    "        -b upper limit of integration in tests (%lg)\n"
	    "        -C Gauss-(T)Chebyshev of the first kind quadrature\n"
	    /* "        -E Gauss-Laguerre (exponential) quadrature\n" */
	    "        -e estimate maximum error in quadrature\n"
	    "        -f # set float in parameter list\n"
	    "        -G generalized Gaussian quadrature (Bremer, Gimbutas and\n"
	    "           Rokhlin algorithm)\n"
	    "        -H Gauss-Hermite quadrature\n"
	    "        -i # set integer in parameter list\n"
	    "        -J Gauss-Jacobi quadrature (experimental, be careful)\n"
	    "        -L Gauss-Legendre quadrature (default)\n"
	    "        -M <order of polynomials to handle>\n"
	    "        -N <number of points in rule>\n"
	    "        -P list pointer names which can be parsed as parameters\n"
	    "        -p # set (parsed) pointer in parameter list\n"	    
	    "        -s <list of singularity orders>\n"	      
	    "        -T Gauss-(T)Chebyshev of the second kind quadrature\n"
	    "        -t # tolerance for error check (%lg)\n"
	    "        -x <x coordinate of singularity>\n"
	    "        -y <y coordinate of singularity>\n"
	    "        -Z singularity at zero (Paget quadrature rule)\n",
	    progname, progname, a, b, tol) ;
    return 0 ;
  }

  rule = baserule ;
  g = gqr_rule_alloc(N) ;

  if ( x != G_MAXDOUBLE ) {
    gqr_parameter_set_int(&p, M) ;
    gqr_parameter_set_double(&p, x) ;
    singularity = GQR_GAUSS_HYPERSINGULAR ;
  }
  if ( y != G_MAXDOUBLE ) {
    if ( p.ni != 1 ) gqr_parameter_set_int(&p, M) ;
    gqr_parameter_set_double(&p, y) ;
    singularity = GQR_GAUSS_HYPERSINGULAR ;
  }
  if ( s != NULL ) {
    for ( i = 0 ; i < ns ; i ++ ) 
      gqr_parameter_set_int(&p, s[i]) ;
    singularity = GQR_GAUSS_MULTISINGULAR ;
  }

  rule = baserule | singularity ;
  gqr_rule_select(g, rule, N, &p) ;
  gqr_rule_write(g, stdout) ;

  if ( !error_check ) return 0 ;

  fprintf(stderr, "%s: testing quadrature rule, tol = %lg\n",
	  progname, tol) ;
  
  error_check = gqr_test_rule(g, a, b, tol, &emax, &Imax, &imax) ;  

  if ( error_check )
    fprintf(stderr,
	    "  PASS: (emax = %lg, Imax = %lg, imax = %d, L_inf = %lg)\n",
	    emax, Imax, imax, emax/Imax) ;
  else
    fprintf(stderr,
	    "  FAIL: (emax = %lg, Imax = %lg, imax = %d, L_inf = %lg)\n",
	    emax, Imax, imax, emax/Imax) ;

  /* /\*check quadratures against analytical results, where available*\/ */
  /* gqr_rule_bgr_check(g, &p, &imax, &emax, analytic_check, stderr) ; */

  /* fprintf(stderr, "maximum error: %lg, basis function %d\n", emax, imax) ; */
  
  return 0 ;
}
