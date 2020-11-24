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
  gdouble x, y, emax ;
  gint N, M, *s, i, ns, imax ;
  gchar ch ;
  gboolean print_help, error_check, analytic_check ;
  gqr_t rule, baserule, singularity ;

  M = 4 ; N = 16 ; x = G_MAXDOUBLE ; y = G_MAXDOUBLE ;
  baserule = GQR_GAUSS_LEGENDRE ; singularity = GQR_GAUSS_REGULAR ;
  s = NULL ;
  print_help = FALSE ; error_check = FALSE ; analytic_check = FALSE ;
  gqr_logging_init(stderr, argv[0], G_LOG_LEVEL_WARNING, NULL) ;

  gqr_parameter_clear(&p) ;

  while ( (ch = getopt(argc, argv, "haCef:GHi:LM:N:p:Ps:Tx:y:")) != EOF ) {
    switch(ch) {
    case 'h':
    default: 
      print_help = TRUE ;
      break ;
    case 'a': analytic_check = TRUE ; break ;
    case 'C': baserule = GQR_GAUSS_CHEBYSHEV_1 ; break ;
    case 'e': error_check = TRUE ; break ;
    case 'f': gqr_parameter_set_double(&p, atof(optarg)) ; break ;
    case 'G': baserule = GQR_GAUSS_GENERALIZED ; break ;
    case 'H': baserule = GQR_GAUSS_HERMITE ; break ;
    case 'i': gqr_parameter_set_int(&p, atoi(optarg)) ; break ;
    case 'L': baserule = GQR_GAUSS_LEGENDRE ; break ;
    case 'M': M = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'T': baserule = GQR_GAUSS_CHEBYSHEV_2 ; break ;
    case 'P': gqr_pointers_list(stderr, TRUE) ; return 0 ; break ;
    case 'p':
      gqr_parameter_set_pointer(&p, gqr_pointer_parse(optarg)) ;
      break ;
    case 's': fill_array(optarg, &s, &ns) ; break ;
    case 'x': x = atof(optarg) ; break ;
    case 'y': y = atof(optarg) ; break ;
    }
  }

  if ( print_help || argc == 1) {
    fprintf(stderr, 
	    "%s: compute quadrature rules\n\n"
	    "Usage: %s <options>\n\n"
	    "Options:\n"
	    "        -h print this message and exit\n"
	    "        -C Gauss-(T)Chebyshev of the first kind quadrature\n"
	    "        -e calculate estimated error in quadrature of basis "
	    "functions\n"
	    "        -f # set float in parameter list\n"
	    "        -G generalized Gaussian quadrature\n"
	    "        -H Gauss-Hermite quadrature\n"
	    "        -i # set integer in parameter list\n"
	    "        -L Gauss-Legendre quadrature (default)\n"
	    "        -M <order of polynomials to handle>\n"
	    "        -N <number of points in rule>\n"
	    "        -P list pointer names which can be parsed as parameters\n"
	    "        -p # set (parsed) pointer in parameter list\n"	    
	    "        -s <list of singularity orders>\n"	      
	    "        -T Gauss-(T)Chebyshev of the second kind quadrature\n"
	    "        -x <x coordinate of singularity>\n"
	    "        -y <y coordinate of singularity>\n",
	    argv[0], argv[0]) ;
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

  if ( baserule != GQR_GAUSS_GENERALIZED ) return 0 ;

  if ( !error_check ) return 0 ;

  /*check quadratures against analytical results, where available*/
  gqr_rule_bgr_check(g, &p, &imax, &emax, analytic_check, stderr) ;

  fprintf(stderr, "maximum error: %lg, basis function %d\n", emax, imax) ;
  
  return 0 ;
}
