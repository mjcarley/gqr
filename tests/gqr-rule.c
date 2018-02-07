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
  gdouble x, y ;
  gint N, M ;
  gint *s, i, ns ;
  gchar ch ;
  gboolean print_help ;
  gqr_t rule, baserule, singularity ;

  M = 4 ; N = 16 ; x = G_MAXDOUBLE ; y = G_MAXDOUBLE ;
  baserule = GQR_GAUSS_LEGENDRE ; singularity = GQR_GAUSS_REGULAR ;
  s = NULL ; print_help = FALSE ;
  gqr_logging_init(stderr, argv[0], G_LOG_LEVEL_WARNING, NULL) ;

  while ( (ch = getopt(argc, argv, "hCHLM:N:s:Tx:y:")) != EOF ) {
    switch(ch) {
    case 'h':
    default: 
      print_help = TRUE ;
      break ;
    case 'C': baserule = GQR_GAUSS_CHEBYSHEV_1 ; break ;
    case 'H': baserule = GQR_GAUSS_HERMITE ; break ;
    case 'L': baserule = GQR_GAUSS_LEGENDRE ; break ;
    case 'M': M = atoi(optarg) ; break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'T': baserule = GQR_GAUSS_CHEBYSHEV_2 ; break ;
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
	    "        -H Gauss-Hermite quadrature\n"
	    "        -L Gauss-Legendre quadrature (default)\n"
	    "        -M <order of polynomials to handle>\n"
	    "        -N <number of points in rule>\n"
	    "        -s <list of singularity orders>\n"	      
	    "        -T Gauss-(T)Chebyshev of the second kind quadrature\n"
	    "        -x <x coordinate of singularity>\n"
	    "        -y <y coordinate of singularity>\n",
	    argv[0], argv[0]) ;
    return 0 ;
  }

  rule = baserule ;
  g = gqr_rule_alloc(N) ;

  gqr_parameter_clear(&p) ;
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

  return 0 ;
}
