/* Copyright (C) 2020 by Michael Carley */

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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <glib.h>

#define SIGN(_x) ((_x < 0) ? -1 : 1) 

extern gdouble WANDZURA_7[], WANDZURA_25[], WANDZURA_54[], WANDZURA_85[],
  WANDZURA_126[], WANDZURA_175[] ;

extern gdouble delta4[], delta8[], delta12[], delta16[] ;
extern gdouble dth_angular4[], dth_angular8[], dth_angular12[],
  dth_angular16[] ;

gint radial_quadrature_select(gint N, gdouble d, gdouble **q, gint *nq) ;
gint angular_quadrature_select(gint N, gdouble r0, gdouble th0,
			       gdouble **q, gint *nq) ;
gint legendre_quadrature_select(gint N, gdouble **q, gint *nq) ;
gint element_shape_3d(gint ne, gdouble s, gdouble t,
		      gdouble *L, gdouble *dLds, gdouble *dLdt) ;
gint element_point_interp_3d(gdouble *xe, gint xstr, gint ne,
			     gdouble *L, gdouble *dLds, gdouble *dLdt,
			     gdouble *y, gdouble *n, gdouble *J) ;
gint element_point_3d(gdouble *xe, gint xstr, gint ne,
		      gdouble s, gdouble t,
		      gdouble *y, gdouble *n, gdouble *J) ;

gint newman_tri_shape (gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		       gdouble imn[], gint hmax,
		       gdouble i[], gdouble j[]) ;

extern void dgesdd_(gchar *JOBZ, gint *M, gint *N, gdouble *A, gint *lda,
		    gdouble *S, gdouble *U, gint *ldu, gdouble *VT, gint *ldvt,
		    gdouble *work, gint *lwork, gint *iwork, gint *info) ;

typedef gint (*tri_quad_func_t)(gdouble *x0, gdouble *x1, gdouble *n,
				gdouble s, gdouble t, gpointer data,
				gdouble *f, gint nf) ;

static void invert2x2(gdouble *A, gdouble *Ai)

{
  gdouble det ;

  det = A[0]*A[3] - A[1]*A[2] ;

  Ai[0] =  A[3]/det ; Ai[1] = -A[1]/det ;
  Ai[2] = -A[2]/det ; Ai[3] =  A[0]/det ;

  return ;
}

static gdouble *wandzura_select(gint n)

{
  switch (n) {
  default: g_assert_not_reached() ; break ;
  case   7: return WANDZURA_7 ; break ;
  case  25: return WANDZURA_25 ; break ;
  case  54: return WANDZURA_54 ; break ;
  case  85: return WANDZURA_85 ; break ;
  case 126: return WANDZURA_126 ; break ;
  case 175: return WANDZURA_175 ; break ;
  }
  
  return NULL ;
}

static gdouble quad_radial_func(gdouble d, gint i, gdouble *q, gint nq)

{
  gdouble f, df, t, w ;
  gint j ;

  f = 0.0 ;
  for ( j = 0 ; j < nq ; j ++ ) {
    t = q[2*j+0] ; w = q[2*j+1] ;
    df = (1.0 - d)*t + d ;
    f += pow(df, i)*w ;
  }
  
  return f ;
}

static gdouble radial_check(gdouble d, gint N, gdouble *q, gint nq,
			    gboolean print)

{
  gdouble emax = 0.0, Iq, Ia, x0, x1 ;
  gint i ;

  x0 = 0.0 ; x1 = 1.0 ;
  Ia = (log((1.0-d)*x1+d) - log((1.0-d)*x0+d))/(1.0-d) ;
  Iq = quad_radial_func(d, -1, q, nq) ;

  emax = fabs(Ia-Iq) ;
  if ( print ) 
    fprintf(stderr, "-1 %lg %lg %e\n", Ia, Iq, fabs(Ia-Iq)) ;

  for ( i = 0 ; i < 2*N-2 ; i ++ ) {
    Iq = quad_radial_func(d, i, q, nq) ;
    Ia = (pow((1.0-d)*x1+d,i+1) -
	  pow((1.0-d)*x0+d,i+1))/(gdouble)(i+1)/(1.0-d) ;

    emax = MAX(emax, fabs(Ia-Iq)) ;
    if ( print ) 
      fprintf(stderr, "%d %lg %lg %e\n", i, Ia, Iq, fabs(Ia-Iq)) ;
  }
  
  return emax ;
}

static gdouble atan2_uw(gdouble y, gdouble x)

{
  gdouble th ;

  th = atan2(y, x) ;

  if ( th >= 0 ) return th ;

  return 2.0*M_PI + th ;
}

static gint triangle_decomp(gdouble *x0, gdouble *x1, gdouble *x2,
			    gdouble *r0, gdouble *th0, gdouble *r, gdouble *th)

{
  gdouble r1, r2 ;

  r1 = sqrt((x0[0] - x1[0])*(x0[0] - x1[0]) +
	    (x0[1] - x1[1])*(x0[1] - x1[1])) ;
  r2 = sqrt((x0[0] - x2[0])*(x0[0] - x2[0]) +
	    (x0[1] - x2[1])*(x0[1] - x2[1])) ;

  if ( r1 > r2 ) {
    *r = r1 ; *r0 = r2/r1 ;

    *th = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) ;

    *th0 = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) - *th ;

    if ( *th0 < -M_PI ) *th0 += 2.0*M_PI ;
    
    return 0 ;
  }

  *r = r2 ; *r0 = r1/r2 ;
  
  *th = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) ;

  *th0 = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) - *th ;

  if ( *th0 > 0 ) *th0 -= 2.0*M_PI ;
  
  return 0 ;
}

static gint triangle_mapping(gdouble *xt, gint xstr, gint nt,
			     gdouble s, gdouble t,
			     gdouble *x0, gdouble *J0, gdouble *A, gdouble *xti)


/*
  mapping information for transformation of unit triangle to physical
  element via remapped triangle with conformal mapping at singularity

  xt:    nodes of physical triangle
  xstr:  stride in xt
  nt:    number of nodes on xt (3 or 6)
  s, t:  coordinates of singular point on unit triangle
  x0:    calculated physical location of singular point
  J0:    reference Jacobian at singularity
  A:     transformation matrix for Affine mapping
  xti:   triangle for integration using singular quadrature rules

  integrate on xti and map xti to coordinates (s1,t1) on unit triangle
  with

  (s1,t1) = [A]*r[\cos\theta \sin\theta]' + (s,t)'

  map to physical triangle using shape functions as usual

  x = L(s1,t1)*[xt]
*/
  
{
  gdouble L[16], dLds[16], dLdt[16], n[3], U[16], S[16], V[16] ;
  gdouble dr[6], work[138], Ai[4] ;
  gint i, two = 2, three = 3, lwork, info, iwork ;
  
  /*physical location of singularity and reference Jacobian*/
  element_shape_3d(nt, s, t, L, dLds, dLdt) ;
  element_point_interp_3d(xt, xstr, nt, L, dLds, dLdt, x0, n, J0) ;

  /*SVD of Jacobian matrix*/
  for ( i = 0 ; i < nt ; i ++ ) {
    dr[0] += dLds[i]*xt[i*xstr+0] ;
    dr[1] += dLds[i]*xt[i*xstr+1] ;
    dr[2] += dLds[i]*xt[i*xstr+2] ;
    dr[3] += dLdt[i]*xt[i*xstr+0] ;
    dr[4] += dLdt[i]*xt[i*xstr+1] ;
    dr[5] += dLdt[i]*xt[i*xstr+2] ;
  }

  /*lwork = 138 from a workspace query: size of problem is fixed so
    lwork is too*/
  lwork = 138 ; info = 1 ; iwork = 1 ;
  dgesdd_("A", &three, &two, dr, &three, S, U, &three, V, &two,
	  work, &lwork, &iwork, &info) ;
  
  /*A = V*inv(S) = V*[1/S(1) 0; 0 1/S(2)*/
  A[0] = V[0]/S[0] ; A[1] = V[1]/S[1] ;
  A[2] = V[2]/S[0] ; A[3] = V[3]/S[1] ;
  invert2x2(A, Ai) ;

  /*remapped unit triangle*/
  xti[0] = Ai[0]*(0 - s) + Ai[1]*(0 - t) ;
  xti[1] = Ai[2]*(0 - s) + Ai[3]*(0 - t) ;
  xti[2] = Ai[0]*(1 - s) + Ai[1]*(0 - t) ;
  xti[3] = Ai[2]*(1 - s) + Ai[3]*(0 - t) ;
  xti[4] = Ai[0]*(0 - s) + Ai[1]*(1 - t) ;
  xti[5] = Ai[2]*(0 - s) + Ai[3]*(1 - t) ;
  
  return 0 ;
}



/* static gint triangle_quad_test_3d(gdouble *xt, gint xstr, gint nt, */
/* 			   gdouble s, gdouble t, gint N) */

/* { */
/*   gdouble x0[3], J0, xti[64], A[4], d, G[3], dG[3], Iq[16], Iw[16] ; */
/*   tri_quad_func_t func ; */

/*   func = triangle_quad_laplace ; */
  
/*   triangle_mapping(xt, xstr, nt, s, t, x0, &J0, A, xti) ; */

/*   fprintf(stderr, "x0 = (%lg,%lg,%lg); J0 = %lg\n", x0[0], x0[1], x0[2], J0) ; */

/*   d = 0.1 ; */

/*   triangle_quad(func, xti, xt, xstr, nt, A, s, t, x0, J0, d, N, Iq, 3) ; */
/*   fprintf(stderr, "single layer: %lg %lg %lg\n", */
/* 	  Iq[0], Iq[1], Iq[2]) ; */

/*   wandzura_quad(func, x0, xt, xstr, nt, Iw, 3) ; */
/*   fprintf(stderr, "Wandzura single layer: %lg %lg %lg\n", */
/*   	  Iw[0], Iw[1], Iw[2]) ; */

/*   if ( nt == 3 ) { */
/*     newman_tri_shape(x0, &(xt[0*xstr]), &(xt[1*xstr]), &(xt[2*xstr]), */
/* 		     NULL, 0.0, G, dG) ; */

/*     G[0]  *= -0.25*M_1_PI ;  G[1] *= -0.25*M_1_PI ;  G[2] *= -0.25*M_1_PI ; */
/*     dG[0] *= -0.25*M_1_PI ; dG[1] *= -0.25*M_1_PI ; dG[2] *= -0.25*M_1_PI ; */
    
/*     fprintf(stderr, "Newman single layer: %lg %lg %lg\n", */
/* 	    G[0], G[1], G[2]) ; */
/*     fprintf(stderr, "error:               %e %e %e\n", */
/*     	    fabs(G[0]-Iq[0]), fabs(G[1]-Iq[1]), fabs(G[2]-Iq[2])) ; */
/*     fprintf(stderr, "Newman double layer: %lg %lg %lg\n", */
/* 	    dG[0], dG[1], dG[2]) ; */
/*   } */
  
/*   return 0 ; */
/* } */

static gint triangle_quad_test(gdouble *xt, gint xstr, gdouble s0, gdouble t0,
			       gdouble d, gint N)

{
  gdouble r0, th0, r, th, *qr, *qt, rr, t, M, sgn, c, f ;
  gdouble x, y, J, A[2048], Aw[2048], x0[3], xti[64], J0 ;
  gint i, j, k, ip, n, m, Np, idx[] = {0, 1, 2, 0}, nqr, nqt, nw ;
  
  triangle_mapping(xt, xstr, 3, s0, t0, x0, &J0, A, xti) ;
  
  memset(A , 0, 2048*sizeof(gdouble)) ;
  memset(Aw, 0, 2048*sizeof(gdouble)) ;
  nw = 175 ;

  Np = 8 ;
  /*quadrature on inner disc*/
  nqt = 3*(N+2) + 2 ;
  nqr = (N+1)/2 + 1 ;

  legendre_quadrature_select(N, &qr, &nqr) ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    t = 2.0*M_PI*i/nqt ;
    for ( j = 0 ; j < nqr ; j ++ ) {
      rr = 0.5*d + 0.5*d*qr[2*j+0] ;
      x = x0[0] + rr*cos(t) ;
      y = x0[1] + rr*sin(t) ;
      for ( ip = 0, n = 0 ; n <= Np ; n ++ ) {
	for ( m = 0 ; m <= n ; m ++, ip ++ ) {
	  f = pow(x,m)*pow(y,n-m) ;
	  A[ip] += rr*0.5*d*qr[2*j+1]*2.0*M_PI/nqt*f ;
	}
      }
    }
  }
  
  for ( k = 0 ; k < 3 ; k ++ ) {
    triangle_decomp(x0, &(xt[idx[k]*xstr]), &(xt[idx[k+1]*xstr]),
		    &r0, &th0, &r, &th) ;
    sgn = SIGN(th0) ;
    th0 = fabs(th0) ;
    /*select th quadrature*/
    angular_quadrature_select(N, r0, th0, &qt, &nqt) ;
    c = 0.0 ;
    for ( i = 0 ; i < nqt ; i ++ ) {
      t = th0*qt[i*2+0] ;
      c += th0*qt[2*i+1] ;
      /*radial limit and quadrature selection*/
      M = r0*sin(th0)/(r0*sin(th0-t) + sin(t)) ;
      radial_quadrature_select(N, d/M/r, &qr, &nqr) ;
      for ( j = 0 ; j < nqr ; j ++ ) {
      	rr = d/r + (M-d/r)*qr[2*j+0] ;
	x = x0[0] + r*rr*cos(sgn*t+th) ;
	y = x0[1] + r*rr*sin(sgn*t+th) ;
	for ( ip = 0, n = 0 ; n <= Np ; n ++ ) {
	  for ( m = 0 ; m <= n ; m ++, ip ++ ) {
	    f = pow(x,m)*pow(y,n-m) ;
	    A[ip] += (M-d/r)*qr[2*j+1]*qt[2*i+1]*rr*r*r*(th0)*f ;
	  }
	}
      }
    }
  }

  /*reference quadrature*/
  qt = wandzura_select(nw) ;
  /*Jacobian is constant on the plane triangle*/
  J =
    (xt[xstr*0+0] - xt[xstr*2+0])*(xt[xstr*1+1] -  xt[xstr*2+1]) -
    (xt[xstr*1+0] - xt[xstr*2+0])*(xt[xstr*0+1] -  xt[xstr*2+1]) ;    
  
  for ( i = 0 ; i < nw ; i ++ ) {
    x = xt[xstr*0+0]*qt[3*i+0] +
      xt[xstr*1+0]*qt[3*i+1] +
      xt[xstr*2+0]*(1.0 - qt[3*i+0] - qt[3*i+1]) ;
    y = xt[xstr*0+1]*qt[3*i+0] +
      xt[xstr*1+1]*qt[3*i+1] +
      xt[xstr*2+1]*(1.0 - qt[3*i+0] - qt[3*i+1]) ;
      for ( ip = 0, n = 0 ; n <= Np ; n ++ ) {
	for ( m = 0 ; m <= n ; m ++, ip ++ ) {
	  f = pow(x,m)*pow(y,n-m) ;
	  Aw[ip] += qt[3*i+2]*J*f ;
	}
      }
  }

  for ( ip = 0, n = 0 ; n <= Np ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++, ip ++ ) {
      fprintf(stderr, "m=%d, n=%d, A = %lg, error = %lg\n",
	      m, n, A[ip], fabs(A[ip]-Aw[ip])) ;
    }
  }
  
  return 0 ;
}


static gint radial_quadrature_test(gint N)

{
  gdouble *delta, eint, emax, emin, d, *q ;
  gint i, j, nj, nq ;
  
  switch ( N ) {
  default: g_assert_not_reached() ;
  case 4:  delta = delta4 ; break ;
  case 8:  delta = delta8 ; break ;
  case 12:  delta = delta12 ; break ;
  case 16:  delta = delta16 ; break ;
  }

  nj = 15 ;
  for ( i = 0 ; delta[i] != 1.0 ; i ++ ) {
    fprintf(stderr, "d = (%lg,%lg)\n", delta[i], delta[i+1]) ;
    emax = 0.0 ; emin = G_MAXDOUBLE ;
    for ( j = 0 ; j < nj ; j ++ ) {
      d = delta[i] + (delta[i+1] - delta[i])/(nj+1)*(j+1) ;
      
      radial_quadrature_select(N, d, &q, &nq) ;
      
      eint = radial_check(d, N, q, nq, FALSE) ;
      
      emax = MAX(eint, emax) ;
      emin = MIN(eint, emin) ;
    }
    fprintf(stderr, "  error(min, max) = (%lg, %lg)\n", emin, emax) ;
  }

  return 0 ;
}

static gdouble angular_check(gdouble d, gint N, gdouble *q, gint nq,
			     gboolean print)

{
  gdouble emax = 0.0, Ic, Is, Ia, t0, t1, t ;
  gint i, j, nj ;

  t0 = 0.0 ; t1 = d ;

  nj = 3*(N+2)+3 ;

  /*sines and cosines*/
  Ic = 0.0 ;
  for ( j = 0 ; j < nq ; j ++ ) {
    Ic += (t1 - t0)*q[j*2+1] ;    
  }

  Ia = t1 - t0 ;
  emax = fabs(Ia-Ic) ;

  if ( print ) 
    fprintf(stderr, "0 %lg %lg %e\n", Ia, Ic, fabs(Ia-Ic)) ;

  for ( i = 1 ; i < nj ; i ++ ) {
    Ic = Is = 0.0 ;
    for ( j = 0 ; j < nq ; j ++ ) {
      t = (t1 - t0)*q[2*j+0] ;
      Ic += cos(i*t)*q[j*2+1]*(t1 - t0) ;
      Is += sin(i*t)*q[j*2+1]*(t1 - t0) ;
    }
    Ia = (sin(i*t1) - sin(i*t0))/i ;
    emax = MAX(emax, fabs(Ia-Ic)) ;
    if ( print ) 
      fprintf(stderr, "%d %lg %lg %e\n", i, Ia, Ic, fabs(Ia-Ic)) ;
    Ia = -(cos(i*t1) - cos(i*t0))/i ;
    emax = MAX(emax, fabs(Ia-Is)) ;
    if ( print ) 
      fprintf(stderr, "%d %lg %lg %e\n", i, Ia, Is, fabs(Ia-Is)) ;
  }
  
  return emax ;
}

static gint angular_quadrature_test(gint N)

{
  gdouble *dth, eint, emax, emin, d, *q, r0 ;
  gint i, j, nj, nq ;
  
  switch ( N ) {
  default: g_assert_not_reached() ;
  case 4:
    dth = dth_angular4 ;
    break ;
  case 8:
    dth = dth_angular8 ;
    break ;
  case 12:
    dth = dth_angular12 ;
    break ;
  case 16:
    dth = dth_angular16 ;
    break ;
  }

  r0 = 0.5 ;
  
  nj = 9 ;
  for ( i = 0 ; dth[i] != 3.14159 ; i ++ ) {
  /* for ( i = 10 ; i < 11 ; i ++ ) { */
    fprintf(stderr, "d = (%lg,%lg)\n", dth[i], dth[i+1]) ;
    emax = 0.0 ; emin = G_MAXDOUBLE ;
    for ( j = 0 ; j < nj ; j ++ ) {
      d = dth[i] + (dth[i+1] - dth[i])/(nj+1)*(j+1) ;
      
      angular_quadrature_select(N, r0, d, &q, &nq) ;
      
      /* emax = angular_check(d, N, q, nq, FALSE) ; */
      eint = angular_check(d, N, q, nq, FALSE) ;
      
      /* fprintf(stderr, "N=%d; d=%lg; emax=%lg\n", N, d, emax) ; */
      emax = MAX(eint, emax) ;
      emin = MIN(eint, emin) ;
    }
    fprintf(stderr, "  error(min, max) = (%lg, %lg)\n", emin, emax) ;
  }

  return 0 ;
}

static void print_help_message(FILE *f, gchar *progname, gint N)

{
  fprintf(f, "%s: testing quadrature rules generated for boundary integral\n\n",
	  progname) ;

  fprintf(f,
	  "Usage: %s [options]\n\n"
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -A check angular quadrature rules\n"
	  "  -N # order of quadrature rule (%d)\n"
	  "  -R check radial quadrature rules\n\n",
	  progname, N) ;

  fprintf(f,
	  "If neither -A nor -R is selected, the test is performed on an\n"
	  "arbitrary triangle, using Wandzura and Xu's high order rules\n"
	  "for comparison\n\n") ;

  fprintf(f,
	  "Angular and radial rule checks are performed by comparing\n"
	  "numerically evaluated and analytical results for integrals.\n"
	  "Results are output as\n\n"
	  "  i Ia Ic err\n\n"
	  "with Ia and Ic analytical and numerical results and"
	  "  err=abs(Ia-Ic).\n\n"
	  "Functions used for evaluation are\n\n"
	  "  [(1-d)*t + d]^i (radial)\n\n"
	  "and\n\n"
	  "  cos(i*t), sin(i*t) (angular)\n\n"
	  "over the appropriate integration range.\n") ;
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gint N, xstr ;
  gdouble d, xt[64], s, t ;
  gchar ch ;
  
  xt[2*0+0] = -1.0 ; xt[2*0+1] = -0.5 ; 
  xt[2*1+0] =  1.0 ; xt[2*1+1] = -0.5 ; 
  xt[2*2+0] =  0.0 ; xt[2*2+1] =  0.75 ; 

  N = 4 ;  d = 0.1 ;

  while ( (ch = getopt(argc, argv, "hAN:R")) != EOF ) {
    switch(ch) {
    case 'h':
      print_help_message(stderr, argv[0], N) ;
      return 0 ;
      break ;
    case 'A':
      angular_quadrature_test(N) ;
      return 0 ;
      break ;
    case 'N': N = atoi(optarg) ; break ;
    case 'R':
      radial_quadrature_test(N) ;
      return 0 ;
      break ;
    }
  }

  xstr = 4 ;
  xt[0*xstr+0] = -0.3  ; xt[0*xstr+1] = 0.5  ; xt[0*xstr+2] =  0.4  ;
  xt[1*xstr+0] =  0.35 ; xt[1*xstr+1] = 0.55 ; xt[1*xstr+2] =  0.2  ;
  xt[2*xstr+0] =  0.0  ; xt[2*xstr+1] = 1.2  ; xt[2*xstr+2] =  0.3  ;
  xt[3*xstr+0] =  0.0  ; xt[3*xstr+1] = 0.5  ; xt[3*xstr+2] =  0.25 ;
  xt[4*xstr+0] =  0.2  ; xt[4*xstr+1] = 0.8  ; xt[4*xstr+2] =  0.25 ;
  xt[5*xstr+0] = -0.1  ; xt[5*xstr+1] = 0.8  ; xt[5*xstr+2] =  0.3  ;

  xt[0*xstr+2] = xt[1*xstr+2] = xt[2*xstr+2] = 0.0 ;
  
  s = 0.3 ; t = 0.2 ;
  
  d = 0.0 ;
  triangle_quad_test(xt, xstr, s, t, d, N) ;

  return 0 ;
}
