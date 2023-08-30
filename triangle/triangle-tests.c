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

gint dgesdd_(gchar *JOBZ, gint *M, gint *N, gdouble *A, gint *lda,
	     gdouble *S, gdouble *U, gint *ldu, gdouble *VT, gint *ldvt,
	     gdouble *work, gint *lwork, gint *iwork, gint *info) ;

typedef gint (*tri_quad_func_t)(gdouble *x0, gdouble *x1, gdouble *n,
				gdouble s, gdouble t, gpointer data,
				gdouble *f, gint nf) ;

gint invert2x2(gdouble *A, gdouble *Ai)

{
  gdouble det ;

  det = A[0]*A[3] - A[1]*A[2] ;

  Ai[0] =  A[3]/det ; Ai[1] = -A[1]/det ;
  Ai[2] = -A[2]/det ; Ai[3] =  A[0]/det ;

  return 0 ;
}

gdouble *wandzura_select(gint n)

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

gint wandzura_quad(tri_quad_func_t func, gdouble *x0,
		   gdouble *xt, gint xstr, gint nt, gdouble *Iq, gint nqi)

{
  gdouble x[3], n[3], J, *q, A, wt, f[16], s, t ;
  gint nw, i, iq ;

  nw = 175 ;

  q = wandzura_select(nw) ;

  memset(Iq, 0, nqi*sizeof(gdouble)) ;
  for ( i = 0 ; i < nw ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ;
    element_point_3d(xt, xstr, nt, s, t, x, n, &J) ;
    func(x0, x, n, s, t, NULL, f, nqi) ;
    wt = J*q[3*i+2] ;
    for ( iq = 0 ; iq < nqi ; iq ++ )
      Iq[iq] += wt*f[iq] ;
  }
  
  return 0 ;
}

gdouble quad_radial_func(gdouble d, gint i, gdouble *q, gint nq)

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

gdouble radial_check(gdouble d, gint N, gdouble *q, gint nq, gboolean print)

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

gdouble atan2_uw(gdouble y, gdouble x)

{
  gdouble th ;

  th = atan2(y, x) ;

  if ( th >= 0 ) return th ;

  return 2.0*M_PI + th ;
}

gint triangle_decomp(gdouble *x0, gdouble *x1, gdouble *x2,
		     gdouble *r0, gdouble *th0, gdouble *r, gdouble *th)

{
  gdouble r1, th1, r2, th2 ;

  r1 = sqrt((x0[0] - x1[0])*(x0[0] - x1[0]) +
	    (x0[1] - x1[1])*(x0[1] - x1[1])) ;
  r2 = sqrt((x0[0] - x2[0])*(x0[0] - x2[0]) +
	    (x0[1] - x2[1])*(x0[1] - x2[1])) ;

  if ( r1 > r2 ) {
    *r = r1 ; *r0 = r2/r1 ;

    *th = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) ;

    *th0 = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) - *th ;

    if ( *th0 < -M_PI ) *th0 += 2.0*M_PI ;
    
    /* fprintf(stderr, "r1 > r2") ; */
    
    return 0 ;
  }
  
  /* fprintf(stderr, "r1 < r2") ; */

  *r = r2 ; *r0 = r1/r2 ;
  
  *th = atan2_uw(x2[1] - x0[1], x2[0] - x0[0]) ;

  *th0 = atan2_uw(x1[1] - x0[1], x1[0] - x0[0]) - *th ;

  if ( *th0 > 0 ) *th0 -= 2.0*M_PI ;
  
  return 0 ;
}

gint triangle_mapping(gdouble *xt, gint xstr, gint nt, gdouble s, gdouble t,
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
  gint i, one = 1, two = 2, three = 3, lwork, info, iwork ;
  
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

gint triangle_quad_laplace(gdouble *x0, gdouble *x1, gdouble *n,
			   gdouble s, gdouble t, gpointer data,
			   gdouble *f, gint nf)

{
  gdouble R, r[3] ;

  r[0] = x1[0] - x0[0] ;
  r[1] = x1[1] - x0[1] ;
  r[2] = x1[2] - x0[2] ;
  
  R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]) ;
  
  f[0] = 0.25*M_1_PI/R*(1-s-t) ;
  f[1] = 0.25*M_1_PI/R*(  s  ) ;
  f[2] = 0.25*M_1_PI/R*(    t) ;

  /* f *= (r[0]*n[0] + r[1]*n[1] + r[2]*n[2])/R/R ; */
  
  return 0 ;
}

gint triangle_quad(tri_quad_func_t func,
		   gdouble *xti,
		   gdouble *xt, gint xstr, gint nt,
		   gdouble *A, gdouble s, gdouble t,
		   gdouble *x0, gdouble J0,
		   gdouble d, gint N,
		   gdouble *Iq, gint nqi)

/*
  xti:  mapped plane triangle
  xt:   physical space triangle
  xstr: stride in xt
  nt:   number of nodes on xt (3 or 6)
  A:    mapping matrix
  s,t:  singularity coordinates on unit triangle
  x0:   physical location of singularity
  J0:   Jacobian at singularity
  d:    radius of central disc
  N:    order of integration
  Iq:   integrals of func over triangle
  nqi:  number of functions to integrate
*/
  
{
  gdouble r0, th0, r, th, *qr, *qt, rr, ti, M, sgn ;
  gdouble x[3], s1, t1, J, n[3], wt, Aq, C, S, zero[2]={0,0} ;
  gdouble Ac, As, f[16] ;
  gint i, j, k, iq, idx[] = {0, 1, 2, 0}, nqr, nqt ;

  memset(Iq, 0, nqi*sizeof(gdouble)) ;
  
  /*quadrature on inner disc*/
  nqt = 3*(N+2) + 2 ;
  nqr = (N+1)/2 + 1 ;

  legendre_quadrature_select(N, &qr, &nqr) ;
  for ( i = 0 ; i < nqt ; i ++ ) {
    ti = 2.0*M_PI*i/nqt ;
    C = cos(ti) ; S = sin(ti) ;
    Ac = A[0]*C + A[1]*S ;
    As = A[2]*C + A[3]*S ;
    for ( j = 0 ; j < nqr ; j ++ ) {
      rr = 0.5*d*(1.0 + qr[2*j+0]) ;
      /*coordinates on unit triangle*/
      s1 = s + rr*Ac ; t1 = t + rr*As ;
      /*map to physical triangle*/
      element_point_3d(xt, xstr, nt, s1, t1, x, n, &J) ;
      /*quadrature weight*/
      wt = J/J0*rr*qr[2*j+1]*0.5*d*2.0*M_PI/nqt ;
      func(x0, x, n, s1, t1, NULL, f, nqi) ;
      for ( iq = 0 ; iq < nqi ; iq ++ )
	Iq[iq] += wt*f[iq] ;
    }
  }

  for ( k = 0 ; k < 3 ; k ++ ) {
    /*singularity is at origin on remapped triangle*/
    triangle_decomp(zero, &(xti[2*idx[k]]), &(xti[2*idx[k+1]]),
		    &r0, &th0, &r, &th) ;
    sgn = SIGN(th0) ; th0 = fabs(th0) ;
    /*select th quadrature*/
    angular_quadrature_select(N, r0, th0, &qt, &nqt) ;
    for ( i = 0 ; i < nqt ; i ++ ) {
      ti = th0*qt[i*2+0] ;
      C = cos(sgn*ti+th) ; S = sin(sgn*ti+th) ;
      Ac = A[0]*C + A[1]*S ;
      As = A[2]*C + A[3]*S ;
      /*radial limit and quadrature selection*/
      M = r0*sin(th0)/(r0*sin(th0-ti) + sin(ti)) ;
      radial_quadrature_select(N, d/M/r, &qr, &nqr) ;
      for ( j = 0 ; j < nqr ; j ++ ) {
      	rr = d + (M*r-d)*qr[2*j+0] ;
	/*coordinates on unit triangle*/
	s1 = s + rr*Ac ; t1 = t + rr*As ;
	/*map to physical triangle*/
	element_point_3d(xt, xstr, nt, s1, t1, x, n, &J) ;
	/*quadrature weight*/
	wt = (M*r-d)*qr[2*j+1]*qt[2*i+1]*rr*th0*J/J0 ;
	func(x0, x, n, s1, t1, NULL, f, nqi) ;
	for ( iq = 0 ; iq < nqi ; iq ++ )
	  Iq[iq] += wt*f[iq] ;	
      }
    }
  }
  
  return 0 ;
}

gint triangle_quad_test_3d(gdouble *xt, gint xstr, gint nt,
			   gdouble s, gdouble t, gint N)

{
  gdouble x0[3], J0, xti[64], A[4], Ai[4], d, G[3], dG[3], Iq[16], Iw[16] ;
  tri_quad_func_t func ;

  func = triangle_quad_laplace ;
  
  triangle_mapping(xt, xstr, nt, s, t, x0, &J0, A, xti) ;

  fprintf(stderr, "x0 = (%lg,%lg,%lg); J0 = %lg\n", x0[0], x0[1], x0[2], J0) ;

  d = 0.1 ;

  triangle_quad(func, xti, xt, xstr, nt, A, s, t, x0, J0, d, N, Iq, 3) ;
  fprintf(stderr, "single layer: %lg %lg %lg\n",
	  Iq[0], Iq[1], Iq[2]) ;

  wandzura_quad(func, x0, xt, xstr, nt, Iw, 3) ;
  fprintf(stderr, "Wandzura single layer: %lg %lg %lg\n",
  	  Iw[0], Iw[1], Iw[2]) ;

  if ( nt == 3 ) {
    newman_tri_shape(x0, &(xt[0*xstr]), &(xt[1*xstr]), &(xt[2*xstr]),
		     NULL, 0.0, G, dG) ;

    G[0]  *= -0.25*M_1_PI ;  G[1] *= -0.25*M_1_PI ;  G[2] *= -0.25*M_1_PI ;
    dG[0] *= -0.25*M_1_PI ; dG[1] *= -0.25*M_1_PI ; dG[2] *= -0.25*M_1_PI ;
    
    fprintf(stderr, "Newman single layer: %lg %lg %lg\n",
	    G[0], G[1], G[2]) ;
    fprintf(stderr, "error:               %e %e %e\n",
    	    fabs(G[0]-Iq[0]), fabs(G[1]-Iq[1]), fabs(G[2]-Iq[2])) ;
    fprintf(stderr, "Newman double layer: %lg %lg %lg\n",
	    dG[0], dG[1], dG[2]) ;
  }
  
  return 0 ;
}

gint triangle_quad_test(gdouble *xt, gint xstr, gdouble s0, gdouble t0,
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

gint triangle_decomp_test(gdouble *xt, gdouble *x0)

{
  gdouble r0, th0, r, th ;
  gint i, idx[] = {0, 1, 2, 0} ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    fprintf(stdout, "%lg %lg 0 0\n", xt[i*2+0], xt[i*2+1]) ;
  }
  
  fprintf(stdout, "%lg %lg 0 0\n", x0[0], x0[1]) ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    fprintf(stderr, "%d: ", i) ;
    triangle_decomp(x0, &(xt[2*idx[i]]), &(xt[2*idx[i+1]]),
		    &r0, &th0, &r, &th) ;
    fprintf(stderr, "\n") ;
    fprintf(stdout, "%lg %lg %lg %lg\n", r, th, r0, th0) ;
  }
  
  return 0 ;
}

gint radial_quadrature_test(gint N)

{
  gdouble *delta, eint, emax, d, *q ;
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
    eint = 0.0 ;
    for ( j = 0 ; j < nj ; j ++ ) {
      d = delta[i] + (delta[i+1] - delta[i])/(nj+1)*(j+1) ;
      
      radial_quadrature_select(N, d, &q, &nq) ;
      
      emax = radial_check(d, N, q, nq, FALSE) ;
      
      /* fprintf(stderr, "N=%d; d=%lg; emax=%lg\n", N, d, emax) ; */
      eint = MAX(eint, emax) ;
    }
    fprintf(stderr, "eint = %lg\n", eint) ;
  }

  return 0 ;
}

gdouble angular_check(gdouble d, gint N, gdouble *q, gint nq, gboolean print)

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

gint angular_quadrature_test(gint N)

{
  gdouble *dth, eint, emax, d, *q, r0 ;
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
  /* for ( i = 0 ; dth[i] != 3.14159 ; i ++ ) { */
  for ( i = 10 ; i < 11 ; i ++ ) {
    fprintf(stderr, "d = (%lg,%lg)\n", dth[i], dth[i+1]) ;
    eint = 0.0 ;
    for ( j = 0 ; j < nj ; j ++ ) {
      d = dth[i] + (dth[i+1] - dth[i])/(nj+1)*(j+1) ;
      
      angular_quadrature_select(N, r0, d, &q, &nq) ;
      
      /* emax = angular_check(d, N, q, nq, FALSE) ; */
      emax = angular_check(d, N, q, nq, TRUE) ;
      
      /* fprintf(stderr, "N=%d; d=%lg; emax=%lg\n", N, d, emax) ; */
      eint = MAX(eint, emax) ;
    }
    fprintf(stderr, "eint = %lg\n", eint) ;
  }

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gint N, xstr, nt ;
  gdouble d, xt[64], x0[2], s, t ;
  gchar ch ;
  
  xt[2*0+0] = -1.0 ; xt[2*0+1] = -0.5 ; 
  xt[2*1+0] =  1.0 ; xt[2*1+1] = -0.5 ; 
  xt[2*2+0] =  0.0 ; xt[2*2+1] =  0.75 ; 

  x0[0] = 0.1 ; x0[1] = -0.4 ;

  N = 4 ;  d = 0.1 ;

  while ( (ch = getopt(argc, argv, "N:")) != EOF ) {
    switch(ch) {
    case 'N': N = atoi(optarg) ; break ;
    }
  }

  /* angular_quadrature_test(N) ; */

  /* return 0 ; */
  
  /* radial_quadrature_test(N) ; */

  /* return 0 ; */

  xstr = 4 ; nt = 3 ;
  xt[0*xstr+0] = -0.3  ; xt[0*xstr+1] = 0.5  ; xt[0*xstr+2] =  0.4  ;
  xt[1*xstr+0] =  0.35 ; xt[1*xstr+1] = 0.55 ; xt[1*xstr+2] =  0.2  ;
  xt[2*xstr+0] =  0.0  ; xt[2*xstr+1] = 1.2  ; xt[2*xstr+2] =  0.3  ;
  xt[3*xstr+0] =  0.0  ; xt[3*xstr+1] = 0.5  ; xt[3*xstr+2] =  0.25 ;
  xt[4*xstr+0] =  0.2  ; xt[4*xstr+1] = 0.8  ; xt[4*xstr+2] =  0.25 ;
  xt[5*xstr+0] = -0.1  ; xt[5*xstr+1] = 0.8  ; xt[5*xstr+2] =  0.3  ;

  xt[0*xstr+2] = xt[1*xstr+2] = xt[2*xstr+2] = 0.0 ;
  
  s = 0.3 ; t = 0.2 ;
  
  /* triangle_quad_test_3d(xt, xstr, nt, s, t, N) ; */
  
  /* return 0 ; */
  
  d = 0.0 ;
  triangle_quad_test(xt, xstr, s, t, d, N) ;

  return 0 ;
  
  triangle_decomp_test(xt, x0) ;

  return 0 ;
}
