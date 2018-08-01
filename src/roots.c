/* Copyright (C) 2007 by  Michael Carley */

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

#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H */

#include "gqr.h"
#include "gqr-private.h"

/**
 * @defgroup roots Root-finding for abscissa of quadrature rules
 * @{
 * 
 */

#if 0
#ifndef M_PI
#define M_PI		3.14159265358979323846  /* pi */
#endif /*M_PI */
#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923  /* pi/2 */
#endif /*M_PI_2 */
#endif

#define GQR_POLYVAL(_P,_x) ((_P)[0] + (_x)*((_P)[1] + (_x)*(_P)[2]))
#define GQR_POLYDER(_P,_x) (_P[1] + 2.0*(_x)*_P[2])
#define GQR_POLYDER2(_P,_x) (2.0*_P[2])

/*this is recommended for double precision*/
#ifndef GQR_TAYLOR_LENGTH
#define GQR_TAYLOR_LENGTH 33
#endif /*GQR_TAYLOR_LENGTH */

const gdouble GQR_INVERSE_FACTORIALS[] =
  { 1.00000000000000000000e+00, 1.00000000000000000000e+00,
  5.00000000000000000000e-01, 1.66666666666666657415e-01,
  4.16666666666666643537e-02, 8.33333333333333321769e-03,
  1.38888888888888894189e-03, 1.98412698412698412526e-04,
  2.48015873015873015658e-05, 2.75573192239858925110e-06,
  2.75573192239858247483e-07, 2.50521083854417003715e-08,
  2.08767569878681291379e-09, 1.60590438368217038137e-10,
  1.14707455977297018891e-11, 7.64716373181982650293e-13,
  4.77947733238738967107e-14, 2.81145725434553558647e-15,
  1.56192069685862400623e-16, 8.22063524662432487331e-18,
  4.11031762331216484407e-19, 1.95729410633911301875e-20,
  8.89679139245060605134e-22, 3.86817017063067119571e-23,
  1.61173757109611912058e-24, 6.44695028438444029913e-26,
  2.47959626322479902454e-27, 9.18368986379561663084e-29,
  3.27988923706982094340e-30, 1.13099628864478144106e-31,
  3.76998762881592933312e-33, 1.21612504155352217019e-34,
  3.80039075485472203590e-36, 1.15163356207718068118e-37,
  3.38715753552115200837e-39, 9.67759295863188659354e-41,
  2.68822026628663792623e-42, 7.26546017915316893590e-44,
  1.91196320504028755368e-45, 4.90246975651360672149e-47
};

static void gqr_polyval(gdouble * P, gint n, gdouble x,
                        gdouble * p, gdouble * dp)
{
  gint i;

  for ((i = n), (*p = P[n]), (*dp = n * P[n]); i >= 1;
       (*p = x * (*p) + P[i]), (*dp = x * (*dp) + i * P[i]), (i--));
  *p = x * (*p) + P[0];

  return;
}

static void gqr_runge_kutta_solve(gdouble th0, gdouble th1, gint ns,
                                  gdouble * P, gdouble * Q, gdouble * R,
                                  gdouble * x)
{
  gdouble p, q, r, dp, dq, dr, k0, k1, h, th;
  gint i;

  h = (th1 - th0) / (gdouble) ns;

  p = GQR_POLYVAL(P, (*x));
  q = GQR_POLYVAL(Q, (*x));
  r = GQR_POLYVAL(R, (*x));
  dp = GQR_POLYDER(P, (*x));
  dq = GQR_POLYDER(Q, (*x));
  dr = GQR_POLYDER(R, (*x));

  k0 =
    -h * p * r / (r * sqrt(p * r) +
                  0.25 * (dr * p - dp * r + 2 * r * q) * sin(2 * th0));
  for ((i = 0), (th = th0 + h); i < ns; (i++), (th += h)) {
    p = GQR_POLYVAL(P, (*x) + k0);
    q = GQR_POLYVAL(Q, (*x) + k0);
    r = GQR_POLYVAL(R, (*x) + k0);
    dp = GQR_POLYDER(P, (*x) + k0);
    dq = GQR_POLYDER(Q, (*x) + k0);
    dr = GQR_POLYDER(R, (*x) + k0);
    k1 =
      -h * p * r / (r * sqrt(p * r) +
                    0.25 * (dr * p - dp * r + 2 * r * q) * sin(2 * th));
    *x += 0.5 * (k0 + k1);
    k1 = k0;
  }

  return;
}

static void gqr_taylor_coefficients(gdouble * U, gint n,
                                    gdouble u, gdouble du,
                                    gdouble p, gdouble q, gdouble r,
                                    gdouble dp, gdouble dq, gdouble dr,
                                    gdouble d2p, gdouble d2q, gdouble d2r)
{
  gint k;

  U[0] = u;
  U[1] = du;
  U[2] = -(q * U[1] + r * U[0]) / p;
  U[3] = -(dp * U[2] + dq * U[1] + q * U[2] + dr * U[0] + r * U[1]) / p;

  for (k = 2; k < n;
       (U[k + 2] =
        -((k * dp + q) * U[k + 1] +
          (0.5 * k * (k - 1) * d2p + k * dq + r) * U[k] +
          (0.5 * k * (k - 1) * d2q + k * dr) * U[k - 1] + 0.5 * k * (k -
                                                                     1) *
          d2r * U[k - 2]) / p), (k++));

  for (k = 0; k <= n; (U[k] *= GQR_INVERSE_FACTORIALS[k]), (k++));

  return;
}

static void gqr_newton_solve(gdouble * U, gint n, gdouble x0,
                             gdouble tol, gdouble * x, gdouble * df)
{
  gint j;
  gdouble f;

  for ((j = 0), (f = 1.0); (f > 0.0 ? f : -f) > tol && j < 16;
       (gqr_polyval(U, n, *x - x0, &f, df)), (*x -= f / (*df)), (j++));

  return;
}

/** 
 * Find the roots of the special function which solves the
 * differential equation \f$p(x)u''(x) + q(x)u'(x) + r(x)u(x) = 0\f$,
 * using the algorithm of Glaser, A., Liu, X. and Rokhlin, V., `A fast
 * algorithm for the calculation of the roots of special functions',
 * SIAM Journal on Scientific Computing, 29(4):1420--1438,
 * http://dx.doi.org/10.1137/06067016X.
 * 
 * @param P three element array containing the coefficients of
 * \f$p(x)\f$, lowest order first;
 * @param Q three element array containing the coefficients of
 * \f$q(x)\f$, lowest order first;
 * @param R three element array containing the coefficients of
 * \f$r(x)\f$, lowest order first;
 * @param x0 initial root;
 * @param du0 derivative of function at \a x0, \f$u'(x_{0})\f$;
 * @param N number of roots and derivatives to be returned;
 * @param x on return, the roots of \f$u(x)\f$;
 * @param du on return, \f$u'(x)\f$ at roots of \f$u(x)\f$.
 * 
 * @return 0 on success.
 */

gint gqr_function_roots(gdouble * P, gdouble * Q, gdouble * R,
                        gdouble x0, gdouble du0, gint N,
                        gdouble * x, gdouble * du)
{
  gdouble U[GQR_TAYLOR_LENGTH + 2], tol = 1e-18;
  gint i;

  g_return_val_if_fail(P != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(Q != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(R != NULL, GQR_NULL_PARAMETER);

  for ((x[0] = x[1] = x0), (du[0] = du0), (i = 1); i < N;
       (i++), (x[i] = x[i - 1])) {
    gqr_runge_kutta_solve(M_PI_2, -M_PI_2, 5, P, Q, R, &(x[i]));
    gqr_taylor_coefficients(U, GQR_TAYLOR_LENGTH, 0.0, du[i - 1],
                            GQR_POLYVAL(P, x[i - 1]), GQR_POLYVAL(Q,
                                                                  x[i -
                                                                    1]),
                            GQR_POLYVAL(R, x[i - 1]), GQR_POLYDER(P,
                                                                  x[i -
                                                                    1]),
                            GQR_POLYDER(Q, x[i - 1]), GQR_POLYDER(R,
                                                                  x[i -
                                                                    1]),
                            GQR_POLYDER2(P, x[i - 1]), GQR_POLYDER2(Q,
                                                                    x[i -
                                                                      1]),
                            GQR_POLYDER2(R, x[i - 1]));
    gqr_newton_solve(U, GQR_TAYLOR_LENGTH, x[i - 1], tol, &x[i], &du[i]);
  }

  return 0;
}

gint gqr_function_nextroot(gdouble * P, gdouble * Q, gdouble * R,
                           gdouble x0, gdouble u0, gdouble * x, gdouble * du)
{
  gdouble U[GQR_TAYLOR_LENGTH + 2], tol = 1e-18;

  g_return_val_if_fail(P != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(Q != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(R != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(x != NULL, GQR_NULL_PARAMETER);
  g_return_val_if_fail(du != NULL, GQR_NULL_PARAMETER);

  *x = x0;
  gqr_runge_kutta_solve(0.0, -M_PI_2, 5, P, Q, R, x);
  gqr_taylor_coefficients(U, GQR_TAYLOR_LENGTH, u0, 0.0,
                          GQR_POLYVAL(P, x0), GQR_POLYVAL(Q, x0),
                          GQR_POLYVAL(R, x0),
                          GQR_POLYDER(P, x0), GQR_POLYDER(Q, x0),
                          GQR_POLYDER(R, x0),
                          GQR_POLYDER2(P, x0), GQR_POLYDER2(Q, x0),
                          GQR_POLYDER2(R, x0));
  gqr_newton_solve(U, GQR_TAYLOR_LENGTH, x0, tol, x, du);

  return 0;
}

/**
 * @}
 * 
 */
