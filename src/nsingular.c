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

#include "gqr.h"
#include "gqr-private.h"

/** 
 * Compute a quadrature rule for the integration of functions with
 * singularities and near-singularities up to second order, where the
 * singularities are, for example, of the form
 * \f$1/((x-t)^{2}+y^{2})\f$, \f$y\neq0\f$.
 * 
 * @param n number of points in quadrature rule
 * @param m maximum order of polynomial to integrate
 * @param x position of singularity
 * @param y position of singularity
 * @param xk abscissae of quadrature rule
 * @param wt weights of quadrature rule
 * 
 * @return 0 on success
 */

gint grule_near_singular(gint n, gint m, gdouble x, gdouble y,
                         gdouble * xk, gdouble * wt)
{
  gdouble *P, *I;
  gdouble R;
  gint i, j;
  gdouble *psi, *M;

  grule_legendre(n, xk, wt);
  P = (gdouble *) g_malloc(n * n * sizeof(gdouble));
  M = (gdouble *) g_malloc(MAX(4 * m, n) * sizeof(gdouble));
  psi = (gdouble *) g_malloc(4 * m * n * sizeof(gdouble));
  I = (gdouble *) g_malloc(2 * m * sizeof(gdouble));

  for (i = 0; i < 4 * m * n; i++)
    psi[i] = 0.0;
  i = 0;
  for (j = 0; j < n; j++) {
    P[i * n + j] = 1.0;
    P[(i + 1) * n + j] = xk[j];
  }

  for (i = 1; i < m - 1; i++) {
    for (j = 0; j < n; j++) {
      P[(i + 1) * n + j] =
        ((gdouble) (2 * i + 1) * xk[j] * P[i * n + j] -
         (gdouble) i * P[(i - 1) * n + j]) / (gdouble) (i + 1);
    }
  }

  M[0 * m + 0] = 2.0;
  for (i = 1; i < m; i++)
    M[0 * m + i] = 0.0;
  gqr_quad_log(m, x, fabs(y), I);
  gqr_legendre_integrals(m, I, &(M[1 * m]));
  gqr_quad_1_r(m, x, fabs(y), I);
  gqr_legendre_integrals(m, I, &(M[2 * m]));
  gqr_quad_2_r(m, x, fabs(y), I);
  gqr_legendre_integrals(m, I, &(M[3 * m]));

  for (i = 0; i < n; i++) {
    R = sqrt((x - xk[i]) * (x - xk[i]) + y * y);
    for (j = 0; j < m; j++) {
      psi[i * 4 * m + j + 0 * m] = P[j * n + i];
      psi[i * 4 * m + j + 1 * m] = psi[i * 4 * m + j + 0 * m] * log(R);
      psi[i * 4 * m + j + 2 * m] = psi[i * 4 * m + j + 0 * m] / R;
      psi[i * 4 * m + j + 3 * m] = psi[i * 4 * m + j + 2 * m] / R;
    }
  }

  gqr_lsqr_min_norm(psi, 4 * m, n, M, MAX(4 * m, n));

  for (i = 0; i < n; i++)
    wt[i] = M[i];

  g_free(P);
  g_free(M);
  g_free(I);

  return 0;
}
