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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include "gqr.h"
#include "gqr-private.h"

/*
  The third argument is a pointer to match the GQRFunc typedef.
 */

gdouble legendre_func(gdouble t, gint d, gint * m)
{
  if (t < -1.0 || t > 1.0)
    g_error("%s: t=%lg out of range", __FUNCTION__, t);

  if (d == 0)
    return gsl_sf_legendre_Pl(*m, t);

  return gqr_legendre_dPdx(t, *m, d);

  if (d > *m)
    return 0.0;

  if (d == 1) {
    if (t == -1.0)
      return (gsl_sf_pow_int(-1, *m - d) * gsl_sf_pow_int(0.5, d) *
              gsl_sf_fact(*m + d) / gsl_sf_fact(d) / gsl_sf_fact(*m - d));

    if (t == 1.0)
      return (0.5 * (*m) * (*m + 1));

    return (t * gsl_sf_legendre_Pl(*m, t) -
            gsl_sf_legendre_Pl(*m - 1, t)) / (t * t - 1.0);
  }

  g_assert_not_reached();

  return 0;
}

/** 
 * Compute a quadrature rule for the integration of functions with
 * singularities and near-singularities of various orders.
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

gint grule_multi_nsingular(gint n, gint m, gdouble x, gdouble y,
                           gint nsing, gint * s, gdouble * xk, gdouble * wt)
{
  gdouble *P, *I;
  gdouble R;
  gint i, j, k, ns;
  gdouble *psi, *M;

  ns = nsing + 1;
  grule_legendre(n, xk, wt);
  P = (gdouble *) g_malloc(n * n * sizeof(gdouble));
  M = (gdouble *) g_malloc(MAX(ns * m, n) * sizeof(gdouble));
  psi = (gdouble *) g_malloc(ns * m * n * sizeof(gdouble));
  I = (gdouble *) g_malloc(2 * m * sizeof(gdouble));

  for (i = 0; i < ns * m * n; i++)
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

  for (k = 0; k < nsing; k++) {
    switch (s[k]) {
    case 0:
      gqr_quad_log(m, x, fabs(y), I);
      gqr_legendre_integrals(m, I, &(M[(k + 1) * m]));
      break;
    case 1:
      gqr_quad_1_r(m, x, fabs(y), I);
      gqr_legendre_integrals(m, I, &(M[(k + 1) * m]));
      break;
    case 2:
      gqr_quad_2_r(m, x, fabs(y), I);
      gqr_legendre_integrals(m, I, &(M[(k + 1) * m]));
      break;
    case 3:
      gqr_quad_3_r(m, x, fabs(y), I);
      gqr_legendre_integrals(m, I, &(M[(k + 1) * m]));
      break;
    default:
      g_error("%s: quadrature of near-singular functions with order %d "
              "singularity not yet implemented ", __FUNCTION__, s[k]);
      break;
    }
  }

  for (i = 0; i < n; i++) {
    R = sqrt((x - xk[i]) * (x - xk[i]) + y * y);
    for (j = 0; j < m; j++) {
      psi[i * ns * m + j + 0 * m] = P[j * n + i];
    }
  }

  for (k = 0; k < nsing; k++) {
    switch (s[k]) {
    case 0:
      for (i = 0; i < n; i++) {
        R = sqrt((x - xk[i]) * (x - xk[i]) + y * y);
        for (j = 0; j < m; j++) {
          psi[i * ns * m + j + (k + 1) * m] =
            psi[i * ns * m + j + 0 * m] * log(R);
        }
      }
      break;
    default:
      for (i = 0; i < n; i++) {
        R = sqrt((x - xk[i]) * (x - xk[i]) + y * y);
        for (j = 0; j < m; j++) {
          psi[i * ns * m + j + (k + 1) * m] =
            psi[i * ns * m + j + 0 * m] / gsl_pow_int(R, s[k]);
        }
      }
      break;
    }
  }

  gqr_lsqr_min_norm(psi, ns * m, n, M, MAX(ns * m, n));

  for (i = 0; i < n; i++)
    wt[i] = M[i];

  g_free(P);
  g_free(M);
  g_free(I);

  return 0;
}

gint grule_multi_singular(gint n, gint m, gdouble x,
                          gint nsing, gint * s, gdouble * xk, gdouble * wt)
{
  gdouble *fp2;
  gdouble *P, *Q, *Px, *dPx;
  gdouble *fp, pv;
  gdouble t;
  gqr_rule_t *g, *h;
  gint i, j, k, ss, ns;
  gint smax;
  gdouble *psi, *M;

  ns = nsing + 1;
  g = gqr_rule_alloc(n);
  gqr_rule_select(g, GQR_GAUSS_LEGENDRE, n, NULL);
  h = gqr_rule_alloc(2 * m);
  gqr_rule_select(h, GQR_GAUSS_LEGENDRE, 2 * m, NULL);
  k = gqr_rule_length(g);
  Px = (gdouble *) g_malloc((m + 2) * sizeof(gdouble));
  dPx = (gdouble *) g_malloc((m + 1) * sizeof(gdouble));
  P = (gdouble *) g_malloc((n + 1) * k * sizeof(gdouble));
  Q = (gdouble *) g_malloc(2 * (n + 2) * sizeof(gdouble));
  M = (gdouble *) g_malloc(MAX(ns * m, n) * sizeof(gdouble));
  psi = (gdouble *) g_malloc(ns * m * n * sizeof(gdouble));
  fp2 = (gdouble *) g_malloc(ns * n * sizeof(gdouble));

  smax = 0;
  for (i = 0; i < nsing; i++)
    smax = MAX(smax, s[i]);
  fp = (gdouble *) g_malloc((smax + 1) * sizeof(gdouble));

  for (i = 0; i < ns * m * n; i++)
    psi[i] = 0.0;
  i = 0;
  for (j = 0; j < k; j++) {
    P[i * k + j] = 1.0;
    P[(i + 1) * k + j] = gqr_rule_abscissa(g, j);
  }
  if ((x != -1.0) && (x != 1.0)) {
    Q[0] = 0.5 * log(fabs((x + 1) / (x - 1)));
    Q[1] = x * Q[0] - 1.0;
  }
  else {
    Q[0] = M_LN2;
    Q[1] = M_LN2 - 0.5;
  }

  for (i = 1; i < n; i++) {
    for (j = 0; j < k; j++) {
      P[(i + 1) * k + j] =
        ((gdouble) (2 * i + 1) * gqr_rule_abscissa(g, j) * P[i * k + j] -
         (gdouble) i * P[(i - 1) * k + j]) / (gdouble) (i + 1);
    }
    Q[i + 1] =
      ((gdouble) (2 * i + 1) * x * Q[i] -
       (gdouble) i * Q[i - 1]) / (gdouble) (i + 1);
  }
  i = n;
  Q[i + 1] =
    ((gdouble) (2 * i + 1) * x * Q[i] -
     (gdouble) i * Q[i - 1]) / (gdouble) (i + 1);

  Px[0] = 1.0;
  Px[1] = x;
  for (i = 1; i <= m; i++) {
    Px[i + 1] =
      ((gdouble) (2 * i + 1) * x * Px[i] -
       (gdouble) i * Px[i - 1]) / (gdouble) (i + 1);
  }

  dPx[0] = 0.0;
  dPx[1] = 1.0;
  if ((x == 1.0) || (x == -1.0))
    for (i = 1; i < m; i++)
      dPx[i + 1] =
        x * gsl_pow_int(x, i + 1) * 0.5 * (gdouble) ((i + 1) * (i + 2));
  else
    for (i = 1; i < m; i++)
      dPx[i + 1] = (i + 1) * (x * Px[i + 1] - Px[i]) / (x * x - 1.0);

  if ((x == 1.0) || (x == -1.0))
    fp2[0] = 2.0 * M_LN2 - 2.0;
  else
    fp2[0] = 2 * Q[1] + log(fabs(x * x - 1.0));

  if ((x == 1.0) || (x == -1.0))
    for (i = 1; i < n; i++)
      fp2[i] = gqr_finite_part_Pn_log(x, i);
  else
    for (i = 1; i < n; i++)
      fp2[i] = 2.0 / (2.0 * i + 1) * (Q[i + 1] - Q[i - 1]);

  pv = gqr_finite_part(-1, 1, x, 1.0);
  for (i = 1; i <= smax; i++)
    fp[i] = gqr_finite_part(-1, 1, x, (gdouble) i);

  /*fill the unweighted Legendre polynomial entries */
  for (i = 0; i < m; i++)
    M[i] = 0.0;
  M[0 * m + 0] = 2.0;

  for (i = 0; i < n; i++) {
    t = gqr_rule_abscissa(g, i);
    for (j = 0; j < m; j++)
      psi[i * ns * m + j + 0 * m] = P[j * k + i];
  }

  for (ss = 0; ss < nsing; ss++) {
    if (s[ss] == 0) {
      for (i = 0; i < m; i++)
        M[(ss + 1) * m + i] = fp2[i];
      for (i = 0; i < n; i++) {
        t = gqr_rule_abscissa(g, i);
        for (j = 0; j < m; j++)
          psi[i * ns * m + j + (ss + 1) * m] =
            psi[i * ns * m + j + 0 * m] * log(fabs(x - t));
      }
    }
    else {
      for (i = 0; i < m; i++)
        M[(ss + 1) * m + i] =
          gqr_finite_part_integral((gqr_func) legendre_func, &i, x,
                                   (gdouble) s[ss], -1, 1, h);

      for (i = 0; i < n; i++) {
        t = gqr_rule_abscissa(g, i);
        for (j = 0; j < m; j++)
          psi[i * ns * m + j + (ss + 1) * m] =
            psi[i * ns * m + j + 0 * m] / gsl_sf_pow_int(x - t, s[ss]);
      }
    }
  }

  gqr_lsqr_min_norm(psi, ns * m, n, M, MAX(ns * m, n));

  for (i = 0; i < n; i++) {
    xk[i] = gqr_rule_abscissa(g, i);
    wt[i] = M[i];
  }

  gqr_rule_free(g);
  gqr_rule_free(h);
  g_free(P);
  g_free(Q);
  g_free(Px);
  g_free(dPx);
  g_free(M);
  g_free(fp2);
  g_free(fp);

  return 0;
}
