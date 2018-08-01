/* Copyright (C) 2007, 2010 by  Michael Carley */

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

/**
 * @defgroup legendre  Integrals of Legendre polynomials
 * @{
 *
 * The functions in this section are not necessary for most users but
 * might be helpful under certain circumstances. They allow the
 * evaluation of integrals of Legendre polynomials with arbitrary
 * weight functions.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H */

#include "gqr.h"
#include "gqr-private.h"

gint gqr_quad_3_r(gint N, gdouble x, gdouble y, gdouble * I)
{
  gdouble Rp1, Rm1, R2, Ir;
  gint i, sgn;

  g_assert(y != 0.0);
  Rm1 = sqrt((x + 1.0) * (x + 1.0) + y * y);
  Rp1 = sqrt((x - 1.0) * (x - 1.0) + y * y);
  R2 = x * x + y * y;
  if (y == 0.0) {
    I[0] = 0.5 / (x - 1.0) / (x - 1.0) - 0.5 / (x + 1.0) / (x + 1.0);
    I[1] = 0.5 / (x - 1.0) / (x - 1.0) + 1.5 / (x + 1.0) / (x + 1.0);
    I[2] = log((x + 1.0) / (x - 1.0)) + 0.5 / (x - 1.0) / (x - 1.0) +
      3.5 / (x + 1.0) / (x + 1.0);
    g_assert(!isnan(I[0]));
    g_assert(!isnan(I[1]));
    g_assert(!isnan(I[2]));
  }
  else {
    I[0] = ((1.0 - x) / Rp1 + (1.0 + x) / Rm1) / y / y;
    I[1] = ((x - x * x - y * y) / Rp1 + (x + x * x + y * y) / Rm1) / y / y;
    Ir = log(2 * Rp1 - 2 * (x - 1)) - log(2 * Rm1 - 2 * (x + 1));
    I[2] = Ir - (y * y * (x + 1.0) + x * x * (x - 1.0)) / Rp1 / y / y +
      (y * y * (x - 1.0) + x * x * (x + 1.0)) / Rm1 / y / y;
  }

  for ((i = 3), (sgn = 1); i <= N; (i++), (sgn = -sgn)) {
    I[i] = (1.0 / Rp1 - (gdouble) sgn / Rm1) / (gdouble) (i - 2);
    I[i] +=
      (2 * i - 3) * x / (i - 2) * I[i - 1] - (i - 1) * (x * x + y * y) / (i -
                                                                          2) *
      I[i - 2];
    g_assert(!isnan(I[i]));
  }

  return 0;
}

gint gqr_quad_2_r(gint N, gdouble x, gdouble y, gdouble * I)
{
  gdouble Rp1, Rm1, R2;
  gint i, sgn;

/*   g_assert(y != 0.0) ; */
  Rm1 = (x + 1.0) * (x + 1.0) + y * y;
  Rp1 = (x - 1.0) * (x - 1.0) + y * y;
  R2 = x * x + y * y;
  if (y == 0.0) {
    I[0] = 2.0 / (x * x - 1.0);
    I[1] = 2.0 / (x * x - 1.0) + log((x - 1.0) / (x + 1.0));
  }
  else {
    I[0] = (atan2(1.0 - x, y) - atan2(-1.0 - x, y)) / y;
    I[1] = 0.5 * (log(Rp1) - log(Rm1)) + x * I[0];
  }

  for ((i = 2), (sgn = 1); i <= N; (i++), (sgn = -sgn)) {
    I[i] = 2.0 * x * I[i - 1] - R2 * I[i - 2];
    I[i] += (1 + sgn) / (gdouble) (i - 1);
    if (isnan(I[i]))
      g_error("%s: NaN error", __FUNCTION__);
  }

  return 0;
}

gint gqr_quad_1_r(gint N, gdouble x, gdouble y, gdouble * I)
{
  gdouble Rm1, Rp1, R2;
  gint i, sgn;

/*   g_assert(y != 0.0) ; */
  Rp1 = sqrt((x - 1.0) * (x - 1.0) + y * y);
  Rm1 = sqrt((x + 1.0) * (x + 1.0) + y * y);
  R2 = x * x + y * y;
  I[0] = log(2 * Rp1 - 2 * (x - 1)) - log(2 * Rm1 - 2 * (x + 1));
  I[1] = Rp1 - Rm1 + x * I[0];

  for ((i = 2), (sgn = 1); i <= N; (i++), (sgn = -sgn)) {
    I[i] = (x * (gdouble) (2 * i - 1) * I[i - 1] -
            R2 * (gdouble) (i - 1) * I[i - 2]) / (gdouble) i;
    I[i] += (Rp1 + sgn * Rm1) / (gdouble) i;
    g_assert(!isnan(I[i]));
  }

  return 0;
}

gint gqr_quad_log(gint N, gdouble x, gdouble y, gdouble * I)
{
  gdouble R;
  gdouble Lm1, Lp1;
  gdouble th, C;
  gdouble *I2;
  gint n, sgn;

/*   g_assert(y != 0.0) ; */

  I2 = (gdouble *) g_malloc((N + 3) * sizeof(gdouble));
  gqr_quad_2_r(N + 2, x, y, I2);
  R = sqrt(x * x + y * y);
  th = atan2(y, x);
  C = R * cos(th);
  Lm1 = log((x + 1) * (x + 1) + y * y);
  Lp1 = log((x - 1) * (x - 1) + y * y);

  for ((n = 0), (sgn = 1); n <= N; (n++), (sgn = -sgn)) {
    I[n] = 0.5 * (Lp1 + sgn * Lm1) - I2[n + 2] + C * I2[n + 1];
    I[n] /= (n + 1);
    g_assert(!isnan(I[n]));
  }

  g_free(I2);

  return 0;
}

/** 
 * Process a list of integrals of powers on a weight function to
 * generate the integrals of the corresponding Legendre polynomials.
 * 
 * @param N number of integrals in list
 * @param I integrals \f$\int w(t)t^{n}\,\mathrm{d}t\f$ for 
 * \f$n=0,\ldots,N-1\f$
 * @param P on output integrals \f$\int w(t)P_{n}(t)\,\mathrm{d}t\f$, 
 * for \f$n=0,\ldots,N-1\f$
 * 
 * @return 0 on success.
 */

gint gqr_legendre_integrals(gint N, gdouble * I, gdouble * P)
{
  gint i, j;
  gdouble Cn, Dn;

  if (N > 32)
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_WARNING,
          "%s: algorithm suffers from overflow for N>32", __FUNCTION__);

  P[0] = I[0];
  P[1] = I[1];
  for (i = 2; i < N; i += 2) {
    Cn = 1.0 / (gdouble) (i + 1);
    Dn = 1.0 / (gdouble) (i + 3);
    P[i] = I[i] - Cn * P[0];
    P[i + 1] = I[i + 1] - 3.0 * Dn * P[1];
    Cn *= (gdouble) (i) / (gdouble) (i + 3);
    Dn *= (gdouble) (i) / (gdouble) (i + 5);
    for (j = 1; j < i / 2; j++) {
      P[i] -= Cn * (4 * j + 1) * P[2 * j];
      Cn *= (gdouble) (i - 2 * j) / (gdouble) (i + 2 * j + 3);
      P[i + 1] -= Dn * (4 * j + 3) * P[2 * j + 1];
      Dn *= (gdouble) (i - 2 * j) / (gdouble) (i + 2 * j + 5);
    }
    P[i] /= Cn * (i * 2 + 1);
    P[i + 1] /= Dn * (i * 2 + 3);
    g_assert(!isnan(P[i]));
    g_assert(!isnan(P[i + 1]));
  }

  return 0;
}

/**
 * @}
 * 
 */
