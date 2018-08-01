/* Copyright (C) 2007, 2008 by  Michael Carley */

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
 * @defgroup finite Evaluation of finite part integrals
 *
 * @{
 *
 */

#include <stdio.h>
#include <math.h>

#include <glib.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "gqr.h"

/**
 * Compute the finite part integral
 * \f$\mathrm{FP}\int_{a}^{b}1/(y-x)^{\gamma}\,\mathrm{d}x\f$ using
 * the method of
 * Brandao, M. P., 1987, `Improper integrals in theoretical
 * aerodynamics: The problem revisited', AIAA Journal,
 * 25(9):1258--1260.
 *
 * @param a lower limit of integration;
 * @param b upper limit of integration;
 * @param y location of singularity;
 * @param g order of singularity, \f$\gamma\f$.
 *
 * @return finite part of
 *  \f$\int_{a}^{b}1/(y-x)^{\gamma}\,\mathrm{d}x\f$.
 */

gdouble gqr_finite_part(gdouble a, gdouble b, gdouble y, gdouble g)
{
  gint n;
  gdouble alpha, fp, bmy, yma;

  if (g >= 1.0) {
    n = (gint) floor(g - 1.0);
    alpha = g - n - 1.0;
  }
  else {
    n = 0;
    alpha = g;
  }

  g_assert(n >= 0);
  g_assert(a < b);

  if ((y > a) && (y < b)) {
    if (alpha != 0.0) {
      bmy = pow(b - y, n + alpha);
      yma = pow(y - a, n + alpha);
      fp = -(bmy + yma * cos((alpha + n + 1) * M_PI)) /
        (alpha + n) / bmy / yma;
    }
    else {
      if (n == 0) {
        fp = log((y - a) / (b - y));
      }
      else {
        bmy = gsl_pow_int(b - y, n);
        yma = gsl_pow_int(y - a, n);
        fp = -(bmy + yma * pow(-1, n + 1)) / n / bmy / yma;
      }
    }
    return fp;
  }

  if (y == a) {
    if (alpha == 0.0) {
      if (n == 0)
        fp = -log(b - a);
      else {
        fp = pow(-1.0, n) / (gdouble) n / gsl_pow_int(b - a, n);
      }
    }
    else {
      fp = -cos((alpha + n + 1) * M_PI) / (alpha + n) / pow(b - a, alpha + n);
    }
    return fp;
  }

  if (y == b) {
    if (alpha == 0.0) {
      if (n == 0)
        fp = log(b - a);        /*unchecked */
      else {
        fp = -1.0 / (gdouble) n / gsl_pow_int(b - a, n);
      }
    }
    else {
      fp = -1.0 / (alpha + n) / gsl_pow_int(b - a, alpha + n);  /*unchecked */
    }
    return fp;
  }

  if ((alpha != 0.0) || (n > 0)) {
    fp =
      1.0 / (alpha + n) * (1.0 / pow(y - b, alpha + n) -
                           1.0 / pow(y - a, alpha + n));
  }
  else {
    fp = log((y - a) / (y - b));
  }

  return fp;
}

/**
 * Compute the finite part integral
 * \f$FP\int_{a}^{b}f(x)/(y-x)^{\gamma}\,\mathrm{d}x\f$  using a
 * Gaussian quadrature for the regular integral. The method is that of
 * Brandao, M. P., 1987, `Improper integrals in theoretical
 * aerodynamics: The problem revisited', AIAA Journal,
 * 25(9):1258--1260.
 *
 * @param f a gqr_func which returns the \f$j\f$th derivative of the
 * function;
 * @param data data to pass to f;
 * @param y position of the singularity;
 * @param gm strength (exponent) of the singularity;
 * @param a lower limit of integration;
 * @param b upper limit of integration;
 * @param g Gaussian quadrature rule for the regular integral.
 *
 * @return the finite part of
 * \f$\int_{a}^{b}f(x)/(y-x)^{\gamma}\,\mathrm{d}x\f$.
 */

gdouble gqr_finite_part_integral(gqr_func f, gpointer data,
                                 gdouble y, gdouble gm,
                                 gdouble a, gdouble b, gqr_rule_t * g)
{
  gdouble xbar, dx, alpha, fp, x, df, cft, *fj;
  gint i, j, n, m;

  gqr_rule_scale(g, a, b, &xbar, &dx);

  n = (gint) floor(gm - 1.0);
  alpha = gm - n - 1.0;

  g_assert(alpha == 0.0);
  if (alpha == 0.0)
    m = n;
  else
    m = n + 1;

  fj = (gdouble *) g_malloc((m + 1) * sizeof(gdouble));
  for ((j = 0), (cft = -1); j <= m; j++) {
    fj[j] = cft * f(y, j, data);
    cft *= -1.0 / (gdouble) (j + 1);
  }

  fp = 0.0;
  for (i = 0; i < gqr_rule_length(g); i++) {
    x = xbar + dx * gqr_rule_abscissa(g, i);
    df = f(x, 0, data);
    for (j = 0; j <= m; j++)
      df += fj[j] * pow(y - x, j);
    fp += df / pow(y - x, gm) * gqr_rule_weight(g, i) * dx;
  }

  for (j = 0; j <= m; j++)
    fp -= gqr_finite_part(a, b, y, alpha + n + 1 - j) * fj[j];

  g_free(fj);

  return fp;
}

/**
 * Return the finite part of
 * \f$\int_{-1}^{1}(1-t)^n/(1+t)^m\,\mathrm{d}t\f$.
 *
 * @param n exponent of \f$(1-t)\f$;
 * @param m exponent of \f$(1+t)\f$;
 *
 * @return finite part of integral.
 */

gdouble gqr_finite_part_1mt_n(gint n, gint m)
{
  gdouble fp, cft;
  gint q;

  switch (m) {
  case 1:
    if (n == 0)
      return M_LN2;
    for ((fp = M_LN2), (q = 1); q <= n; q++)
      fp -= 1.0 / (gdouble) q;
    fp *= gsl_pow_int(2.0, n);
    return fp;
    break;
  case 2:
    if (n == 0)
      return 0.5;
    for ((fp = -n * M_LN2 - 1.0), (q = 2), (cft = -1.0); q <= n; q++) {
      fp -= cft / (gdouble) (q - 1) * gsl_sf_choose(n, q);
      cft *= -1;
    }
    fp *= gsl_pow_int(2.0, n - 1);
    return fp;
    break;
  default:
    g_assert_not_reached();
    break;
  }

  return 0.0;
}

/**
 * Compute the \f$m\f$th derivative of the \f$n\f$th Legendre
 * polynomial \f$P_{n}(x)\f$ by direct evaluation (i.e. without using
 * the recursion relation).
 *
 * @param x \f$x\f$;
 * @param n order of Legendre polynomial;
 * @param m order of derivative.
 *
 * @return \f$P_{n}^{(m)}(x)\f$.
 */

gdouble gqr_legendre_dPdx(gdouble x, gint n, gint m)
{
  gdouble P, x1, x2, t1, t2, cft;
  gint k;

  if (m > n)
    return 0.0;

  x1 = 1.0 - x;
  x2 = 1.0 + x;
  t1 = gsl_pow_int(-1.0, m);
  t2 = gsl_pow_int(-1.0, n);
  cft = gsl_pow_int(-1.0, m) / gsl_sf_fact(m) * gsl_sf_fact(n + m) *
    gsl_pow_int(0.5, m + 1) / gsl_sf_fact(n - m);
  for ((P = 0.0), (k = m); k <= n; k++) {
    P += (t1 + t2) * cft;
    t1 *= x1;
    t2 *= x2;
    cft *=
      -0.5 / (gdouble) (k + 1) * (gdouble) (n + k + 1) / (gdouble) (k - m +
                                                                    1) * (n -
                                                                          k);
  }

  return P;
}

static gdouble quad_1pt(gint m)
{
  if (m == -1)
    return M_LN2;
  if (m == 0)
    return 2.0;
  return gsl_pow_int(2, m + 1) / (gdouble) (m + 1);
}

/**
 * Compute the finite part of
 * \f$\int_{-1}^{1}P_{n}(x)/(1+t)^{m}\,\mathrm{d}x\f$.
 *
 * @param n order of Legendre polynomial;
 * @param m order of singularity.
 *
 * @return \f$\int_{-1}^{1}P_{n}(x)/(1+t)^{m}\,\mathrm{d}x\f$.
 */

gdouble gqr_finite_part_Pn1(gint n, gint m)
{
  gdouble fp, cft;
  gint k;

  cft = 0.5;
  for ((fp = 0.0), (k = 0); k <= n; k++) {
    fp += cft * (gqr_finite_part_1mt_n(k, m) +
                 gsl_pow_int(-1, n) * quad_1pt(k - m));
    cft *=
      -0.5 / (gdouble) (k + 1) * (gdouble) (n + k + 1) / (gdouble) (k +
                                                                    1) * (n -
                                                                          k);
  }

  return fp;
}

gdouble gqr_finite_part_Pn_log(gdouble x, gint n)
{
  gint m, q, k;
  gdouble cft, fp, tmp;

  fp = 0;
  m = (gint) floor(0.5 * (gdouble) n);

  if (fabs(x) == 1.0) {
    cft = gsl_sf_lnfact(n);
    cft = gsl_sf_lnfact(2 * n) - 2.0 * cft;
    cft = exp(cft);
    if (n == 2 * m) {
      for (k = 0; k <= m; k++) {
        for ((q = 0), (tmp = M_LN2); q <= m - k; q++)
          tmp -= 1.0 / (gdouble) (2 * q + 1);
        fp += cft * tmp / (gdouble) (2 * m - 2 * k + 1);
        cft *=
          -1.0 / (2 * n - 2 * k) / (2 * n - 2 * k - 1) / (k + 1) * (n -
                                                                    k) * (n -
                                                                          2 *
                                                                          k) *
          (n - 2 * k - 1);
      }
      fp *= 2.0 * gsl_pow_int(0.5, n);
    }
    else {
      for (k = 0; k <= m; k++) {
        for ((q = 0), (tmp = 0.0); q <= m - k; q++)
          tmp -= 1.0 / (gdouble) (2 * q + 1);
        fp += cft * tmp / (gdouble) (2 * m - 2 * k + 2);
        cft *=
          -1.0 / (2 * n - 2 * k) / (2 * n - 2 * k - 1) / (k + 1) * (n -
                                                                    k) * (n -
                                                                          2 *
                                                                          k) *
          (n - 2 * k - 1);
      }
      fp *= 2.0 * gsl_pow_int(0.5, n);
      if (x == -1)
        fp = -fp;
    }
    return fp;
  }

  return fp;
}

/**
 * Return the value of
 * \f$\int_{-1}^{1}t^{n}\log|x-t|\,\mathrm{d}t\f$.
 *
 * @@param n power of \f$t\f$;
 * @@param x position of singularity.
 *
 * @@return value of integral.
 */

gdouble gqr_finite_part_tn_log(gdouble x, gint n)
{
  gint k, i;
  gdouble xn, cft, pv;

  pv = 0.0;
  if (x == -1.0)
    pv = (1 - gsl_pow_int(-1, n + 1)) * M_LN2;
  if (x == 1.0)
    pv = (1 + gsl_pow_int(-1, n)) * M_LN2;

  if ((x > -1.0) && (x < 1.0)) {
    xn = gsl_pow_int(x, n + 1);
    pv = (1 - xn) * log(1 - x) + (gsl_pow_int(-1, n) + xn) * log(1 + x);
  }

  if (x < -1.0) {
    xn = gsl_pow_int(x, n + 1);
    pv = (1 - xn) * log(1.0 - x) + (gsl_pow_int(-1, n) + xn) * log(1 + x);
  }

  if (x > 1.0) {
    xn = gsl_pow_int(x, n + 1);
    pv = (1 - xn) * log(x - 1.0) + (gsl_pow_int(-1, n) + xn) * log(1 + x);
  }

  for ((cft = 1.0), (i = gsl_pow_int(-1, n)), (k = 1); k <= n + 1; k++) {
    pv -= (1 + i) * cft / (gdouble) (n - k + 2);
    cft *= x;
    i = -i;
  }
  pv /= (gdouble) (n + 1);

  return pv;
}

/**
 * @}
 *
 */
