/* Copyright (C) 2007, 2013 by  Michael Carley */

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
#include <glib.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H */

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#include "gqr.h"
#include "gqr-private.h"

static inline gdouble grule_lpoly(gint n, gdouble x)
{
  gdouble Pnp1, Pn, Pnm1;
  gint i;

  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;

  Pnm1 = 1.0;
  Pn = x;
  for (i = 1; i < n; i++) {
    Pnp1 = ((2 * i + 1) * x * Pn - i * Pnm1) / (gdouble) (i + 1);
    Pnm1 = Pn;
    Pn = Pnp1;
  }

  return Pn;
}

static inline gdouble grule_diff_lpoly(gint n, gdouble x)
{
  gdouble Pnp1, Pn, Pnm1;
  gdouble dPnp1, dPn, dPnm1;
  gint i;

  if (n == 0)
    return 0.0;
  if (n == 1)
    return 1.0;

  Pnm1 = 1.0;
  Pn = x;
  dPnm1 = 0.0;
  dPn = 1.0;
  for (i = 1; i < n; i++) {
    Pnp1 = ((2 * i + 1) * x * Pn - i * Pnm1) / (gdouble) (i + 1);
    dPnp1 =
      ((2 * i + 1) * x * dPn - i * dPnm1 + (2 * i + 1) * Pn) / (gdouble) (i +
                                                                          1);
    Pnm1 = Pn;
    Pn = Pnp1;
    dPnm1 = dPn;
    dPn = dPnp1;
  }

  return dPn;
}

gint grule_legendre(gint n, gdouble * x, gdouble * w)
{
  gdouble p[] = { 1, 0, -1 }
  , q[] = {
  0, -2, 0}
  , r[] = {
  0, 0, 0};
  gdouble x0, du;
  gint i;
  static GHashTable *cache = NULL;
  gdouble *xwc;

  g_assert(n > 0);

  if (cache == NULL)
    cache = g_hash_table_new(NULL, NULL);

  if ((xwc = g_hash_table_lookup(cache, GINT_TO_POINTER(n))) != NULL) {
    g_debug("%s: recovering cached rule for n=%d", __FUNCTION__, n);
    memcpy(x, xwc, n * sizeof(gdouble));
    memcpy(w, &(xwc[n]), n * sizeof(gdouble));
    return 0;
  }

  r[0] = n * (n + 1);
  if (2 * (n / 2) != n) {
    du = grule_diff_lpoly(n, 0.0);
    gqr_function_roots(p, q, r, 0, du, n / 2 + 1, &(x[n / 2]), &(w[n / 2]));
    for (i = n; i > n / 2; i--) {
      w[n - i - 1] = w[i] = 2.0 / (1.0 - x[i] * x[i]) / (w[i] * w[i]);
      x[n - i - 1] = -x[i];
    }
    w[n / 2] = 2.0 / (w[n / 2] * w[n / 2]);
    return 0;
  }

  gqr_function_nextroot(p, q, r, 0.0, grule_lpoly(n, 0.0), &x0, &du);
  gqr_function_roots(p, q, r, x0, du, n / 2, &(x[n / 2]), &(w[n / 2]));

  for (i = n - 1; i >= n / 2; i--) {
    w[n - i - 1] = w[i] = 2.0 / (1.0 - x[i] * x[i]) / w[i] / w[i];
    x[n - i - 1] = -x[i];
  }

  g_debug("%s: caching rule for n=%d", __FUNCTION__, n);
  xwc = (gdouble *) g_malloc(2 * n * sizeof(gdouble));
  memcpy(xwc, x, n * sizeof(gdouble));
  memcpy(&(xwc[n]), w, n * sizeof(gdouble));

  g_hash_table_insert(cache, GINT_TO_POINTER(n), xwc);

  return 0;
}
