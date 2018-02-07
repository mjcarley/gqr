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

#include <glib.h>

#include "gqr.h"
#include "gqr-private.h"

gint grule_logarithmic_smith(gint n, gdouble *x, gdouble *w)

{
  switch (n) {
  default: g_assert_not_reached() ; break ;
  case 16:
    x[0] = 2.3700206864932238e-03 ; x[1] = 3.0080510979786590e-02 ;
    x[2] =  1.1241533899502559e-01 ; x[3] = 2.3245061212753959e-01 ;
    x[4] = 3.2336545369335223e-01 ; x[5] = 3.6693498940855235e-01 ;
    x[6] =  4.3534313165564342e-01 ; x[7] = 4.9035421190908879e-01 ;
    x[8] =  5.0964578809091121e-01 ; x[9] = 5.6465686834435658e-01 ;
    x[10] =  6.3306501059144771e-01 ; x[11] = 6.7663454630664777e-01 ;
    x[12] =  7.6754938787246041e-01 ; x[13] = 8.8758466100497446e-01 ;
    x[14] =  9.6991948902021341e-01 ; x[15] = 9.9762997931350683e-01 ;
    w[0] = 8.7872775457998945e-03 ; w[1] = 5.2618273180791153e-02 ; 
    w[2] = 1.0931182786191369e-01 ; w[3] = 1.1851054984142320e-01 ; 
    w[4] = 5.4055166914867690e-02 ; w[5] = 5.5428789868546179e-02 ; 
    w[6] = 7.0833300620546866e-02 ; w[7] = 3.0454814166111351e-02 ; 
    w[8] = 3.0454814166111351e-02 ; w[9] = 7.0833300620546866e-02 ; 
    w[10] = 5.5428789868546179e-02 ; w[11] = 5.4055166914867690e-02 ; 
    w[12] = 1.1851054984142320e-01 ; w[13] = 1.0931182786191369e-01 ; 
    w[14] = 5.2618273180791153e-02 ; w[15] = 8.7872775457998945e-03 ;
    break ;
  case 12:
    x[0] = 3.6579100979781121e-03 ; x[1] = 4.3699309230310013e-02 ; 
    x[2] = 1.4890582891382192e-01 ; x[3] = 2.7373190858291296e-01 ; 
    x[4] = 3.4203076696372625e-01 ; x[5] = 4.3460004680260811e-01 ; 
    x[6] = 5.6539995319739189e-01 ; x[7] = 6.5796923303627375e-01 ; 
    x[8] = 7.2626809141708704e-01 ; x[9] = 8.5109417108617813e-01 ; 
    x[10] = 9.5630069076969004e-01 ; x[11] = 9.9634208990202189e-01 ;
    w[0] = 1.3406672201072389e-02 ; w[1] = 7.3037063871955160e-02 ;
    w[2] = 1.2907331563757191e-01 ; w[3] = 1.0344140876084500e-01 ;
    w[4] = 5.6876165283215312e-02 ; w[5] = 1.2416537424534030e-01 ;
    w[6] = 1.2416537424534030e-01 ; w[7] = 5.6876165283215312e-02 ;
    w[8] = 1.0344140876084500e-01 ; w[9] = 1.2907331563757191e-01 ;
    w[10] = 7.3037063871955160e-02 ; w[11] = 1.3406672201072389e-02 ;
    break ;
  }

  return 0 ;
}
