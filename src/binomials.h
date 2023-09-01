#ifndef _BINOMIALS_H_INCLUDED_
#define _BINOMIALS_H_INCLUDED_

extern const gdouble _BINOMIALS[] ;
extern const gdouble GQR_FACTORIALS[] ;
#define gqr_binomial(_m,_k) (_BINOMIALS[(_m)*((_m)+1)/2+(_k)])

#define gqr_binomial_list(_m) &(_BINOMIALS[(_m)*((_m)+1)/2])

#define gqr_factorial(_m) (GQR_FACTORIALS[(_m)])

#endif
