* Brief notes on installation requirements

GQR is a library and executable for generating quadrature rules of
Gaussian type. It can be used inside other codes where quadrature
rules are needed, or to precompute rules for future use. 

* Installation

You will need the autotools suite, BLAS and LAPACK installed, as well
as the wrappers in the blaswrap library, available from

https://github.com/mjcarley/blaswrap

If the rank revealing QR code of Bischof and Quintana-Orti is
available, GQR will use it to generate generalized quadratures. It is
available from

https://dl.acm.org/doi/10.1145/290200.287638

Otherwise, a built-in function will be used. 

After downloading and extracting the source code, generate the
configuration files by

  . autogen.sh

To configure and install GQR,

  ./configure [options]
  make
  make install

The configure script takes a number of options and is influenced by
certain environment variables. Use

  /configure --help

to see what these are.

This will generate the library, an executable called gqr-rule, and the
documentation for using the library. To precompute rules, use
gqr-rule. For example:

  gqr-rule -T -N 7

generates the seven point Gauss-Chebyshev quadrature of the second
kind. To see what else gqr-rule can do, view the help using

  gqr-rule -h

