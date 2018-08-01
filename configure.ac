AC_INIT([gqr], 0.1.3)
AC_C_CONST
AC_C_INLINE
AC_PROG_LIBTOOL
AC_CONFIG_SRCDIR([src/gqr.h])
AC_CONFIG_MACRO_DIR([m4])

# Making releases:
#   GQR_MICRO_VERSION += 1;
#   GQR_INTERFACE_AGE += 1;
#   GQR_BINARY_AGE += 1;
# if any functions have been added, set GQR_INTERFACE_AGE to 0.
# if backwards compatibility has been broken,
# set GQR_BINARY_AGE and GQR_INTERFACE_AGE to 0.
#
GQR_MAJOR_VERSION=0
GQR_MINOR_VERSION=1
GQR_MICRO_VERSION=3
GQR_INTERFACE_AGE=1
GQR_BINARY_AGE=1
GQR_VERSION=$GQR_MAJOR_VERSION.$GQR_MINOR_VERSION.$GQR_MICRO_VERSION

AC_SUBST(GQR_MAJOR_VERSION)
AC_SUBST(GQR_MINOR_VERSION)
AC_SUBST(GQR_MICRO_VERSION)
AC_SUBST(GQR_VERSION)

# libtool versioning
LT_RELEASE=$GQR_MAJOR_VERSION.$GQR_MINOR_VERSION
LT_CURRENT=`expr $GQR_MICRO_VERSION - $GQR_INTERFACE_AGE`
LT_REVISION=$GQR_INTERFACE_AGE
LT_AGE=`expr $GQR_BINARY_AGE - $GQR_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$GQR_VERSION
PACKAGE=gqr

AM_INIT_AUTOMAKE([foreign])

# Specify a configuration file
AC_CONFIG_HEADERS([config.h])

AC_DEFINE_UNQUOTED(GQR_MAJOR_VERSION, $GQR_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(GQR_MINOR_VERSION, $GQR_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(GQR_MICRO_VERSION, $GQR_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(GQR_INTERFACE_AGE, $GQR_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(GQR_BINARY_AGE, $GQR_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
AM_MAINTAINER_MODE

AC_PROG_CC

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

AC_CHECK_LIB(m, cos)

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS
ACX_BLAS
ACX_LAPACK


AC_CHECK_PROG(GSL, gsl-config, yes, no)
if test x$GSL = xyes ; then
   gsl_cflags=`gsl-config --cflags`
   gsl_libs=`gsl-config --libs`
else
  echo \#################################################################
  echo
  echo \# GNU Scientific Library not available or not properly installed.
  echo \# You can get the most recent version from
  echo \#      http://www.gnu.org/software/gsl/
  echo
  echo \#################################################################
  exit
fi

PKG_CHECK_MODULES([GLIB], [glib-2.0])

AC_SUBST(GLIB_CFLAGS)
AC_SUBST(GLIB_LIBS)
AC_SUBST(gsl_cflags)
AC_SUBST(gsl_libs)

CFLAGS="$CFLAGS $GLIB_CFLAGS $gsl_cflags"
LDFLAGS="$LDFLAGS $GLIB_LIBS $gsl_libs"

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

dnl header file checks
#AC_CHECK_HEADERS([f2c.h], [AC_DEFINE(HAVE_F2C_H)])

#AC_CHECK_LIB(m, sqrt)
#AC_CHECK_LIB(F77, s_cmp)
#AC_CHECK_LIB(blas, f2c_snrm2)
#AC_CHECK_LIB(lapack, dgelsd_,,[AC_MSG_NOTICE(
#             [** LAPACK not found, some features will not be available])])

LIBS="$glib_libs $gsl_libs $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_CONFIG_FILES([
Makefile
src/Makefile
src/gqr-config
tests/Makefile
doc/Makefile
doc/html/Makefile
])
AC_OUTPUT
