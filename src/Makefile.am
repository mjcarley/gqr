## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"GQR\"

bin_SCRIPTS=gqr-config

BUILT_SOURCES= \
	gqr-config

gqr-config: gqr-config.in

lib_LTLIBRARIES = libgqr.la

libgqr_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libgqr_la_SOURCES = \
	ghermite.c \
        grule.c  \
	fpart.c \
	glegendre.c \
	korsunsky.c  \
	smith.c \
	gchebyshev.c  \
	kolrokh.c  \
	lsqr.c \
	msingular.c \
	nsingular.c \
	integrals.c \
	roots.c \
	gqr-logging.c \
	gqr-private.h

include_HEADERS = \
	gqr.h 

CLEANFILES = $(BUILT_SOURCES)
