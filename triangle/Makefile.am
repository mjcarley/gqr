AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

BUILT_SOURCES = \
	radial-rules.c \
	angular-rules.c \
	legendre-rules.c

CLEANFILES = \
	radial-rules.c \
	angular-rules.c \
	legendre-rules.c

bin_PROGRAMS = triangle-tests

radial-rules.c:
	./radial-rules

angular-rules.c:
	./angular-rules

legendre-rules.c:
	./legendre-rules

triangle_tests_SOURCES = \
	triangle-tests.c \
	wandzura.c \
	radial-rules.c \
	angular-rules.c \
	legendre-rules.c \
	shape.c
