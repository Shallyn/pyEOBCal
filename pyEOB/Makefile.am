# Makefile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
pkgpyexecdir = ${HOME}/usr/local

pkgpyexec_LTLIBRARIES = PyUtils.la

PyUtils_la_SOURCES = pyFunctions.c

PyUtils_la_LIBADD = ../csrc/libEOB.la
PyUtils_la_LDFLAGS = -module -avoid-version

all: $(pkgpyexec_LTLIBRARIES:.la=.so)

.la.so:
	ln -s .libs/$@ $@ || true
