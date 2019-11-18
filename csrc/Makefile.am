# Makefile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
#pkgpyexecdir = ${HOME}/usr/local

#pkgpyexec_LTLIBRARIES = PyEOBCal.la
lib_LTLIBRARIES = libEOB.la
bin_PROGRAMS = MAIN

libEOB_la_SOURCES = \
dyEvolution.c dyFactorizedWaveform.c \
dyHamiltonian.c dyNewtonianMultipole.c \
dyIntegrator.c dyInitialCondition.c \
myOptParser.c dyHcapNumericalDerivative.c\
dyUtils.c dyFactorizedFlux.c\
myIO.c dyNQCCorrection.c \
mySpherical.c dyBHRingdown.c \
myUtils.c dyHybridRingdown.c

#PyEOBCal_la_LIBADD = -lgsl -lgslcblas
#PyEOBCal_la_LDFLAGS = -module -avoid-version

#all: $(pkgpyexec_LTLIBRARIES:.la=.so)

#.la.so:
	ln -s .libs/$@ $@ || true


libEOB_la_LIBADD = -lgsl -lgslcblas
libEOB_la_LDFLAGS = -module -avoid-version
MAIN_SOURCES = dyMain.c
MAIN_LDADD = libEOB.la
all:$(bin_PROGRAMS)
$(bin_PROGRAMS): $(lib_LTLIBRARIES) 	

