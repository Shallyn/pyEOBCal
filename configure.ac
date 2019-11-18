#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pyEOBCal], [0.0.1], [shallyn.liu@foxmail.com])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([csrc/dyEvolution.c])
AC_CONFIG_HEADERS([csrc/config.h])
AC_CONFIG_FILES([ \
	Makefile \
	csrc/Makefile])

AM_INIT_AUTOMAKE([ \
  1.11 \
  foreign \
  subdir-objects \
  color-tests \
  parallel-tests \
  tar-ustar \
  dist-xz \
  no-dist-gzip \
])

# Checks for programs.
AC_PROG_CC
AC_PROG_LN_S

# Checks for libraries.
LT_INIT
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL(libtool)

# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-lintl':
AC_CHECK_LIB([intl], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])


dnl python checks (you can change the required python version bellow)
AM_PATH_PYTHON(3.6)

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow sqrt strchr])
AC_OUTPUT
