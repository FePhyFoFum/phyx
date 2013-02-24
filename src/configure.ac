#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([phyx],[1.0],[eebsmith@umich.edu])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

AC_CHECK_LIB([m], [fabs])

have_nlopt=no
AC_SEARCH_LIBS([nlopt_urand],[nlopt_cxx], [have_nlopt=yes])

if test "x${have_nlopt}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The nlopt library is required for optimization.
phyx was developed with nlopt compiled with these flags
 ./configure --with-cxx --without-octave --without-matlab
--------------------------------------------------------])
fi

have_arma=no
AC_SEARCH_LIBS([wrapper_sgesv_],[armadillo], [have_arma=yes])

if test "x${have_arma}" = xno; then
   AC_MSG_ERROR([
--------------------------------------------------------
The armadillo library is required for some linear algebra
operations.
phyx was developed with armadillo compiled with these flags
 ./configure
--------------------------------------------------------])
fi

# Checks for header files.
AC_CHECK_HEADERS([fenv.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
