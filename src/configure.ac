#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([phyx],[1.0],[eebsmith@umich.edu])

AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

#taking out fortran
#AC_PROG_FC

AC_CHECK_LIB([m], [fabs])

AC_PREFIX_DEFAULT(/usr/local)

have_nlopt=no
AC_SUBST([HNLOPT],[Y])
AC_CHECK_LIB(nlopt,nlopt_version,[have_nlopt=yes],[have_nlopt=no])
AC_CHECK_HEADER(nlopt.h,[have_nlopt=yes],[have_nlopt=no])
if test "x${have_nlopt}" = xno; then
   AC_MSG_WARN([
--------------------------------------------------------
The nlopt library is required for optimization for some
programs. You can install from standard repos or you 
can compile with these flags
 ./configure --without-octave --without-matlab
--------------------------------------------------------])
   AC_SUBST([HNLOPT],[N])
fi

#need to find a better lib to search for
have_arma=no
AC_SUBST([HARM],[Y])
AC_SEARCH_LIBS([wrapper_sgesv_], [armadillo], [have_arma=yes])

if test "x${have_arma}" = xno; then
   AC_MSG_WARN([
--------------------------------------------------------
The armadillo library is required for some linear algebra
operations.
phyx was developed with armadillo compiled with these flags
 ./configure
--------------------------------------------------------])
    AC_SUBST([HARM],[N])
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

#############################
# Platform specific setup
#############################

OPT_FLAGS="-O3 -std=c++11"

AC_SUBST([HOMP],[N])
if `$CXX -v 2>&1 | grep 'gcc version' >/dev/null 2>&1` ; then
   OPT_FLAGS+=" -fopenmp"
   AC_SUBST([HOMP],[Y])
fi

#AC_CANONICAL_HOST
## Check for which host we are on and setup a few things
## specifically based on the host
#case $host_os in
#  darwin* )
#        # Do something specific for mac
#        ;;
#  linux*)
#        # Do something specific for linux
#        OPT_FLAGS="$OPT_FLAGS -march=native"
#        ;;
#   *BSD*)
#        # Do something specific for BSD
#        ;;
#    *)
#        #Default Case
#        AC_MSG_ERROR([Your platform is not currently supported])
#        ;;
#esac

AC_SUBST(OPT_FLAGS)

#############################

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
