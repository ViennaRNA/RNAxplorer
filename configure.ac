AC_INIT([RNAxplorer],[0.9],[ronny@tbi.univie.ac.at],[RNAxplorer])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)

AC_SUBST(VERSION)

AC_CONFIG_MACRO_DIR([m4])

AM_SILENT_RULES([yes])

AM_INIT_AUTOMAKE(-Wall -Werror foreign)


#Create debug conditional
# prevent setting xFLAGS to default of -g -O2
if test x"$CFLAGS" = x""; then
  AC_SUBST(CFLAGS, [ ])
fi
if test x"$CPPFLAGS" = x""; then
  AC_SUBST(CPPFLAGS, [ ])
fi

AC_ARG_ENABLE([debug],
              AC_HELP_STRING([--enable-debug],
              [compile with debug symbols @<:@default=no@:>@]),
              [want_debug="$enableval"], [want_debug=no])

if test "x$want_debug" = "xyes"; then
  CFLAGS="$CFLAGS -O0 -ggdb"
  CPPFLAGS="$CPPFLAGS -O0 -ggdb"
  AC_DEFINE([DEBUG], 1, [Define for debugging])
else
  CFLAGS="$CFLAGS -O2"
  CPPFLAGS="$CPPFLAGS -O2"
fi

#add lapack, (needed for multidimensional sampling)
AX_LAPACK([
      AC_CHECK_HEADERS([openblas/lapacke.h lapacke/lapacke.h lapacke.h],[rx_lapacke_header=yes;break;],)
      AS_IF([test "x$rx_lapacke_header" = "xyes"],[
        LDFLAGS="$LDFLAGS -llapack"
      ],[
        AC_MSG_ERROR([could not find lapacke.h])
      ])
      ],[AC_MSG_ERROR([could not find required version of Lapack library])])

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AM_PROG_CC_C_O

#Output this variables to the makefiles

AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

#Check headers
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
AC_C_INLINE   #check if C compiler uses INLINE

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config.h])

AM_PROG_AR

AC_PROG_LIBTOOL
LT_INIT([pic-only, static])
AC_SUBST(LIBTOOL_DEPS)

dnl Check for OpenMP and prepare flags
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"        

#Search for the RNAlib2 and set variables for Makefiles if found.  _CFLAGS and _LIBS can be overriden with the modue Name from PKG_CHECK_MODULES
PKG_CHECK_MODULES([RNALIB], [RNAlib2 >= 2.2.0],
                  [ AC_SUBST([RNALIB_CFLAGS])
                    AC_SUBST([RNALIB_LIBS])
                  ],
                  [])

dnl Checks for libraries.
AC_CHECK_LIB([m], [exp])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h float.h limits.h stdlib.h string.h strings.h unistd.h unistd.h math.h])

dnl Checks for funtions
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor strdup strstr strchr strrchr strstr strtol strtoul pow rint sqrt erand48 memset memmove erand48])

AC_RXP_INIT

AC_OUTPUT

AC_RXP_NOTICE

