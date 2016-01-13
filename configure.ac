AC_INIT(RNAxplorer,1.0)
AM_INIT_AUTOMAKE(-Wall -Werror foreign)

AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AM_PROG_CC_C_O

#Output this variables to the makefiles

AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

#Check headers
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
AC_C_INLINE   #check if C compiler uses INLINE

AC_CONFIG_HEADERS(config.h) 

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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

