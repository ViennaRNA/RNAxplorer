
AC_DEFUN([RXP_ENABLE_SWIG_INTERFACES],[

  AX_REQUIRE_DEFINED([AX_PKG_SWIG])

  RXP_ADD_PACKAGE([swig],
                  [SWIG scripting language interfaces],
                  [yes],
                  [with_swig=no],
                  [with_swig=yes],
                  [${srcdir}/interfaces/Makefile.am])

  AS_IF([test "x$with_swig" != "xno"],[
    AX_PKG_SWIG(2.0.0, [has_swig="yes"], [has_swig="no"])
  ])

  RXP_ENABLE_SWIG_PYTHON

])

AC_DEFUN([RXP_ENABLE_SWIG_PYTHON],[

  RXP_ADD_PACKAGE([python],
                  [Python interface],
                  [yes],
                  [with_python=no],
                  [with_python=yes],
                  [${srcdir}/interfaces/Python/Makefile.am])


  ## check for python requirements
  AS_IF([test "x$with_python" != "xno"],[
    ## if swig is not available, check whether we already have swig generated sources
    if test "x$has_swig" != "xyes"
    then
      AC_RXP_TEST_FILE([${srcdir}/interfaces/Python/RNAxplorer_wrap.c],[],[
        with_python="no"
      ])
      AC_RXP_TEST_FILE([${srcdir}/interfaces/Python/RNAxplorer.py],[],[
        with_python="no"
      ])
    fi
  ])

  AS_IF([test "x$with_python" != "xno"],[
    AX_PYTHON_DEVEL([< '3.0.0'])
    AM_PATH_PYTHON
    AX_SWIG_PYTHON
##    pythondir=$PYTHON_SITE_PKG
##    pyexecdir=$PYTHON_SITE_PKG_EXEC

    AC_SUBST(PYTHONDIR,$pythondir)
    AC_SUBST(PKGPYTHONDIR,$pkgpythondir)
    AC_SUBST(PYEXECDIR,$pyexecdir)
    AC_SUBST(PKGPYEXECDIR,$pkgpyexecdir)

    AC_DEFINE([WITH_PYTHON_INTERFACE], [1], [Create the python interface to RNAxplorerlib])
    AC_SUBST([PYTHON_INTERFACE], [Python])
    AC_CONFIG_FILES([interfaces/Python/Makefile])
  ])
])
