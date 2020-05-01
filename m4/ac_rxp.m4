# This file is part of Autoconf.                       -*- Autoconf -*-

# RNAxplorer 2016 Ronny Lorenz, Gregor Entzian
#

##----------------##
## Public macros. ##
##----------------##

AC_DEFUN([AC_RXP_INIT],[

AX_COMPILER_VENDOR
AC_CANONICAL_HOST

##--------------------##
## Enable scripting   ##
## language interface ##
##--------------------##

RXP_ENABLE_SWIG_INTERFACES


##------------------##
## Prepare files    ##
##------------------##

AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([Makefile src/Makefile])
])

AC_DEFUN([AC_RXP_NOTICE],[

# get directory paths

eval _bindir=$(eval printf "%s" $bindir)


AS_IF([test $with_python = "yes"],[
  eval _python3_arch_dir=$(eval printf "%s" ${py3execdir})
  eval _python3_lib_dir=$(eval printf "%s" ${python3dir})
  ],[
    _python3_arch_dir=""
    _python3_lib_dir=""
    _python3_install="Not to be installed"
])

# Notify the user

AC_MSG_NOTICE([

Configured successful with the following options:

RNAxplorerlib Scripting Interfaces:
  Python3 Interface:         ${with_python:-yes}     $python3_enabled_but_failed


Files will be installed in the following directories:

  Executables:       $_bindir
  Python3 Interface: $_python3_install
    (binaries):      $_python3_arch_dir
    (scripts):       $_python3_lib_dir
])

])



