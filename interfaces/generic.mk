
SWIG_main_src = $(srcdir)/../RNAxplorer.i

#SWIG_tmaps = \
#  $(srcdir)/../typemaps.i
#  $(srcdir)/../ptr2array.i \
#  $(srcdir)/../md_globals_tmaps.i

#SWIG_misc_src = \
#  $(srcdir)/../compare.i \
#  $(srcdir)/../constraints.i \
#  $(srcdir)/../eval.i \
#  $(srcdir)/../fold_compound.i \
#  $(srcdir)/../inverse.i \
#  $(srcdir)/../mfe.i \
#  $(srcdir)/../model_details.i \
#  $(srcdir)/../params.i \
#  $(srcdir)/../part_func.i \
#  $(srcdir)/../plotting.i \
#  $(srcdir)/../subopt.i \
#  $(srcdir)/../utils.i

SWIG_module_name = RNAxplorer

SWIG_wrapper = RNAxplorer_wrap.cpp

SWIG_src = $(SWIG_main_src) #$(SWIG_misc_src) $(SWIG_tmaps)
