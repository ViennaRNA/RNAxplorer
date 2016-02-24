include ../generic.mk

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py


pkgpyexecdir = $(pyexecdir)/RNAxplorer

pkgpyexec_LTLIBRARIES = _RNAxplorer.la
pkgpyexec_PYTHON = RNAxplorer/__init__.py

_RNAxplorer_la_SOURCES = $(top_srcdir)/interfaces/RNAxplorer.i $(SWIG_wrapper)

_RNAxplorer_la_CPPFLAGS = -I$(top_srcdir)/src $(PYTHON_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS)
_RNAxplorer_la_LIBADD = $(LDADD) $(RNALIB_LIBS) $(top_builddir)/src/libRNAxplorer_conv.la $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS)
#if WITH_LIBSVM
#_RNAxplorer_la_LIBADD += -lstdc++
#endif

_RNAxplorer_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PYTHON_EXTRA_LDFLAGS)
#if WITH_LTO
#_RNAxplorer_la_LDFLAGS += -Wc,$(LTO_LDFLAGS)
#endif

$(SWIG_wrapper): $(SWIG_src)
	$(AM_V_GEN)$(SWIG) -I$(top_srcdir)/src \
                  $(AX_SWIG_PYTHON_OPT) \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

# We rename RNAxplorer.py to __init__.py so that it can be installed
# in an "RNAxplorer" directory and be treated as a package.
RNAxplorer/__init__.py: $(SWIG_module_name).py
	-test -d RNAxplorer || $(MKDIR_P) RNAxplorer
	cp `test -f RNAxplorer.py || echo '$(srcdir)/'`RNAxplorer.py RNAxplorer/__init__.py

clean-local:
	-rm -rf RNAxplorer

# Remove the .la file - RNAxplorer.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpyexecdir)/_RNAxplorer.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pyexec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpyexecdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

