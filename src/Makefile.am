
AM_CPPFLAGS = $(CFLAGS) $(RNALIB_CFLAGS) 

bin_PROGRAMS = RNAxplorer

# convenience library
noinst_LTLIBRARIES = libRNAxplorer_conv.la

LDADD = $(top_builddir)/src/libRNAxplorer_conv.la $(RNALIB_LIBS)

GENGETOPT_CMDL =  RNAxplorer_cmdl.c RNAxplorer_cmdl.h
GENGETOPT_FILES =  RNAxplorer.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated content
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
%_cmdl.h %_cmdl.c: %.ggo
	$(AM_V_GEN)gengetopt -i $<

libRNAxplorer_conv_la_SOURCES = \
  RNAwalk.c \
  RNAwalk.h \
  meshpoint.c \
  meshpoint.h \
  distorted_sampling.c \
  distorted_sampling.h \
  distorted_samplingMD.c \
  distorted_samplingMD.h \
  barrier_lower_bound.c \
  barrier_lower_bound.h \
  paths.c \
  paths.h

RNAxplorer_SOURCES = ${GENGETOPT_CMDL} RNAxplorer.c
