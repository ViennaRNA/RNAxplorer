
AM_CPPFLAGS = $(CFLAGS) $(RNALIB_CFLAGS) 

LDADD = $(RNALIB_LIBS)

bin_PROGRAMS = RNAxplorer

GENGETOPT_CMDL =  RNAxplorer_cmdl.c RNAxplorer_cmdl.h
GENGETOPT_FILES =  RNAxplorer.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated content
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
%_cmdl.h %_cmdl.c: %.ggo
	$(AM_V_GEN)gengetopt -i $<

RNAxplorer_SOURCES =  RNAxplorer_cmdl.c \
                      RNAxplorer.c \
                      RNAwalk.c \
                      RNAwalk.h \
                      meshpoint.c \
                      meshpoint.h \
                      distorted_sampling.c \
                      distorted_sampling.h \
                      barrier_lower_bound.c \
                      barrier_lower_bound.h
