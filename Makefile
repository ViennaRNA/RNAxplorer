CC = gcc
VRNA_INC = $(shell pkg-config --cflags "RNAlib2 >= 2.2")
VRNA_LIB = $(shell pkg-config --libs "RNAlib2 >= 2.2")
#VRNA_INC = -I/home/mescalin/ronny/WORK/ViennaRNA/src
#VRNA_LIB = -L/home/mescalin/ronny/WORK/ViennaRNA/src/ViennaRNA -lRNA -lgomp

INCLUDES = $(VRNA_INC) -DLOOP_EN
LIBS    = ${VRNA_LIB} -lm
CFLAGS = -Wall -O3

# rules

OFILES =        RNAxplorer.o RNAwalk.o meshpoint.o RNAxplorer_cmdl.o barrier_lower_bound.o distorted_sampling.o
EXEFILE =       RNAxplorer




all:                            $(EXEFILE)

debug:				CFLAGS := $(CFLAGS) -g3

debug:				all

profile:			CFLAGS := $(CFLAGS) -pg -g3

profile:			all

all:                            $(EXEFILE)

clean:				
				rm -f $(OFILES) $(EXEFILE)

barrier_lower_bound.o:		barrier_lower_bound.c barrier_lower_bound.h
				$(CC) -c $(CFLAGS) $(INCLUDES) barrier_lower_bound.c

distorted_sampling.o:		distorted_sampling.c distorted_sampling.h
				$(CC) -c $(CFLAGS) $(INCLUDES) distorted_sampling.c

RNAxplorer.o:			RNAxplorer.c RNAxplorer.h meshpoint.h RNAwalk.h RNAxplorer_cmdl.h
				$(CC) -c $(CFLAGS) $(INCLUDES) RNAxplorer.c

meshpoint.o:			meshpoint.c meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) meshpoint.c

RNAwalk.o:			RNAwalk.c RNAwalk.h meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) RNAwalk.c

RNAxplorer_cmdl.h:		RNAxplorer.ggo
				gengetopt -i RNAxplorer.ggo

RNAxplorer_cmdl.o:		RNAxplorer_cmdl.h RNAxplorer.ggo
				$(CC) -c $(CFLAGS) $(INCLUDES) RNAxplorer_cmdl.c

RNAxplorer: 			$(OFILES)
				$(CC) $(CFLAGS) -o $(EXEFILE) $(OFILES) $(LIBS)

# End of file
