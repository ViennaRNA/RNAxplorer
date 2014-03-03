CC = gcc
VRNA_INC = $(shell pkg-config --cflags "RNAlib2 >= 2.0")
VRNA_LIB = $(shell pkg-config --libs "RNAlib2 >= 2.0")

INCLUDES = $(VRNA_INC) -DLOOP_EN
LIBS    = ${VRNA_LIB} -lm
CFLAGS = -Wall -O3

# rules

OFILES =        RNAxplorer.o RNAwalk.o meshpoint.o RNAxplorer_cmdl.o
EXEFILE =       RNAxplorer




all:                            $(EXEFILE)

debug:				CFLAGS := $(CFLAGS) -g3

debug:				all

profile:			CFLAGS := $(CFLAGS) -pg -g3

profile:			all

all:                            $(EXEFILE)

clean:				
				rm -f $(OFILES) $(EXEFILE)

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
