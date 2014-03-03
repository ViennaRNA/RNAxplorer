CC = gcc
VRNA_INC = $(shell pkg-config --cflags "RNAlib2 >= 2.0")
VRNA_LIB = $(shell pkg-config --libs "RNAlib2 >= 2.0")

INCLUDES = $(VRNA_INC) -DLOOP_EN
LIBS    = ${VRNA_LIB} -lm
CFLAGS = -Wall -O3

# rules

OFILES =        PathFinder.o RNAwalk.o meshpoint.o Pathfinder_cmdl.o
EXEFILE =       PathFinder




all:                            $(EXEFILE)

debug:				CFLAGS := $(CFLAGS) -g3

debug:				all

profile:			CFLAGS := $(CFLAGS) -pg -g3

profile:			all

all:                            $(EXEFILE)

clean:				
				rm -f $(OFILES) $(EXEFILE)

PathFinder.o:			PathFinder.c PathFinder.h meshpoint.h RNAwalk.h Pathfinder_cmdl.h
				$(CC) -c $(CFLAGS) $(INCLUDES) PathFinder.c

meshpoint.o:			meshpoint.c meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) meshpoint.c

RNAwalk.o:			RNAwalk.c RNAwalk.h meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) RNAwalk.c

Pathfinder_cmdl.h:		Pathfinder.ggo
				gengetopt -i Pathfinder.ggo

Pathfinder_cmdl.o:		Pathfinder_cmdl.h Pathfinder.ggo
				$(CC) -c $(CFLAGS) $(INCLUDES) Pathfinder_cmdl.c

PathFinder: 			$(OFILES)
				$(CC) $(CFLAGS) -o $(EXEFILE) $(OFILES) $(LIBS)

# End of file
