CC = gcc
VRNADIR = ${HOME}/WORK/VRNA-git
INCLUDES = -I$(VRNADIR)/H -I$(VRNADIR)/lib -DLOOP_EN
LIBS    = -L$(VRNADIR)/lib -lRNA -lm
CFLAGS = -O3 -Wall
# rules

OFILES =        PathFinder.o RNAwalk.o meshpoint.o
EXEFILE =       PathFinder

all:                            $(EXEFILE)

debug:				CFLAGS := $(CFLAGS) -g3

debug:				all

debug_parallel:			CFLAGS := $(CFLAGS) -g3 -DUSE_OPENMP -fopenmp

debug_parallel:			all

profile:			CFLAGS := $(CFLAGS) -pg -g3

profile:			all

profile_parallel:		CFLAGS := $(CFLAGS) -g3 -pg -DUSE_OPENMP -fopenmp

profile_parallel:		all

parallel:			CFLAGS := $(CFLAGS) -DUSE_OPENMP -fopenmp

parallel:			all

all:                            $(EXEFILE)

clean:				
				rm -f $(OFILES) $(EXEFILE)

PathFinder.o:			PathFinder.c PathFinder.h meshpoint.h RNAwalk.h
				$(CC) -c $(CFLAGS) $(INCLUDES) PathFinder.c

meshpoint.o:			meshpoint.c meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) meshpoint.c

RNAwalk.o:			RNAwalk.c RNAwalk.h meshpoint.h
				$(CC) -c $(CFLAGS) $(INCLUDES) RNAwalk.c


PathFinder: 			$(OFILES)
				$(CC) $(CFLAGS) -o $(EXEFILE) $(OFILES) $(LIBS)

# End of file
