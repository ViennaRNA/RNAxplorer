#ifndef __RNA_WALK__
#define __RNA_WALK__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fold_vars.h"

#define DEBUG
//#define _DEBUG_LOOPIDX

#define MIN2(A, B)        ((A) < (B) ? (A) : (B))
#define MAX2(A, B)        ((A) > (B) ? (A) : (B))

/* simulated annealing related variables */
extern int        simulatedAnnealing;
extern FLT_OR_DBL tstart;
extern FLT_OR_DBL tstop;
extern float      treduction;

/* Monte Carlo related variables */
extern int rememberStructures;
extern int maxRest;
extern int backWalkPenalty;

/* init function that initializes the random number generator
*  and also constructs the S and S1 sequence encoding arrays
*  and the pair table for further usage
*/
void initRNAWalk(char *seq, int circ);

/* free all arrays that were allocated by the init function
*/
void freeRNAWalkArrays(void);

/* the actual walking function that performs a walk on the
*  structure landscape according to the defined method
*  The returned secondary structure in dot bracket notation
*  is the one where the walk ended...
*  Implemented walks are defined below...
*/
#define GRADIENT_WALK   0
#define MC_METROPOLIS   1
char *structureWalk(char *seq, char *structure, int method, int circ);

/* a simple position finding function that searches
*  for the first entry in a sorted array of floats
* that exceeds a given value
*/
int getPosition(float *array, float value, int array_size);

/* this function is "stolen" from the pathfinder.c
*  of Ivo and Xtof
*/
static int *pair_table_to_loop_index (short *pt);

#endif
