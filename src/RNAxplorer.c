#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/findpath.h>
#include <ViennaRNA/2Dfold.h>
#include <ViennaRNA/2Dpfold.h>
#include <ViennaRNA/mm.h>

#include "RNAxplorer_cmdl.h"
#include "RNAwalk.h"
#include "meshpoint.h"
#include "barrier_lower_bound.h"
#include "distorted_sampling.h"
#include "distorted_samplingMD.h"
#include "paths.h"
#include "RNAxplorer.h"

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

/**
 *** \file RNAxplorer.c
 **/

int whatToDo = FIND_BEST_FOLDINGPATH;

static char *seq;
int maxKeep = 100;
int maxIterations = 1;
int maxStorage = 10;
int maximum_distance1 = 5;
int maximum_distance2 = 5;

float betascale = 1.0;

int method;
extern int circ;
static char scale1[] = "....,....1....,....2....,....3....,....4";
static char scale2[] = "....,....5....,....6....,....7....,....8";

static char *extended_options = NULL;

int main(int argc, char *argv[]) {
  struct RNAxplorer_args_info args_info;

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if(RNAxplorer_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if(args_info.temp_given)
    temperature = args_info.temp_arg;

  /* method */
  if(args_info.method_given){
    char *m = args_info.method_arg;
    if(!strcmp(m, "MC"))
      method = MC_METROPOLIS;
    else if(!strcmp(m, "MC-SA")){
      method = MC_METROPOLIS;
      simulatedAnnealing = 1;
    }
    else if(!strcmp(m, "GW"))
      method = GRADIENT_WALK;
    else if(!strcmp(m, "DB-MFE"))
      whatToDo = FIND_DISTANCE_BASED_MFE_PATH;
    else if(!strcmp(m, "BLUBB"))
      whatToDo = FIND_2D_BARRIER_ESTIMATE;
    else if(!strcmp(m, "SM"))
      whatToDo = FIND_2D_LANDSCAPE_ESTIMATE;
  }

  /* maximum number of simulations / iterations */
  if(args_info.extended_opt_given)
    extended_options = strdup(args_info.extended_opt_arg);

  /* maximum number of simulations / iterations */
  if(args_info.iterations_given)
    maxIterations = args_info.iterations_arg;

  /* maxkeep for Flamm et al. direct path heuristics */
  if(args_info.maxKeep_given)
    maxKeep = args_info.maxKeep_arg;

  /* Amount of best solutions to store per iteration */
  if(args_info.maxStore_given)
    maxStorage = args_info.maxStore_arg;

  if(args_info.circ_given)
    circ = 1;

  if(args_info.cooling_rate_given)
    treduction = args_info.cooling_rate_arg;

  if(args_info.tstart_given)
    tstart = args_info.tstart_arg + K0;

  if(args_info.tstop_given)
    tstop = args_info.tstop_arg + K0;

  if(args_info.penalizeBackWalks_given)
    backWalkPenalty = 1;

  if(args_info.basinStructure_given)
    whatToDo = FIND_BASIN_STRUCTURE;

  if(args_info.maxDist_given)
    maximum_distance1 = maximum_distance2 = args_info.maxDist_arg;

  if(args_info.maxDist1_given)
    maximum_distance1 = args_info.maxDist1_arg;

  if(args_info.maxDist2_given)
    maximum_distance2 = args_info.maxDist2_arg;

  if(args_info.betaScale_given)
    betascale = args_info.betaScale_arg;

  /* free allocated memory of command line data structure */
  RNAxplorer_cmdline_parser_free(&args_info);

  switch (whatToDo) {
  case FIND_BASIN_STRUCTURE:
    GetBasinStructure();
    break;

  default:
    RNAxplorer();
    break;
  }
  return (EXIT_SUCCESS);
}

void GetBasinStructure(void) {
  char *s1;
  seq = get_line(stdin);
  s1 = get_line(stdin);

  vrna_md_t md;
  vrna_md_set_default(&md);
  /* set user-defined model details */
  md.circ = circ;
  md.uniq_ML = 1;

  initRNAWalk(seq, &md);
  char *basinStructure = structureWalk(seq, s1, GRADIENT_WALK);
  fprintf(stdout, "%s\n", basinStructure);
}

void RNAxplorer() {
  char *s1 = NULL, *s2 = NULL, *line, *start_struct = NULL, *target_struct = NULL;
  int istty, n;

  istty = isatty(fileno(stdout)) && isatty(fileno(stdin));
  do{
    if(istty){
      printf(
          "\nInput strings\n1st line: sequence (upper or lower case)\n2nd + 3rd line: start and target structure (dot bracket notation)\n@ to quit\n");
      printf("%s%s\n", scale1, scale2);
    }
    if((line = get_line(stdin)) == NULL)
      break;

    /* skip comment lines and get filenames */
    while((*line == '*') || (*line == '\0') || (*line == '>')){
      if(*line == '>')
        printf("%s\n", line);
      free(line);
      if((line = get_line(stdin)) == NULL)
        break;
    }

    if((line == NULL) || (strcmp(line, "@") == 0))
      break;

    seq = (char *) vrna_alloc(strlen(line) + 1);
    (void) sscanf(line, "%s", seq);
    free(line);
    n = (int) strlen(seq);

    s1 = (char *) vrna_alloc(sizeof(char) * (n + 1));
    s2 = (char *) vrna_alloc(sizeof(char) * (n + 1));

    if((start_struct = get_line(stdin)) == NULL)
      vrna_message_error("1st structure missing\n");
    else if(strlen(start_struct) != n)
      vrna_message_error("sequence and 1st structure have unequal length");
    strcpy(s1, start_struct);

    if((target_struct = get_line(stdin)) == NULL)
      vrna_message_error("2nd structure missing\n");
    else if(strlen(target_struct) != n)
      vrna_message_error("sequence and 2nd et structure have unequal length");
    strcpy(s2, target_struct);
    if(istty)
      printf("length = %d\n", n);

    size_t numberOfReferences = 2;
    char ** totalReferences = (char**) vrna_alloc(numberOfReferences * sizeof(char*) + 1);
    totalReferences[0] = s1;
    totalReferences[1] = s2;
    //scan stdin for additional references
    char * tmpStruct;
    while((tmpStruct = get_line(stdin)) != NULL){
      if(strlen(tmpStruct) == n){
        char * newReference = (char *) vrna_alloc(sizeof(char) * (n + 1));
        strcpy(newReference, tmpStruct);
        numberOfReferences++;
        totalReferences = vrna_realloc(totalReferences, numberOfReferences * sizeof(char*));
        totalReferences[numberOfReferences - 1] = newReference;
      }
      free(tmpStruct);
    }

    /* get fold compound for MFE and PF computation */
    //double mfe;
    vrna_md_t md;
    vrna_md_set_default(&md);
    /* set user-defined model details */
    md.circ = circ;
    md.uniq_ML = 1; /* in case we need M1 arrays */
    md.compute_bpp = 0;
    md.betaScale = betascale;

    vrna_fold_compound_t *vc = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

    fprintf(stdout, "%s\n", seq);
    for(int i = 0; i < numberOfReferences; i++){
      fprintf(stdout, "%s %6.2f\n", totalReferences[i], vrna_eval_structure(vc, totalReferences[i]));
    }

    vrna_path_t *foldingPath, *Saddle, *r;

    switch (whatToDo) {
    case FIND_2D_BARRIER_ESTIMATE: {
      barrier_estimate_2D(seq, &md, s1, s2, maximum_distance1, maximum_distance2);
    }
      break;

    case FIND_2D_LANDSCAPE_ESTIMATE: {
      //gridLandscapeT* grid = estimate_landscape(vc, totalReferences, numberOfReferences, maxIterations, extended_options);
      gridLandscapeT* grid = estimate_landscapeMD(vc, totalReferences, numberOfReferences, maxIterations,
          extended_options);
      printLandscape(grid, vc);
      free_gridLandscape(grid);
    }
      break;

    case FIND_DISTANCE_BASED_MFE_PATH: {
      foldingPath = get_path(seq, s1, s2, maxKeep/*, &numSteps, circ*/);
      Saddle = getSaddlePoint(foldingPath/*, numSteps*/);

      fprintf(stdout, "# direct Path:\n# barrier: %6.2f\n\n", Saddle->en);
      for(r = foldingPath; r->s; r++){
        fprintf(stdout, "%s %6.2f\n", r->s, r->en);
      }
      free_path(foldingPath);
      foldingPath = NULL;

      /* this was the old way to compute the folding path... now the new one is following */
      fprintf(stdout, "\n# searching for alternative paths\n# ");
      fflush(stdout);
      foldingPath = levelSaddlePoint2(seq, s1, s2/*, &numSteps*/, 0, maxIterations, maxKeep, maxStorage,
          maximum_distance1, maximum_distance2);
      fprintf(stdout, "\n# done\n\n# Path with detours:\n# barrier: %6.2f\n\n",
          getSaddlePoint(foldingPath/*, numSteps*/)->en);
      for(r = foldingPath; r->s; r++){
        fprintf(stdout, "%s %6.2f\n", r->s, r->en);
      }
      free_path(foldingPath);
      foldingPath = NULL;
      fflush(stdout);
    }
      break;
    default:
      vrna_init_rand();
      levelSaddlePoint(seq, s1, s2, maxIterations, maxKeep, method, maxStorage);
      break;
    }

    vrna_fold_compound_free(vc);

    free(seq);
    for(int i = 0; i < numberOfReferences; i++){
      free(totalReferences[i]);
    }
    free(totalReferences);
    free(start_struct);
    free(target_struct);
  } while(1);
}

/* taken from ivo's "neighbor.c" */
void print_structure(short* pt, int E) {
  int i;
  for(i = 1; i <= pt[0]; i++){
    if(pt[i] == 0){
      fprintf(stderr, ".");
      continue;
    }
    if(pt[i] > i){
      fprintf(stderr, "(");
      continue;
    }
    if(pt[i] < i){
      fprintf(stderr, ")");
      continue;
    }
  }
  fprintf(stderr, " %4d\n", E);
  fflush(stderr);
}
