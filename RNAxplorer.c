#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/findpath.h>
#include <ViennaRNA/2Dfold.h>
#include <ViennaRNA/2Dpfold.h>
#include <ViennaRNA/mm.h>

#include "RNAxplorer_cmdl.h"
#include "RNAwalk.h"
#include "meshpoint.h"
#include "RNAxplorer.h"

#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

/**
*** \file RNAxplorer.c
**/


int         whatToDo = FIND_BEST_FOLDINGPATH;

static char *seq;
int maxKeep = 100;
int maxIterations = 1;
int maxStorage = 10;
int maximum_distance1 = 5;
int maximum_distance2 = 5;
  
int curr_iteration = 0;  

int method;
extern  int circ;
static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

typedef struct nb_t{

  float B;
  int   n_lu;
  int   n_ld;
  int   n_ru;
  int   n_rd;
  int   k;
  int   l;
  int   prev_idx;

  float en;
  char *s;
} nb_t;


typedef struct position{
  int i;
  int j;
} position;

void barrier_estimate_2D(char *seq,char *s1, char *s2);
void printBarrier(float B, float E, char *s);


int main(int argc, char *argv[]) {
  
  struct RNAxplorer_args_info args_info;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAxplorer_cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  /* temperature */
  if(args_info.temp_given) temperature = args_info.temp_arg;

  /* method */
  if(args_info.method_given){
    char *m = args_info.method_arg;
    if(!strcmp(m, "MC"))
      method = MC_METROPOLIS;
    else if(!strcmp(m, "MC-SA")){
      method = MC_METROPOLIS;
      simulatedAnnealing = 1;
    } else if(!strcmp(m, "GW"))
      method = GRADIENT_WALK;
    else if(!strcmp(m, "DB-MFE"))
      whatToDo = FIND_DISTANCE_BASED_MFE_PATH;
    else if(!strcmp(m, "BLUBB"))
      whatToDo = FIND_2D_BARRIER_ESTIMATE;
  }

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



  /* free allocated memory of command line data structure */
  RNAxplorer_cmdline_parser_free (&args_info);

  switch(whatToDo){
    case FIND_BASIN_STRUCTURE:  GetBasinStructure();
                                break;
    
    default:                    RNAxplorer();
                                break;
  }
  return(EXIT_SUCCESS);
}


void barrier_estimate_2D(char *seq,char *s1, char *s2){
  short *pt1, *pt2;
  pt1 = vrna_pt_get(s1);
  pt2 = vrna_pt_get(s2);
  int i, n = pt1[0];
  /* compute symmetrical difference between both structures */
  int a = 0;
  int b = 0;
  int c = 0;
  for(i=1; i<=n; i++){
    if(pt1[i] == pt2[i]){
      if(pt1[i] != 0) c++;
    }
    else{
      if(i < pt1[i]) a++;
      if(i < pt2[i]) b++;
    }  
  }
  printf("%d:%d ... %d:%d\n", a+c, a+maximum_distance1, b+c, b+maximum_distance2);
  int maxD1 = (a+c < a+maximum_distance1) ? a+maximum_distance1 : a+c;
  int maxD2 = (b+c < b+maximum_distance2) ? b+maximum_distance2 : a+c;
  
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.circ     = circ;
  md.uniq_ML  = 1;

  vrna_fold_compound *vc = vrna_get_fold_compound_2D(seq, s1, s2, &md, VRNA_OPTION_MFE | VRNA_OPTION_DIST_CLASS);
  vrna_sol_TwoD_t *mfe_s = vrna_TwoD_fold(vc, maxD1, maxD2);

  nb_t **neighbors = (nb_t **)vrna_alloc((maxD1+1) * sizeof(nb_t *));
  int number_of_states = 0;

  /* make a lucky guess for the real max distancies */
  for(i=0;i<(maxD1+1);i++){
    neighbors[i] = (nb_t *)vrna_alloc((maxD2+1) * sizeof(nb_t));
  }
  float mfe_s1, mfe_s2;
  int map_s1, map_s2;

  int max_k = 0;
  int min_k = INF;
  int *max_l = (int *)vrna_alloc(sizeof(int) * (maxD1 + 1));
  int *min_l = (int *)vrna_alloc(sizeof(int) * (maxD1 + 1));
  for(i=0;i<(maxD1+1);i++){
    max_l[i] = 0;
    min_l[i] = INF;
  }
  int **mapping = (int **)vrna_alloc(sizeof(int *) * (maxD1 + 1));
  for(i=0;i<(maxD1+1);i++){
    mapping[i] = (int *)vrna_alloc(sizeof(int) * (maxD2 + 1));
  }

  /* get some statistics of the 2D fold output */
  for(i = number_of_states = 0; mfe_s[i].k != INF; i++){
    int k = mfe_s[i].k;
    int l = mfe_s[i].l;

    if(k == -1) continue;

    max_k = MAX2(max_k, k);
    min_k = MIN2(min_k, k);
    max_l[k] = MAX2(max_l[k], l);
    min_l[k] = MIN2(min_l[k], l);

    if(k == 0){ map_s1 = i; mfe_s1 = mfe_s[i].en;}
    if(l == 0){ map_s2 = i; mfe_s2 = mfe_s[i].en;}

    mapping[k][l] = i;
    number_of_states++;
  }

  /* begin actual graph node construction */
  nb_t *nodes = (nb_t *)vrna_alloc(sizeof(nb_t) * (number_of_states + 1));
  int map_source = (mfe_s1 < mfe_s2) ? map_s1 : map_s2;
  int map_target = (mfe_s2 <= mfe_s1) ? map_s1 : map_s2;
  float min12 = MIN2(mfe_s1, mfe_s2);

  for(i = 0; mfe_s[i].k != INF; i++){
    int k = mfe_s[i].k;
    int l = mfe_s[i].l;

    if(k == -1) continue;

    nodes[i].en       = mfe_s[i].en;
    nodes[i].s        = mfe_s[i].s;
    nodes[i].k        = k;
    nodes[i].l        = l;
    nodes[i].n_lu     = nodes[i].n_ld = nodes[i].n_ru = nodes[i].n_rd = -1;
    nodes[i].B        = (float)INF/100.;
    nodes[i].prev_idx = -1;

    /* set the four neighbor indices, if they are within range */
    if(k - 1 >= min_k){
      if((l - 1 >= min_l[k-1]) && (l - 1 <= max_l[k-1])){
        nodes[i].n_ld = mapping[k-1][l-1];
      }
      if((l + 1 <= max_l[k-1]) && (l + 1 >= min_l[k-1])){
        nodes[i].n_lu = mapping[k-1][l+1];
      }
    }
    if(k + 1 <= max_k){
    
      if((l - 1 >= min_l[k+1]) && (l - 1 <= max_l[k+1])){
        nodes[i].n_rd = mapping[k+1][l-1];
      }
      if((l + 1 <= max_l[k+1]) && (l + 1 >= min_l[k+1])){
        nodes[i].n_ru = mapping[k+1][l+1];
      }
    }
  }

  nb_t *vertices = nodes;
  int vertex_count = number_of_states;
  int shift = 0;

  vertices[map_source].B = 0.;

  /* search for optimal path using Dijkstra algorithm */
  do{
    /* put node with smallest barrier to front */
    int min_idx = 0;
    int min_b   = vertices[0].B;
    for(i = 1; i < vertex_count; i++){
      if(vertices[i].B < min_b){
        min_idx = i;
        min_b = vertices[i].B;
      }
    }
    if(min_idx != 0){ /* swap entries if vertex with smallest barrier was not in front */
      nb_t tmp;
      memcpy(&tmp, vertices, sizeof(nb_t));
      memcpy(vertices, vertices + min_idx, sizeof(nb_t));
      memcpy(vertices + min_idx, &tmp, sizeof(nb_t));

      /* update all vertices that had either of the swapped nodes as neighbors */
      for(i = 0; i < vertex_count; i++){
        if(vertices[i].n_ld == 0)
          vertices[i].n_ld = min_idx;
        else if(vertices[i].n_ld == min_idx)
          vertices[i].n_ld = 0;

        if(vertices[i].n_lu == 0)
          vertices[i].n_lu = min_idx;
        else if(vertices[i].n_lu == min_idx)
          vertices[i].n_lu = 0;

        if(vertices[i].n_rd == 0)
          vertices[i].n_rd = min_idx;
        else if(vertices[i].n_rd == min_idx)
          vertices[i].n_rd = 0;

        if(vertices[i].n_ru == 0)
          vertices[i].n_ru = min_idx;
        else if(vertices[i].n_ru == min_idx)
          vertices[i].n_ru = 0;
      }

      /* did we move the source or target vertex? */
      if(min_idx == map_target){
        map_target = 0;
        break;
      }
      if(min_idx == map_source){
        map_source = 0;
      }
      if(map_target == 0){
        map_target = min_idx;
      }
    }

    /* now lets update the neighboring vertices */
    if(vertices[0].n_ld > 0){
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_ld].en - min12);
      if(tmp_b < vertices[vertices[0].n_ld].B){
        vertices[vertices[0].n_ld].B = tmp_b;
        vertices[vertices[0].n_ld].prev_idx = shift;
      }
    }

    if(vertices[0].n_lu > 0){
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_lu].en - min12);
      if(tmp_b < vertices[vertices[0].n_lu].B){
        vertices[vertices[0].n_lu].B = tmp_b;
        vertices[vertices[0].n_lu].prev_idx = shift;
      }
    }

    if(vertices[0].n_rd > 0){
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_rd].en - min12);
      if(tmp_b < vertices[vertices[0].n_rd].B){
        vertices[vertices[0].n_rd].B = tmp_b;
        vertices[vertices[0].n_rd].prev_idx = shift;
      }
    }

    if(vertices[0].n_ru > 0){
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_ru].en - min12);
      if(tmp_b < vertices[vertices[0].n_ru].B){
        vertices[vertices[0].n_ru].B = tmp_b;
        vertices[vertices[0].n_ru].prev_idx = shift;
      }
    }

    vertices++;
    vertex_count--;
    map_target--;
    shift++;
    /* decrease vertex index of neighboring nodes for remaining set of vertices */
    for(i = 0; i < vertex_count; i++){
      vertices[i].n_ld--;
      vertices[i].n_lu--;
      vertices[i].n_rd--;
      vertices[i].n_ru--;
    }
  } while(vertex_count > 0);

  /* now backtrack the optimal path */
  int curr        = map_target + shift;
  float B         = nodes[curr].B;
  float E_saddle  = nodes[curr].en;
  char  *s_saddle  = nodes[curr].s;

  while(curr != 0){
    if(nodes[curr].en > E_saddle){
      E_saddle = nodes[curr].en;
      s_saddle = nodes[curr].s;
    }
    printf("%d %d -> %6.2f\n", nodes[curr].k, nodes[curr].l, nodes[curr].B);
    curr = nodes[curr].prev_idx;
  }
  printf("%d %d -> %6.2f\n", nodes[0].k, nodes[0].l, nodes[0].B);

  printBarrier(B, E_saddle, s_saddle);

  vrna_free_fold_compound(vc);
  free(nodes);
  for(i=0;i<(maxD1+1);i++){
    free(mapping[i]);
  }
  free(mapping);
}

void printBarrier(float B, float E, char *s){
  if(s)
    printf("Estimated energy barrier is %6.2f\nSaddle: %s (%6.2f)\n", B, s, E);
  else
    printf("Estimated energy barrier is %6.2f\n", B);
}

void GetBasinStructure(void){
  char *s1;
  seq = get_line(stdin);
  s1 = get_line(stdin);
  initRNAWalk(seq, circ);
  char *basinStructure = structureWalk(seq, s1, GRADIENT_WALK, circ);
  fprintf(stdout, "%s\n", basinStructure);
}

void RNAxplorer(){
  char *s1 = NULL, *s2=NULL, *line, *start_struct=NULL, *target_struct=NULL;
  int istty, n;

  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
  do {
    if (istty) {
      printf("\nInput strings\n1st line: sequence (upper or lower case)\n2nd + 3rd line: start and target structure (dot bracket notation)\n@ to quit\n");
      printf("%s%s\n", scale1, scale2);
    }
    if ((line = get_line(stdin))==NULL) break;

    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      if (*line=='>')
        printf("%s\n", line);
        free(line);
        if ((line = get_line(stdin))==NULL) break;
      }

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    seq = (char *) vrna_alloc(strlen(line)+1);
    (void) sscanf(line,"%s",seq);
    free(line);
    n = (int) strlen(seq);

    s1 = (char *) vrna_alloc((unsigned) n+1);
    s2 = (char *) vrna_alloc((unsigned) n+1);

    if ((start_struct = get_line(stdin))==NULL)
      vrna_message_error("1st structure missing\n");
    else if(strlen(start_struct) != n)
      vrna_message_error("sequence and 1st structure have unequal length");
    strncpy(s2, start_struct, n);

    if ((target_struct = get_line(stdin))==NULL)
      vrna_message_error("2nd structure missing\n");
    else if(strlen(target_struct) != n)
      vrna_message_error("sequence and 2nd et structure have unequal length");
      strncpy(s1, target_struct, n);
    if (istty)
      printf("length = %d\n", n);

    /* get fold compound for MFE and PF computation */
    double mfe;
    vrna_md_t md;
    vrna_md_set_default(&md);
    /* set user-defined model details */
    md.circ     = circ;
    md.uniq_ML  = 1; /* in case we need M1 arrays */

    vrna_fold_compound *vc = vrna_get_fold_compound(seq, &md, VRNA_OPTION_MFE | VRNA_OPTION_PF);

    mfe = vrna_fold(vc, NULL);
    vrna_rescale_pf_params(vc, &mfe);
    (void)vrna_pf_fold(vc, NULL);

    fprintf(stdout, "%s\n", seq);
    fprintf(stdout, "%s %6.2f\n", s1, vrna_eval_structure(vc, s1));
    fprintf(stdout, "%s %6.2f\n", s2, vrna_eval_structure(vc, s2));

    path_t *foldingPath, *Saddle, *r;
    curr_iteration = 0;
    switch(whatToDo){
      case FIND_2D_BARRIER_ESTIMATE:      {
                                            barrier_estimate_2D(seq, s1, s2);
                                          }
                                          break;

      case FIND_DISTANCE_BASED_MFE_PATH:  {
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
                                            foldingPath = levelSaddlePoint2(s1, s2/*, &numSteps*/, 0);
                                            fprintf(stdout, "\n# done\n\n# Path with detours:\n# barrier: %6.2f\n\n", getSaddlePoint(foldingPath/*, numSteps*/)->en);
                                            for(r = foldingPath; r->s; r++){
                                              fprintf(stdout, "%s %6.2f\n", r->s, r->en);
                                             }
                                            free_path(foldingPath);
                                            foldingPath=NULL;
                                            fflush(stdout);
                                          }
                                          break;
      default:                            vrna_init_rand();
                                          levelSaddlePoint(s1, s2);
                                          break;
    }
  
    free(seq);
    free(s1);
    free(s2);
    free(start_struct);
    free(target_struct);
  } while(1);
}

void levelSaddlePoint(char *s1, char *s2){

  int iterator = maxIterations;
  path_t *foldingPath = get_path(seq, s1, s2, maxKeep/*, &steps, circ*/);
  path_t *Saddle = getSaddlePoint(foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  fprintf(stdout, "old Path:\nbarrier: %6.2f\n\n", Saddle->en);
  int d;
  for(d=0; foldingPath[d].s; d++){
    fprintf(stdout, "%s %6.2f\n", foldingPath[d].s, foldingPath[d].en);
  }
  path_t *newLeftSaddle, *newRightSaddle;
  path_t *path_left, *path_right;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list(&bestMeshPoints);
  fprintf(stdout, "\nsearching for alternative paths...");
  fflush(stdout);
  initRNAWalk(seq, circ);
  while(iterator > 0){
    char *newSaddle = structureWalk(seq, Saddle->s, method, circ);
    path_left = get_path(seq, s1, newSaddle, maxKeep/*, &steps1, circ*/);
    path_right = get_path(seq, newSaddle, s2, maxKeep/*, &steps2, circ*/);
  
    newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
    newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
    newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);
    iterator--;
    
    if(newSaddleEn < oldSaddleEn){
      insert_meshpoint(newSaddle, newSaddleEn, &bestMeshPoints, maxStorage);
    }
    free_path(path_left);
    free_path(path_right);
    fprintf(stdout, "\rsearching for alternative paths...(%2.1f %% done)", (float)(maxIterations-iterator)/(float)maxIterations * 100.);
    fflush(stdout);
  }
  if(bestMeshPoints.count > 0){
    path_left = get_path(seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
    path_right = get_path(seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);
 
    newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
    newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
    newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);

    fprintf(stdout, "\nnewPath with barrier: %6.2f\n\n", newSaddleEn);
    for(d=0; path_left[d].s; d++){
      fprintf(stdout, "%s %6.2f\n", path_left[d].s, path_left[d].en);
    }
    for(d=1; path_right[d].s; d++){
      fprintf(stdout, "%s %6.2f\n", path_right[d].s, path_right[d].en);
    }
  }
  else{
    fprintf(stdout, "no better path found...\n :-/\n");
  }
  clear_meshpoints(&bestMeshPoints);
}

path_t *levelSaddlePoint2(char *s1, char *s2/*, int *num_entry*/, int iteration){


  int i;
  path_t *foldingPath = get_path(seq, s1, s2, maxKeep/*, &steps, circ*/);
  path_t *Saddle = getSaddlePoint(foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  path_t *newLeftSaddle, *newRightSaddle;
  path_t *path_left = NULL;
  path_t *path_right = NULL;
  path_t *newPath = NULL;
  int steps1, steps2;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list(&bestMeshPoints);

  /* begin the nice pathfinding routine */
  short *pt1, *pt2;
  pt1 = vrna_pt_get(s1);
  pt2 = vrna_pt_get(s2);
  int n = pt1[0];
  
  short *intersect = (short *) vrna_alloc(sizeof(short)*(n+2));
  intersect[0] = n;
  /* compute symmetrical difference between both structures */
  int a = 0;
  int b = 0;
  for(i=1; i<=n; i++){
    if(pt1[i] == pt2[i]){
      intersect[i] = pt1[i];
      pt1[i] = pt2[i] = 0;
    }
    else{
      if(i < pt1[i]) a++;
      if(i < pt2[i]) b++;
      intersect[i] = 0;
    }  
  }

 
  /* first collect all meshpoints where we get an initially better path */
  if(a+b > 1){
    vrna_md_t md;
    vrna_md_set_default(&md);
    md.circ = circ;
    md.uniq_ML  = 1;

    vrna_fold_compound *vc = vrna_get_fold_compound_2D(seq, s1, s2, &md, VRNA_OPTION_MFE | VRNA_OPTION_DIST_CLASS);
    vrna_sol_TwoD_t *mfe_s = vrna_TwoD_fold(vc, a+maximum_distance1, b+maximum_distance2);

    vrna_free_fold_compound(vc);

    for(i=0;  mfe_s[i].k != INF; i++){
      if(mfe_s[i].k == -1){
        free(mfe_s[i].s);
        continue;
      }
      path_left = get_path(seq, s1, mfe_s[i].s, maxKeep/*, &steps1, circ*/);
      path_right = get_path(seq, mfe_s[i].s, s2, maxKeep/*, &steps2, circ*/);

      newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
      newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
      newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);
      if(newSaddleEn <= oldSaddleEn){
        insert_meshpoint_with_struct_energy(mfe_s[i].s, newSaddleEn, mfe_s[i].en, &bestMeshPoints, maxStorage);
      }
      free(mfe_s[i].s);
      free_path(path_left); path_left=NULL;
      free_path(path_right); path_right=NULL;
    }
    free(mfe_s);

#if 0
    for(i = 0; i<= a + maximum_distance1; i++){
      for(j = 0; j<= b + maximum_distance2; j++){
        if(dfold_structs[i][j].en != (float)INF/100.){
          path_left = get_path(seq, s1, dfold_structs[i][j].s, maxKeep/*, &steps1, circ*/);
          path_right = get_path(seq, dfold_structs[i][j].s, s2, maxKeep/*, &steps2, circ*/);
  
          newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
          newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
          newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);
          if(newSaddleEn <= oldSaddleEn){
            insert_meshpoint_with_struct_energy(dfold_structs[i][j].s, newSaddleEn, dfold_structs[i][j].en, &bestMeshPoints, maxStorage);
          }
          free_path(path_left); path_left=NULL;
          free_path(path_right); path_right=NULL;
        }
      }
      free(dfold_structs[i]);
    }
#endif
  }
  if(bestMeshPoints.count > 0){
    /*
    now as we know n better SaddlePoints, we can iterate deeper to maybe obtain an even better saddle
    */
    
    if(iteration < maxIterations){
      meshpoint *cur;
      int t_steps1, t_steps2;
      path_t *t_path_left = NULL, *t_path_right =  NULL;
      path_left = NULL;
      path_right = NULL;
      for(cur = bestMeshPoints.first; cur != NULL; cur = cur->next){
        t_path_left = levelSaddlePoint2(s1, cur->s, /*&t_steps1,*/ iteration+1);
        newLeftSaddle = getSaddlePoint(t_path_left/*, t_steps1*/);

        t_path_right = levelSaddlePoint2(cur->s, s2, /*&t_steps2,*/ iteration+1);
        newRightSaddle = getSaddlePoint(t_path_right/*, t_steps2*/);

        newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);

        if(newSaddleEn < oldSaddleEn){
          if(path_left) free_path(path_left);
          path_left = t_path_left;
          if(path_right) free_path(path_right);
          path_right = t_path_right;
          steps1 = t_steps1;
          steps2 = t_steps2;
          oldSaddleEn = newSaddleEn;
        }
        else{
          free_path(t_path_left);t_path_left=NULL;
          free_path(t_path_right);t_path_right=NULL;
        }
      }
      if(path_left == NULL || path_right == NULL){
        path_left = get_path(seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
        path_right = get_path(seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);
      
      }
    }
    /* if we are in the last iteration step, we just take the best found in this round... */
    else{
      path_left = get_path(seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
      path_right = get_path(seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);
 
      newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
      newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
      newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);
    
    }
    clear_meshpoints(&bestMeshPoints);
    for(steps1=0; path_left[steps1].s; steps1++);
    for(steps2=0; path_right[steps2].s; steps2++);
    newPath = (path_t *)vrna_alloc((steps1+steps2) * sizeof(path_t));
    memcpy((path_t *)newPath, (path_t *)path_left, steps1*sizeof(path_t));
    memcpy(((path_t *)newPath)+(steps1), ((path_t *)path_right) + 1, (steps2-1)*sizeof(path_t));
    if(steps2>0){
      free(path_right[0].s); /* since we skipped this entry and it never would be free'd */
    }
  }
  else{
    free(pt1); free(pt2); free(intersect);
    return foldingPath;
  }
  
  free(pt1); free(pt2); free(intersect);
  free_path(foldingPath);
  if(path_left) free(path_left); /* do not free the structures, since they are further uses in newPath */
  if(path_right) free(path_right); /* do not free the structures, since they are further uses in newPath */

  curr_iteration++;
  //fprintf(stdout, "\rsearching for alternative paths...(%2.1f %% done %d %d)", (float)(curr_iteration)/(float)pow(2*maxStorage, maxIterations) * 100.,curr_iteration, (int)pow(2*maxStorage, maxIterations)  );
  fprintf(stdout, ".");
  fflush(stdout);
  
  
  return newPath;
}

path_t *getSaddlePoint(path_t *foldingPath){
  path_t *r, *saddle;
  saddle = foldingPath;
  for(r = foldingPath; r->s; r++)
    if(saddle->en < r->en)
      saddle = r;
  //if(saddle->en == foldingPath->en) return --r;
  return saddle;
}


/* taken from ivo's "neighbor.c" */
void print_structure(short* pt, int E){
  int i;
  for (i=1; i<=pt[0]; i++) {
    if (pt[i] == 0) {
      fprintf(stderr, ".");
      continue;
    }
    if (pt[i] > i) {
      fprintf(stderr, "(");
      continue;
    }
    if (pt[i] < i) {
      fprintf(stderr, ")");
      continue;
    }
  }
  fprintf(stderr," %4d\n", E);
  fflush(stderr);
}
