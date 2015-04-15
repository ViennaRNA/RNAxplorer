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
    else if(!strcmp(m, "SM"))
      whatToDo = FIND_2D_LANDSCAPE_ESTIMATE;
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


/* landscape sampling via distortion */

typedef struct {
  double kT;
  int *idx;
  char *ref1;
  char *ref2;
  double x;
  double y;
  double *xs;
  double *ys;
  unsigned int *referenceBPs1;
  unsigned int *referenceBPs2;
  short *reference_pt1;
  short *reference_pt2;
} kl_soft_constraints;


kl_soft_constraints *kl_init_datastructures(vrna_fold_compound *vc, const char *s1, const char *s2, double x, double y){
  kl_soft_constraints *data;
  unsigned int i, n;
  char *s = vc->sequence;
  
  n = strlen(s);
  
  /* alloc all memory */
  data = (kl_soft_constraints *)vrna_alloc(sizeof(kl_soft_constraints));
  data->kT            = vc->exp_params->kT;
  data->idx           = vrna_get_iindx(n);
  data->ref1          = strdup(s1);
  data->ref2          = strdup(s2);
  data->reference_pt1 = vrna_pt_get(data->ref1);
  data->reference_pt2 = vrna_pt_get(data->ref2);
  data->referenceBPs1 = vrna_refBPcnt_matrix(data->reference_pt1, TURN); /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
  data->referenceBPs2 = vrna_refBPcnt_matrix(data->reference_pt2, TURN); /* matrix containing number of basepairs of reference structure2 in interval [i,j] */
  data->x             = x;
  data->y             = y;

  return data;
}

FLT_OR_DBL kl_pseudo_energy(int i, int j, int k, int l, char decomp, void *data){

  int d1, d2, ij, kl;
  int *idx                    = ((kl_soft_constraints *)data)->idx;
  short *reference_pt1        = ((kl_soft_constraints *)data)->reference_pt1;
  short *reference_pt2        = ((kl_soft_constraints *)data)->reference_pt2;
  unsigned int *referenceBPs1 = ((kl_soft_constraints *)data)->referenceBPs1;
  unsigned int *referenceBPs2 = ((kl_soft_constraints *)data)->referenceBPs2;
  double x                    = ((kl_soft_constraints *)data)->x;
  double y                    = ((kl_soft_constraints *)data)->y;

  int base_da = (reference_pt1[i] != (unsigned int)j) ? 1 : -1;
  int base_db = (reference_pt2[i] != (unsigned int)j) ? 1 : -1;
  ij = idx[i]-j;
  kl = idx[k]-l;
  d1 = d2 = 0;

  switch(decomp){
    case VRNA_DECOMP_PAIR_HP:     d1 = base_da + referenceBPs1[ij];
                                  d2 = base_db + referenceBPs2[ij];
                                  break;
    case VRNA_DECOMP_PAIR_IL:     d1 = base_da + referenceBPs1[ij] - referenceBPs1[kl];
                                  d2 = base_db + referenceBPs2[ij] - referenceBPs2[kl];
                                  break;
    case VRNA_DECOMP_PAIR_ML:     d1 = base_da + referenceBPs1[ij] - referenceBPs1[idx[i+1]-k+1] - referenceBPs1[idx[k]-j+1];
                                  d2 = base_db + referenceBPs2[ij] - referenceBPs2[idx[i+1]-k+1] - referenceBPs2[idx[k]-j+1];
                                  break;
    case VRNA_DECOMP_ML_UP_5:     d1 = referenceBPs1[ij] - referenceBPs1[idx[k]-j];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[k]-j];
                                  break;
    case VRNA_DECOMP_ML_UP_3:     d1 = referenceBPs1[ij] - referenceBPs1[idx[i]-k];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[i]-k];
                                  break;
    case VRNA_DECOMP_ML_UP:       d1 = referenceBPs1[ij];
                                  d2 = referenceBPs2[ij];
                                  break;
    case VRNA_DECOMP_ML_ML_ML:    d1 = referenceBPs1[ij] - referenceBPs1[idx[i]-k+1] - referenceBPs1[idx[k]-j] ;
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[i]-k+1] - referenceBPs2[idx[k]-j] ;
                                  break;
    case VRNA_DECOMP_EXT_UP_3:    d1 = referenceBPs1[ij] - referenceBPs1[idx[i]-k];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[i]-k];
                                  break;
    case VRNA_DECOMP_EXT_UP_5:    d1 = referenceBPs1[ij] - referenceBPs1[idx[k]-j];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[k]-j];
                                  break;
    case VRNA_DECOMP_EXT_UP:      d1 = referenceBPs1[ij];
                                  d2 = referenceBPs2[ij];
                                  break;
    case VRNA_DECOMP_EXT_EXT:     d1 = referenceBPs1[ij] - referenceBPs1[idx[i]-k] - referenceBPs1[idx[k+1]-j];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[i]-k] - referenceBPs2[idx[k+1]-j];
                                  break;
    case VRNA_DECOMP_EXT_STEM_UP: d1 = referenceBPs1[ij] - referenceBPs1[idx[i]-k];
                                  d2 = referenceBPs2[ij] - referenceBPs2[idx[i]-k];
                                  break;
  }

  return (x * d1 + y * d2)*100;
}

FLT_OR_DBL kl_exp_pseudo_energy(int i, int j, int k, int l, char decomp, void *data){

  double kT                   = ((kl_soft_constraints *)data)->kT;
//  printf("adding %g\n", kl_pseudo_energy(i,j,k,l,decomp,data));
  return exp((-10.*(double)kl_pseudo_energy(i,j,k,l,decomp,data))/kT);
}

typedef struct {
  int k;
  int l;
  int num_structs;
  int max_structs;
  char **structures;
  double pf;
  double mfe;
} gridpointT;


void
estimate_landscape( vrna_fold_compound *vc,
                    const char *s1,
                    const char *s2){

  unsigned int      *mm1, *mm2;
  int               i, j, MAX_k, MAX_l, length, *my_iindx;
  int               idx_1n;
  int               bp_ref1, bp_ref2, bp_mfe, bp_dist, bp_dist_mfe_ref1, bp_dist_mfe_ref2;
  float             e_ref1, e_ref2;
  double            mmfe, rescale;
  short             *pt_ref1, *pt_ref2, *pt_mfe;
  char              *s, *mfe_struct;
  double            kT, distortion_x, distortion_y, exp_ref1, exp_ref2;
  pf_paramT           *pf_parameters;

  s         = vc->sequence;
  length    = vc->length;
  my_iindx  = vc->iindx;
  idx_1n    = my_iindx[1] - length;
  pf_parameters   = vc->exp_params;
  kT              = pf_parameters->kT;


  /* get mfe for this sequence */
  mfe_struct  = (char *)vrna_alloc(sizeof(char)*(length+1));
  mmfe        = (double)vrna_fold(vc, mfe_struct);

  /* get free energies of the reference structures */
  e_ref1  = vrna_eval_structure(vc, s1);
  e_ref2  = vrna_eval_structure(vc, s2);

    vrna_init_rand();

  /* ######################################################################### */
  /* #### generate a pseudo 2D map via distortion of the energy landscape #### */
  /* ######################################################################### */

  /* first get pair table and loop index of the reference structures */
  pt_ref1       = vrna_pt_get(s1);
  pt_ref2       = vrna_pt_get(s2);
  pt_mfe        = vrna_pt_get(mfe_struct);

  /* then compute maximum matching with disallowed reference structures */
  mm1 = maximumMatchingConstraint(s, pt_ref1);
  mm2 = maximumMatchingConstraint(s, pt_ref2);

  /* get number of bp in reference structures */
  for(bp_ref1 = bp_ref2 = bp_mfe = 0, i = 1; i < length; i++){
    if(pt_ref1[i] > i) bp_ref1++;
    if(pt_ref2[i] > i) bp_ref2++;
    if(pt_mfe[i] > i) bp_mfe++;
  }

  /* compute maximum d1 and maximum d2 */
  MAX_k = bp_ref1 + mm1[idx_1n];
  MAX_l = bp_ref2 + mm2[idx_1n];

  /* get base pair distance between both references */
  bp_dist           = vrna_bp_distance(s1, s2);
  bp_dist_mfe_ref1  = vrna_bp_distance(mfe_struct, s1);
  bp_dist_mfe_ref2  = vrna_bp_distance(mfe_struct, s2);

  /* compute the weight that distorts the landscape such that pr(s1) = pr(s2) = pr(mfe) */
  exp_ref1  = exp(-e_ref1/kT);
  exp_ref2  = exp(-e_ref2/kT);

  if(mmfe == e_ref1){
/*
    distortion_x = 0;
    distortion_y = distortion_x - (e_ref1 - e_ref2) / bp_dist;
    we use the Wolfram alpha solution below ;)
*/
    distortion_x = 0;
    distortion_y = (distortion_x * bp_dist_mfe_ref2 - mmfe + e_ref2) / bp_dist_mfe_ref2;
  } else if(mmfe == e_ref2){
/*
    distortion_y = 0;
    distortion_x = (e_ref1 - e_ref2) / bp_dist + distortion_y;
    we use the Wolfram alpha solution below ;)
*/
    distortion_y = 0;
    distortion_x = (distortion_y * bp_dist_mfe_ref1 + e_ref1 - mmfe) / bp_dist_mfe_ref1;
  } else {
/*
    distortion_x = ((e_ref1 * bp_dist_mfe_ref2) / (bp_dist * bp_dist_mfe_ref1)) - (mmfe / bp_dist_mfe_ref1);
    distortion_y = ((e_ref2 * bp_dist_mfe_ref1) / (bp_dist * bp_dist_mfe_ref2)) - (mmfe / bp_dist_mfe_ref2);
    we use the Wolfram alpha solution below ;)
*/
    double nenner = bp_dist * (bp_dist_mfe_ref1 + bp_dist_mfe_ref2 - bp_dist);
    distortion_x = ( bp_dist_mfe_ref2 * e_ref1 - bp_dist_mfe_ref2 * e_ref2 - bp_dist * mmfe + bp_dist * e_ref2) / nenner;
    distortion_y = (-bp_dist_mfe_ref1 * e_ref1 + bp_dist_mfe_ref1 * e_ref2 - bp_dist * mmfe + bp_dist * e_ref1) / nenner;
  }

  printf("d_x = %1.10f, d_y = %1.10f\n", distortion_x, distortion_y);

  /* create the 2D landscape data structure */
  gridpointT **landscape;
  landscape = (gridpointT **)vrna_alloc(sizeof(gridpointT) * (MAX_k + 1));
  for(i = 0; i <= MAX_k; i++)
    landscape[i] = (gridpointT *)vrna_alloc(sizeof(gridpointT) * (MAX_l + 1));

  /* alloc memory for 1000 structures per partition in the landscape */
  for(i = 0; i <= MAX_k; i++)
    for(j = 0; j <= MAX_l; j++){
      landscape[i][j].max_structs = 1000;
      landscape[i][j].structures  = (char **)vrna_alloc(sizeof(char *) * landscape[i][j].max_structs);
    }

  /* prepare pf fold */
  rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y);
  vrna_rescale_pf_params(vc, &rescale);

  /* apply distortion soft constraints */
  vrna_sc_init(vc);
  kl_soft_constraints *data = kl_init_datastructures( vc,
                                                      (const char *)s1,
                                                      (const char *)s2,
                                                      distortion_x,
                                                      distortion_y);
//  vrna_sc_add_f(vc, &kl_pseudo_energy, NULL);
  vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy,(void *)data);

  /* compute distorted partition function */

  (void)vrna_pf_fold(vc, NULL);

  /* prefill the landscape with 1000 sample structures according to current distortion values */
  for(i = 0; i < maxIterations; i++){
    char *sample  = vrna_pbacktrack(vc);
    /* get k,l coords and free energy */
    int k     = vrna_bp_distance(s1, sample);
    int l     = vrna_bp_distance(s2, sample);
    double fe = (double)vrna_eval_structure(vc, sample);

    /* check if we have sufficient memory allocated and alloc more if necessary */
    if(landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs){
      landscape[k][l].max_structs *= 2;
      landscape[k][l].structures = (char **)vrna_realloc(landscape[k][l].structures, sizeof(char *) * landscape[k][l].max_structs);
    }

    /* insert structure */
    landscape[k][l].structures[landscape[k][l].num_structs] = sample;
    landscape[k][l].num_structs++;
    if(landscape[k][l].mfe > fe)
      landscape[k][l].mfe = fe;

  }

#if 1
  /* now we try to fill the landscape with more values */
  if(0)
  for(j = 1; j <= MAX_k; j++){
    distortion_x -= .2;
    distortion_y -= .2;
    printf("d_x = %1.10f, d_y = %1.10f\n", distortion_x, distortion_y);
    
    /* prepare pf fold */
    rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y);
    vrna_rescale_pf_params(vc, &rescale);

    printf("pf_scale = %6.10f (%6.10f)\n", vc->exp_params->pf_scale, rescale);

    /* apply distortion soft constraints */
    vrna_sc_remove(vc);
    kl_soft_constraints *data = kl_init_datastructures( vc,
                                                      (const char *)s1,
                                                      (const char *)s2,
                                                      distortion_x,
                                                      distortion_y);
    vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy,(void *)data);

    (void) vrna_pf_fold(vc, NULL);

    for(i = 0; i < maxIterations; i++){
      char *sample  = vrna_pbacktrack(vc);
      /* get k,l coords and free energy */
      int k     = vrna_bp_distance(s1, sample);
      int l     = vrna_bp_distance(s2, sample);
      double fe = (double)vrna_eval_structure(vc, sample);

      /* check if we have sufficient memory allocated and alloc more if necessary */
      if(landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs){
        landscape[k][l].max_structs *= 2;
        landscape[k][l].structures = (char **)vrna_realloc(landscape[k][l].structures, sizeof(char *) * landscape[k][l].max_structs);
      }

      /* insert structure */
      landscape[k][l].structures[landscape[k][l].num_structs] = sample;
      landscape[k][l].num_structs++;
      if(landscape[k][l].mfe > fe)
        landscape[k][l].mfe = fe;

    }

  }
#else
  /* go along the direct path diagonal */
  int k_new = bp_dist;
  int l_new = 0;
  for(i = 1; i < bp_dist; i++){
    double fe_diag = 0;
    k_new--;
    l_new++;
    if(landscape[k_new][l_new].num_structs > 0){
      fe_diag = landscape[k_new][l_new].mfe;
    } else {
      /* if we haven't seen any structure here before, we just generate one out of the structures we know */
      for(j = 0; j < landscape[k_new+1][l_new-1].num_structs; j++){
        char  *tmp_struct   = landscape[k_new+1][l_new-1].structures[j];
        short *pt_tmp       = make_pair_table(tmp_struct);
        int   *loopidx_tmp  = pair_table_to_loop_index(pt_tmp);
        /* remove a base pair in the current structure such that we get closer to ref1 and further away from ref2 */
        for(p = 1; p < length; p++){
          if(pt_ref2[p] < p) continue;
          if(pt_ref2[p] == pt_tmp[p]){
            char *tmp2_struct = strdup(tmp_struct);
            tmp2_struct[p-1] = tmp2_struct[pt_tmp[p]-1] = '.';
            double tmp_fe = (double)energy_of_struct_par(s, tmp2_struct, mfe_parameters, 0);
            if(tmp_fe < fe_diag) fe_diag = tmp_fe;
            free(tmp2_struct);
          }
        }

        /* insert a base pair into the current structure such that we get closer to ref1 and further away from ref2 */
        for(p = 1; p < length; p++){
          if(pt_ref1[p] < p) continue;
          if((pt_tmp[p] == 0) && (pt_tmp[pt_ref1[p]] == 0)){
            if(loopidx_tmp[p] == loopidx_tmp[pt_ref1[p]]){
              char *tmp_struct = strdup(tmp_struct);
              tmp2_struct[p-1] = '('; tmp2_struct[pt_ref1[p]-1] = ')';
              double tmp_fe = (double)energy_of_struct_par(s, tmp2_struct, mfe_parameters, 0);
              if(tmp_fe < fe_diag) fe_diag = tmp_fe;
              free(tmp2_struct);
            }
          }
        }

        /* clean up */
        free(pt_tmp);
        free(loopidx_tmp);
      }
    }
    
    /* rescale the distortion values */
    
  }

#endif


  /* debug output */
  for(i = 0; i <= MAX_k; i++)
    for(j = 0; j <= MAX_l; j++){
      if(landscape[i][j].num_structs > 0)
        printf("%d\t%d\t%6.2f\t(%d)\n", i, j, landscape[i][j].mfe, landscape[i][j].num_structs);
    }
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

      case FIND_2D_LANDSCAPE_ESTIMATE:    {
                                            estimate_landscape(vc, s1, s2);
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
