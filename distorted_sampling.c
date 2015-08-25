/*
    Secondary structure landscape sampling via distortion of the
    partition function

    (c) Ronny Lorenz
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "distorted_sampling.h"

/*
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/energy_const.h>
#include <ViennaRNA/findpath.h>
#include <ViennaRNA/2Dpfold.h>
*/
#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mm.h>

typedef struct {
  double kT;
  int *idx;
  char *ref1;
  char *ref2;
  double x;
  double y;
  unsigned int *referenceBPs1;
  unsigned int *referenceBPs2;
  short *reference_pt1;
  short *reference_pt2;
} kl_soft_constraints;


kl_soft_constraints *kl_init_datastructures(vrna_fold_compound *vc, const char *s1, const char *s2, double x, double y){
  kl_soft_constraints *data;
  unsigned int n;
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
                    const char *s2,
                    int maxIterations){

  unsigned int      *mm1, *mm2;
  int               i, j, MAX_k, MAX_l, length, *my_iindx;
  int               idx_1n;
  int               bp_ref1, bp_ref2, bp_mfe, bp_dist, bp_dist_mfe_ref1, bp_dist_mfe_ref2;
  float             e_ref1, e_ref2;
  double            mmfe, rescale;
  short             *pt_ref1, *pt_ref2, *pt_mfe;
  char              *s, *mfe_struct;
  double            distortion_x, distortion_y;

  s         = vc->sequence;
  length    = vc->length;
  my_iindx  = vc->iindx;
  idx_1n    = my_iindx[1] - length;


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
  vrna_exp_params_rescale(vc, &rescale);

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
    vrna_exp_params_rescale(vc, &rescale);

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

