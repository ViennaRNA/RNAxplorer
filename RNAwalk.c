
/* this must be included first due to the dependency of pair_mat.h from fold_vars.h */
#include "RNAwalk.h"


#include "fold.h"
#include "utils.h"
#include "energy_const.h"
#include "pair_mat.h"

#include "meshpoint.h"



int         simulatedAnnealing = 1;
FLT_OR_DBL  tstart = 37.0 + K0;
FLT_OR_DBL  tstop = 0.0 + K0;
float       treduction = 0.9998;

int         rememberStructures = 10;
int         maxRest = 100;

int         backWalkPenalty = 0;

static short *S = NULL, *S1 = NULL;

void initRNAWalk(char *seq, int circ){
  int i;
  init_rand();
  update_fold_params();
  make_pair_matrix();

  /* nummerically encode sequence */
  S = (short *) space(sizeof(short)*(strlen(seq)+2));
  S1 = (short *) space(sizeof(short)*(strlen(seq)+2));
  S[0] = S1[0] = strlen(seq);
  for (i=0; i< strlen(seq); i++) {
    S[i+1] = encode_char(seq[i]);
    S1[i+1] = alias[S[i+1]];
  }
  if(circ){
    S1[0] = S1[S[0]];
    S1[S[0]+1] = S1[1];
  }
}

void freeRNAWalkArrays(void){
  if(S) free(S);
  if(S1) free(S1);
}

char *structureWalk(char *seq, char *structure, int method, int circ){
  char *minStruct = strdup(structure);
  float curr_en = (circ) ? energy_of_circ_struct(seq,structure) : energy_of_struct(seq, structure);
  short *curr_pt = make_pair_table(structure);

  /* generate all neighbors that have 1 bp more */
  int i, j;
  float min_e = curr_en;
  float bottom_en = min_e;
  int curr_loop;
  FLT_OR_DBL tcurr = tstart;
  structure_queue *lastStructures = NULL;
  if(rememberStructures > 0){
    lastStructures = (structure_queue *)space(sizeof(structure_queue));
    init_structure_queue(lastStructures);  
  }
  while(1){
    short *curr_pt = make_pair_table(minStruct);
    int *loopidx = pair_table_to_loop_index(curr_pt);
    int neighbor_cnt = 0;
    neighbor *neighbors = NULL;
    
    float lowestNeighborEn = INF;
    int lowestNeighborIdx = 0;
    int min_i = 0, min_j = 0;
    min_e = bottom_en;
    meshpoint *state = NULL;
    if((state = is_in_queue(minStruct, lastStructures)) != NULL){
      neighbors = state->neighbor_list;
      neighbor_cnt = state->neighbor_cnt;
      lowestNeighborEn = state->bestNeighborEn;
    }
    else{

      for(i=1; i <= curr_pt[0]; i++){
        curr_loop = loopidx[i];
        if(!curr_pt[i])
          for(j=i+TURN+1; j <= curr_pt[0]; j++){
            if(!curr_pt[j] && (curr_loop == loopidx[j]) && pair[S1[i]][S1[j]]){
              char *s = strdup(minStruct);
              s[i-1] = '(';
              s[j-1] = ')';
              neighbors = (neighbor *)realloc(neighbors, ++neighbor_cnt*sizeof(neighbor));
              neighbors[neighbor_cnt-1].i = i;
              neighbors[neighbor_cnt-1].j = j;
              neighbors[neighbor_cnt-1].en = (circ) ? energy_of_circ_struct(seq, s) : energy_of_struct(seq, s);
              if(lowestNeighborEn > neighbors[neighbor_cnt-1].en){
                lowestNeighborEn = neighbors[neighbor_cnt-1].en;
                lowestNeighborIdx = neighbor_cnt-1;
              }
              free(s);
            }
          }
      }
      /* generate all neighbors that have 1 bp less */
      for(i=1; i <= curr_pt[0]; i++){
        if(i<curr_pt[i]){
          char *s = strdup(minStruct);
          s[i-1] = s[curr_pt[i]-1] = '.';
          neighbors = (neighbor *)realloc(neighbors, ++neighbor_cnt*sizeof(neighbor));
          neighbors[neighbor_cnt-1].i = -i;
          neighbors[neighbor_cnt-1].j = -curr_pt[i];
          neighbors[neighbor_cnt-1].en = (circ) ? energy_of_circ_struct(seq, s) : energy_of_struct(seq, s);
          if(lowestNeighborEn > neighbors[neighbor_cnt-1].en){
            lowestNeighborEn = neighbors[neighbor_cnt-1].en;
            lowestNeighborIdx = neighbor_cnt-1;
          }
          free(s);
        }
      }
      qsort(neighbors, neighbor_cnt, sizeof(neighbor), sort_neighbors_by_energy_asc);
    }
    if(neighbor_cnt > 0) min_e = MIN2(bottom_en, neighbors[0].en);
    else min_e = bottom_en;

    int ThisIsTheEnd = 0;
    if(rememberStructures > 0){
      insert_structure_in_queue(lastStructures, minStruct, bottom_en, neighbors, neighbor_cnt, lowestNeighborEn,lowestNeighborIdx, rememberStructures);
    }
    switch(method){

/* original metrolpolis rule implementation of the monte carlo walk...
*  this is not needed anymore as we can perform a rejection less walk
*  with the implementation below...
*      case METROPOLIS:      {
*                              int found_path = 0;
*                              float kT = (tcurr)*1.98717/1000.;
*                              if(min_i == 0){
*                                float prob = exp((bottom_en-lowestNeighborEn)/kT);
*                                if(urn() < prob){
*                                  ThisIsTheEnd = 1;
*                                  bottom_en = min_e;
*                                  break;
*                                }
*                              }
*                              while(!found_path){
*                                int path = int_urn(0, neighbor_cnt-1);
*                                float deltaG = bottom_en-neighbor_en[path];
*                                float prob = MIN2(1., exp(deltaG/kT));
*                                float pr = urn();
*                                if(pr < prob){
*                                  if(neighbor_i[path] > 0){
*                                    minStruct[neighbor_i[path]-1] = '(';
*                                    minStruct[neighbor_j[path]-1] = ')';
*                                  }
*                                  else{
*                                    minStruct[-neighbor_i[path]-1] = minStruct[-neighbor_j[path]-1] = '.';
*                                  }
*                                  found_path = 1;
*                                  bottom_en = neighbor_en[path];
*                                }
*                              }
*                            }
*                            break;
*/
      case MC_METROPOLIS:  {
                              float kT = (tcurr)*GASCONST/1000.;
                              float penalty = (backWalkPenalty) ? 2*kT : 0.; /* penalty of 2kT for recently visited states */
                              float lambda;
                              char *s = strdup(minStruct);
                              if(neighbors[0].i < 0)
                                s[-neighbors[0].i-1] = s[-neighbors[0].j-1] = '.';
                              else{
                                s[neighbors[0].i-1] = '(';
                                s[neighbors[0].j-1] = ')';
                              }
                              float deltaG = bottom_en-neighbors[0].en + ((is_in_queue(s, lastStructures) != NULL) ? penalty : 0.);
                              float *pathProbs = (float *)space(neighbor_cnt*sizeof(float));
                              pathProbs[0] = (MIN2(1., exp(deltaG/kT)))/neighbor_cnt;
                              free(s);
                              int i;
                              for(i=1; i<neighbor_cnt; i++){
                                s = strdup(minStruct);
                                if(neighbors[i].i < 0)
                                  s[-neighbors[i].i-1] = s[-neighbors[i].j-1] = '.';
                                else{
                                  s[neighbors[i].i-1] = '(';
                                  s[neighbors[i].j-1] = ')';
                                }
                                deltaG = bottom_en-neighbors[i].en + ((is_in_queue(s, lastStructures) != NULL) ? penalty : 0.);
                                pathProbs[i] = pathProbs[i-1] + (MIN2(1., exp(deltaG/kT)))/neighbor_cnt;
                                free(s);
                              }
                              lambda = pathProbs[neighbor_cnt-1];
                              float rnd = urn()*lambda;
                              int found_path = getPosition(pathProbs, rnd, neighbor_cnt-1);
                              
                              if(min_i == 0){
                                float prob = exp((bottom_en-neighbors[0].en)/kT);
                                if(urn() > prob){
                                  ThisIsTheEnd = 1;
                                  bottom_en = min_e;
                                  break;
                                }
                              }
                              /*
                              if((int)(log(urn())/log(1-lambda)) > maxRest){
                                ThisIsTheEnd = 1;
                                bottom_en = min_e;
                                break;
                              }
                              */
                              if(neighbors[found_path].i > 0){
                                minStruct[neighbors[found_path].i-1] = '(';
                                minStruct[neighbors[found_path].j-1] = ')';
                              }
                              else{
                                minStruct[-neighbors[found_path].i-1] = minStruct[-neighbors[found_path].j-1] = '.';
                              }
                              bottom_en = neighbors[found_path].en;
                              free(pathProbs);
                              if(simulatedAnnealing){
                                tcurr *= treduction;
                                if(tcurr<=tstop){
                                  ThisIsTheEnd = 1;
                                  bottom_en = min_e;
                                  break;
                                }
                              }
                            }
                            break;
      /* default is a gradient walk */
      default:              if(bottom_en > min_e){
                              ThisIsTheEnd = 1;
                              break;
                            }
                            if(neighbors[0].i>0){
                              minStruct[neighbors[0].i-1] = '(';
                              minStruct[neighbors[0].j-1] = ')';
                            }
                            else if(neighbors[0].i<0){
                              minStruct[-neighbors[0].i-1] = minStruct[-neighbors[0].j-1] = '.';
                            }
                            else ThisIsTheEnd = 1;
                            bottom_en = min_e;
                            break;
    }                              
    free(curr_pt);
    free(loopidx);
    if(rememberStructures < 1) free(neighbors);
    if(ThisIsTheEnd) break;
  }
  if(rememberStructures > 0){
    clear_structure_queue(lastStructures);  
  }
  return minStruct;
}


int getPosition(float *array, float value, int max_idx){
  if(max_idx == 0) return 0;
  int center_idx = (max_idx%2==1) ? (max_idx-1)/2 :max_idx/2; 
  if(value>array[center_idx])
    return getPosition(array+center_idx+1, value, max_idx-center_idx-1) + center_idx+1;
  else
    return getPosition(array, value, center_idx);
}

/* this comes from findpath.c ... */
static int *pair_table_to_loop_index (short *pt)
{
  /* number each position by which loop it belongs to (positions structure
     at 1) */
  int i,hx,l,nl;
  int length;
  int *stack = NULL;
  int *loop = NULL;

  length = pt[0];
  stack  = (int *) space(sizeof(int)*(length+1));
  loop   = (int *) space(sizeof(int)*(length+2));
  hx=l=nl=0;

  for (i=1; i<=length; i++) {
    if ((pt[i] != 0) && (i < pt[i])) { /* ( */
      nl++; l=nl;
      stack[hx++]=i;
    }
    loop[i]=l;

    if ((pt[i] != 0) && (i > pt[i])) { /* ) */
      --hx;
      if (hx>0)
        l = loop[stack[hx-1]];  /* index of enclosing loop   */
      else l=0;                 /* external loop has index 0 */
      if (hx<0) {
        nrerror("unbalanced brackets in make_pair_table");
      }
    }
  }
  loop[0] = nl;
  free(stack);
  return (loop);
}

