#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "utils.h"
#include "more_utils.h"
#include "pair_mat.h"
#include "energy_const.h"

#include "findpath.h"
#include "meshpoint.h"
#include "RNAwalk.h"
#include "2Dfold.h"
#include "2Dpfold.h"
#include "PathFinder.h"

int         whatToDo = FIND_BEST_FOLDINGPATH;

static char *seq;
static short *S, *S1;
extern int   st_back;

static FLT_OR_DBL *qb_p, *qm_p, *q1k_p, *qln_p, qo_p, *qm1_p, qio_p, qho_p, qmo_p, *qm2_p;
static short *S_p, *S1_p;
static char *ptype_p;


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

void testParallelFold(void);
void transition_rates(const char *seq, const char *s1, const char *s2);
void fake_barriers_output(TwoDpfold_vars *vars);


void usage(void){
  fprintf(stdout, "usage: PathFinder [OPTIONS]\n");
  fprintf(stdout, "\n[OPTIONS]\n");
  fprintf(stdout, "-M [METHOD]  set method\n");
  fprintf(stdout, "   available methods are:\n");
  fprintf(stdout, "   GW        Gradient Walk             (default)\n");
  fprintf(stdout, "   MC        Monte Carlo walk\n");
  fprintf(stdout, "   MC-SA     Monte Carlo Walk\n             with simulated Annealing\n");
  fprintf(stdout, "   DB-MFE    Distance based MFE structure\n             meshpoints");
  fprintf(stdout, "   TRATES    Transition rate computation\n");
  fprintf(stdout, "\n...useful for all methods:\n");
  fprintf(stdout, "-T <float>   temperature\n");
  fprintf(stdout, "-m <float>   maxKeep for direct path search\n");
  fprintf(stdout, "-s <int>     amount of best solutions to hold per iteration\n             (default: 10)\n");
  fprintf(stdout, "-circ        assume RNA molecule to be circular\n");
  fprintf(stdout, "\n...to use with Monte Carlo simulation (MC):\n");
  fprintf(stdout, "-i <int>     number of simulations\n");
  fprintf(stdout, "-r <int>     remember last <int> structure states\n");
  fprintf(stdout, "--penalizeBackWalks\n");
  fprintf(stdout, "\n...to use with simulated annealing:\n");
  fprintf(stdout, "--cooling-rate <float> the cooling factor\n");
  fprintf(stdout, "--tstart       <float> start temperature in deg. Celcius\n");
  fprintf(stdout, "--tstop        <float> stop temperature in deg. Celcius\n");
  fprintf(stdout, "\n...to use with distance based structure meshpoints (DB-MFE):\n");
  fprintf(stdout, "-i <int>     number of iterations (default: 1)\n");
  fprintf(stdout, "-D <int>     max distance exploration for both structures\n             (default: 5)\n");
  fprintf(stdout, "-D1 <int>    max distance exploration for start structure\n");
  fprintf(stdout, "-D2 <int>    max. distance exploration for target structure\n");
  fprintf(stdout, "\n...misc:\n");
  fprintf(stdout, "--basinStructure   just perform a gradient walk starting\n                   at a given structure\n");
  fprintf(stdout, "-h | --help  show this help\n");
  exit(EXIT_SUCCESS);
}

void argument_check(int argc, char *argv[]){
  int i, r;
  circ = 0;
  for(i=1; i<argc; i++){
    if(argv[i][0]=='-')
      switch(argv[i][1]){
        case 'm':   if(argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();
                    r = sscanf(argv[++i], "%d", &maxKeep);
                    if(!r) usage();
                    break;
        case 'i':   if(argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();
                    r = sscanf(argv[++i], "%d", &maxIterations);
                    if(!r) usage();
                    break;
                    
        case 's':   if(argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();
                    r = sscanf(argv[++i], "%d", &maxStorage);
                    if(!r) usage();
                    break;
                    
        case 'r':   if(argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();
                    r = sscanf(argv[++i], "%d", &rememberStructures);
                    if(!r) usage();
                    break;

        case 'M':   if(argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();    
                    if(!strcmp(argv[++i], "MC"))
                      method = MC_METROPOLIS;
                    else if(!strcmp(argv[i], "MC-SA")){
                      method = MC_METROPOLIS;
                      simulatedAnnealing = 1;
                    }
                    else if(!strcmp(argv[i], "GW"))
                      method = GRADIENT_WALK;
                    else if(!strcmp(argv[i], "DB-MFE"))
                      whatToDo = FIND_DISTANCE_BASED_MFE_PATH;
                    else if(!strcmp(argv[i], "BLUBB"))
                      whatToDo = KLKIN;
                    else if(!strcmp(argv[i], "TRATES"))
                      whatToDo = TRANSITION_RATES;
                    break;

        case 'D':   if(i==argc-1) usage();
                    if(argv[i][2]=='\0'){
                      r = sscanf(argv[++i], "%d", &maximum_distance1);
                      if(!r) usage();
                      maximum_distance2 = maximum_distance1;
                    }             
                    else if(argv[i][2]=='1'){
                      r=sscanf(argv[i+1], "%d", &maximum_distance1);
                      if(!r) usage();
                    }
                    else if(argv[i][2]=='2'){
                      r=sscanf(argv[i+1], "%d", &maximum_distance2);
                      if(!r) usage();
                    }
                    else usage();
                    break;

        case 'c':   if(!strcmp(argv[i], "-circ")) circ = 1;
                    else usage();
                    break;

        case '-':   if(!strcmp(argv[i], "--cooling-rate")){
                      r = sscanf(argv[++i], "%f", &treduction);
                    }
                    else if(!strcmp(argv[i], "--tstart")){
                      r = sscanf(argv[++i], "%lf", &tstart);
                      tstart += K0;
                    }
                    else if(!strcmp(argv[i], "--tstop")){
                      r = sscanf(argv[++i], "%lf", &tstop);
                      tstop += K0;
                    }
                    if(!r) usage();
                    
                    if(!strcmp(argv[i], "--penalizeBackWalks")){
                      backWalkPenalty = 1;r=1;
                    }
                    else if(!strcmp(argv[i], "--basinStructure")){
                      whatToDo = FIND_BASIN_STRUCTURE;r=1;
                    }
                    if(!r) usage();
                    break;

        case 'T':   if (argv[i][2]!='\0') usage();
                    if(i==argc-1) usage();
                    r=sscanf(argv[++i], "%lf", &temperature);
                    if (!r) usage();
                    break;
                    
        default:    usage();
    }
  }
}


int main(int argc, char *argv[]) {
  argument_check(argc, argv);
  switch(whatToDo){
    case FIND_BASIN_STRUCTURE:  GetBasinStructure();
                                break;
    
    case 99:                    testParallelFold();
                                break;

    default:                    PathFinder();
                                break;
  }
  return(EXIT_SUCCESS);
}

typedef struct nb_t{
  path_t *path_u;
  path_t *path_d;
  path_t *path_l;
  path_t *path_r;
  path_t *path_lu;
  path_t *path_ld;
  path_t *path_ru;
  path_t *path_rd;

  float barrier_u;
  float barrier_d;
  float barrier_r;
  float barrier_l;
  float barrier_lu;
  float barrier_ld;
  float barrier_ru;
  float barrier_rd;

  struct nb_t *prev;      
  int dir;
    
  float en;
  char *s;
  double pf;
} nb_t;

#define DIR_U   0
#define DIR_D   1
#define DIR_R   2
#define DIR_L   3
#define DIR_LU  4
#define DIR_RU  5
#define DIR_LD  6
#define DIR_RD  7

typedef struct position{
  int i;
  int j;
} position;






typedef struct PQ_t{
  double pf;
  double w;
  int k;
  int l;
} PQ_t;

typedef struct heap{
  PQ_t **heap;
  unsigned int size;
  unsigned int elem;
} heap;

heap *initPQ(unsigned int size);
void PQ_push(heap *h, PQ_t *item);
PQ_t *PQ_pop(heap *h);
void PQ_reheap_down(heap *h, unsigned int pos);
void PQ_reheap_up(heap *h, unsigned int pos);

heap *initPQ(unsigned int size){
  heap *h = (heap *)space(sizeof(heap));
  h->size = ceil(log(size));
  h->elem = 0;
  h->heap = (PQ_t **)space(sizeof(PQ_t *) * h->size);
  return h;
}

void PQ_push(heap *h, PQ_t *item){
  h->heap[h->elem] = item; /* out the item at last position in the heap */
  PQ_reheap_up(h, h->elem++); /* restore heap condition by looking to parental node */
}

PQ_t *PQ_pop(heap *h){
  PQ_t *item = h->heap[0];
  h->heap[0] = h->heap[h->elem-1];
  h->heap[--h->elem] = NULL;
  PQ_reheap_down(h, 0);
  return item;
}

void PQ_reheap_up(heap *h, unsigned int pos){
  unsigned int p = pos;
  unsigned int new_p = p / 2;
  /* parent node of pos is at position pos/2 */
  while (p > 1 && h->heap[new_p]->w < h->heap[p]->w){
    PQ_t *tmp = h->heap[new_p];
    h->heap[new_p] = h->heap[p];
    h->heap[p] = tmp;
    p = p / 2;
    new_p = new_p / 2;
  }
}

void PQ_reheap_down(heap *h, unsigned int pos){
  unsigned int new_p;
  if(2 * pos >= h->elem){
    return;
  }
  if(2 * pos + 1 <= h->elem && h->heap[2*pos+1]->w < h->heap[2*pos]->w){
    new_p = 2 * pos + 1;
  }
  else{
    new_p = 2 * pos;
  }
  if (h->heap[pos]->w > h->heap[new_p]->w){
    PQ_t *tmp = h->heap[pos];
    h->heap[pos] = h->heap[new_p];
    h->heap[new_p] = tmp;
    PQ_reheap_down(h, new_p);
  }
}





/* define a vertex numbering */
#define VN(a,b,dimb)  (a*dimb + b)
#define l_of_VN(pos,dimb)  (pos%dimb)
#define k_of_VN(pos,dimb)  ((pos-(pos%dimb))/dimb)

position *dijkstra(nb_t **weight_matrix, position start, position stop, int dimX, int dimY);

position *dijkstra(nb_t **weight_matrix, position start, position stop, int dimX, int dimY){
  printf("finding best path from (%d,%d) to (%d,%d)\n", start.i, start.j, stop.i, stop.j);
  /* weights are taken from double pf in the matrix */
  int i,j,k,d1,d2, d;
  int dimMax = MAX2(dimX,dimY);

  heap *h = initPQ(dimX * dimY);
  


  /*
  * wt is the vertex indexed array storing  the shortest path from start to each vertex
  * in the beginning, everything is filled with 0.0
  * fortunately as we use partition function as edge weights, we are not seeking for a minimization
  * of the total edge weight sum, but for a maximization, so 0.0 is a good initialization! ;-)
  * weights of 0 are recognized as paths of infinite costs
  */
  double *wt = (double *)space((dimMax * dimMax) * sizeof(double));
  /*
  * st is the vertex indexed array storing the previous node along the shortest path from
  * start... we use this array to obtain the actual path...
  */
  int *st = (int *)space((dimMax * dimMax) * sizeof(int));
  
  /*
  * for all vertices do relaxation as following:
  * v is neighbor of w
  * we visit nodes w in ascending order of the distance to the source node start 
  *
  * if(wt[w] < wt[v] + MAX(pf(v) - pf(w), 0)){
  *   wt[w] = wt[v] + MAX(pf(v) - pf(w), 0);
  *   st[w] = v;
  *
  */
  
  /* the weight from start node to itself is priorised over all others */
  wt[VN(start.i, start.j, dimX)] = 10000.0;
  st[VN(start.i, start.j, dimX)] = VN(start.i, start.j, dimX);
  
  
  
  /* d denotes the distance from the start node along the graph */
  for(d=1;d<=dimMax;d++){
    /* check node (start.i + 2d, start.j) */
    /* for all neighbors of (start.i + 2d, start.j) check weights */
    /* neighbors may only be the nodes (start.i + 2d - 2, start.j), (start.i + 2d - 1, start.j + 1) and (start.i + 2d - 1, start.j - 1) */
    


    
    for(k = d; k > 0; k--){
      /*
      * check nodes (start.i+2*d-k, start.j + k) and
      * (start.i+2*d-k, start.j - k)
      */
      
    }
  }


  free(wt);
  free(st);
}


static void backtrack(nb_t *final, int dir);
void klkin(char *seq,char *s1, char *s2, int maxKeep);

void klkin(char *seq,char *s1, char *s2, int maxkeep){
  short *pt1, *pt2;
  pt1 = make_pair_table(s1);
  pt2 = make_pair_table(s2);
  int i, j, k, n = pt1[0];
  int length = strlen(seq);
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
  int bp_dist = a+b;
  printf("%d:%d ... %d:%d\n", a+c, a+maximum_distance1, b+c, b+maximum_distance2);
  int maxD1 = (a+c < a+maximum_distance1) ? a+maximum_distance1 : a+c;
  int maxD2 = (b+c < b+maximum_distance2) ? b+maximum_distance2 : a+c;
  
  
  TwoDfold_vars *twoD_vars = get_TwoDfold_variables(seq, s1, s2, circ);
  TwoDfold_solution **dfold_structs = (circ) ? TwoDfold_circ_bound(twoD_vars, maxD1, maxD2) : TwoDfold_bound(twoD_vars, maxD1, maxD2);
  int d1,d2;
  float mmfe = INF;
  nb_t **neighbors = (nb_t **)space((maxD1+1) * sizeof(nb_t *));
  int number_of_states = 0;
  float mfe_s1 = dfold_structs[0][bp_dist].en;
  float mfe_s2 = dfold_structs[bp_dist][0].en;
  for(d1 = 0; d1 <= maxD1;d1++){
    neighbors[d1] = (nb_t *)space((maxD2+1) * sizeof(nb_t));
    for(d2 = 0; d2 <= maxD2;d2++){
      neighbors[d1][d2].s = (dfold_structs[d1][d2].en != (float)INF/100.) ? strdup(dfold_structs[d1][d2].s) : NULL;
      neighbors[d1][d2].en = dfold_structs[d1][d2].en;
      neighbors[d1][d2].path_u = NULL;
      neighbors[d1][d2].path_d = NULL;
      neighbors[d1][d2].path_r = NULL;
      neighbors[d1][d2].path_l = NULL;
      neighbors[d1][d2].path_lu = NULL;
      neighbors[d1][d2].path_ld = NULL;
      neighbors[d1][d2].path_ru = NULL;
      neighbors[d1][d2].path_rd = NULL;
      neighbors[d1][d2].barrier_u = (float)INF/100.;
      neighbors[d1][d2].barrier_d = (float)INF/100.;
      neighbors[d1][d2].barrier_r = (float)INF/100.;
      neighbors[d1][d2].barrier_l = (float)INF/100.;
      neighbors[d1][d2].barrier_lu = (float)INF/100.;
      neighbors[d1][d2].barrier_ld = (float)INF/100.;
      neighbors[d1][d2].barrier_ru = (float)INF/100.;
      neighbors[d1][d2].barrier_rd = (float)INF/100.;
      neighbors[d1][d2].prev = NULL;
      neighbors[d1][d2].dir = -1;
      
      
      if(dfold_structs[d1][d2].en != (float)INF/100.){
        number_of_states++;
        mmfe = (mmfe < dfold_structs[d1][d2].en) ? mmfe : dfold_structs[d1][d2].en;
      }
    }
  }
  double kT, sfact=1.07;
  kT = (temperature+K0)*GASCONST/1000.0; /* in Kcal */
  pf_scale = exp(-(sfact*mmfe)/kT/length);
  if (length>2000)
    fprintf(stdout, "scaling factor %f\n", pf_scale);
  TwoDpfold_vars *q_vars = get_TwoDpfold_variables(seq, s1, s2, circ);
  FLT_OR_DBL **pf_s;
  if(circ)
    pf_s = TwoDpfold_circ_bound(q_vars, maxD1, maxD2);
  else
    pf_s = TwoDpfold_bound(q_vars, maxD1, maxD2);

  for(d1 = 0; d1 <= maxD1;d1++){
    for(d2 = 0; d2 <= maxD2;d2++){
      neighbors[d1][d2].pf = pf_s[d1][d2];
    }
  }
  position startPos, stopPos;
  startPos.i = 0;
  startPos.j = bp_dist;
  stopPos.i = bp_dist;
  stopPos.j = 0;
  position *path = dijkstra(neighbors, startPos, stopPos, maxD1+1, maxD2+1);
  
  float **state = (float **)space(number_of_states * sizeof(float *));
  for(i = 0; i<number_of_states; i++) state[i] = (float *)space(number_of_states * sizeof(float));
  float **b_min = (float **)space((maxD1 + 1) * sizeof(float *));
  for(i = 0; i<=maxD1; i++){
    b_min[i] = (float *)space(sizeof(float) * (maxD2 + 1));
    for(j = 0; j <= maxD2; j++){
      b_min[i][j] = (float)INF/100.;
    }
  }
  /* fix point is x_{0,bp_dist(s1,s2)} */
  b_min[0][bp_dist] = 0.;

  printf("k,l neighborhood created...\n");  
  int it;
  for(it = 0; it <= maxIterations; it++)
  for(i=0;i<=maxD1;i++){
    for(k = 0; k<= i+1; k++){
      d2 = b + i + 1;
      if(d2 > maxD2) break;
      if(k > maxD1) break;
      if(!neighbors[k][d2].s) continue;
      /*
      * for all neighbors of x[k, d2]
      * if there is a non-infinite barrier between the fixpoint and our neighbor,
      * we calculate the direct path between the neighbor and our current position
      * to obtain the barrier from the fixpoint to the current position ...
      */
      
#if 0
      /* is the r neighbor inside definition space */
      if(k + 2 <= maxD1 && d2 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[k+2][d2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k+2][d2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            /* as we investigate our right neighbor, the current state is its left neighbor */
            if(!neighbors[k+2][d2].path_l){
              neighbors[k+2][d2].path_l = get_path(seq, neighbors[k+2][d2].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k+2][d2].barrier_l = getSaddlePoint(neighbors[k+2][d2].path_l)->en - neighbors[k+2][d2].en;
              neighbors[k+2][d2].barrier_l = neighbors[k][d2].en - neighbors[k+2][d2].en;
            }
            /* if it seems to be that we have a better barrier if the path is along our right neighbor, */
            /* we just set the pointer prev to it */
            float b_n = MAX2(b_min[k+2][d2], neighbors[k+2][d2].en - mfe_s1 + neighbors[k+2][d2].barrier_l);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k+2][d2];
              neighbors[k][d2].dir = DIR_R;
            }
          }
        }
      }
      /* is the l neighbor inside definition space */
      if(k - 2 >= 0 && d2 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[k-2][d2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k-2][d2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k-2][d2].path_r){
              neighbors[k-2][d2].path_r = get_path(seq, neighbors[k-2][d2].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k-2][d2].barrier_r = getSaddlePoint(neighbors[k-2][d2].path_r)->en - neighbors[k-2][d2].en;
              neighbors[k-2][d2].barrier_r = neighbors[k][d2].en - neighbors[k-2][d2].en;
            }
            float b_n = MAX2(b_min[k-2][d2], neighbors[k-2][d2].en - mfe_s1 + neighbors[k-2][d2].barrier_r);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k-2][d2];
              neighbors[k][d2].dir = DIR_L;
            }
          }
        }
      }
      /* is the u neighbor inside definition space */
      if(k <= maxD1 && d2 + 2<= maxD2){
        /* is this neighbor really existent */
        if(neighbors[k][d2+2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k][d2+2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k][d2+2].path_d){
              neighbors[k][d2+2].path_d = get_path(seq, neighbors[k][d2+2].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k][d2+2].barrier_d = getSaddlePoint(neighbors[k][d2+2].path_d)->en - neighbors[k][d2+2].en;
              neighbors[k][d2+2].barrier_d = neighbors[k][d2].en - neighbors[k][d2+2].en;
            }
            float b_n = MAX2(b_min[k][d2+2], neighbors[k][d2+2].en - mfe_s1 + neighbors[k][d2+2].barrier_d);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k][d2+2];
              neighbors[k][d2].dir = DIR_U;
            }
          }
        }
      }
      /* is the d neighbor inside definition space */
      if(k <= maxD1 && d2 -2 >= 0){
        /* is this neighbor really existent */
        if(neighbors[k][d2-2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k][d2-2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k][d2-2].path_u){
              neighbors[k][d2-2].path_u = get_path(seq, neighbors[k][d2-2].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k][d2-2].barrier_u = getSaddlePoint(neighbors[k][d2-2].path_u)->en - neighbors[k][d2-2].en;
              neighbors[k][d2-2].barrier_u = neighbors[k][d2].en - neighbors[k][d2-2].en;
            }
            float b_n = MAX2(b_min[k][d2-2], neighbors[k][d2-2].en - mfe_s1 + neighbors[k][d2-2].barrier_u);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k][d2-2];
              neighbors[k][d2].dir = DIR_D;
            }
          }
        }
      }
#endif
      /* is the ru neighbor inside definition space */
      if(k + 1 <= maxD1 && d2 + 1 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[k+1][d2+1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k+1][d2+1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k+1][d2+1].path_ld){
              neighbors[k+1][d2+1].path_ld = get_path(seq, neighbors[k+1][d2+1].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k+1][d2+1].barrier_ld = getSaddlePoint(neighbors[k+1][d2+1].path_ld)->en - neighbors[k+1][d2+1].en;
              neighbors[k+1][d2+1].barrier_ld = neighbors[k][d2].en - neighbors[k+1][d2+1].en;
            }
            float b_n = MAX2(b_min[k+1][d2+1], neighbors[k+1][d2+1].en - mfe_s1 + neighbors[k+1][d2+1].barrier_ld);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k+1][d2+1];
              neighbors[k][d2].dir = DIR_RU;
            }
          }
        }
      }
      /* is the rd neighbor inside definition space */
      if(k + 1 <= maxD1 && d2 - 1 >= 0){
        /* is this neighbor really existent */
        if(neighbors[k+1][d2-1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k+1][d2-1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k+1][d2-1].path_lu){
              neighbors[k+1][d2-1].path_lu = get_path(seq, neighbors[k+1][d2-1].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k+1][d2-1].barrier_lu = getSaddlePoint(neighbors[k+1][d2-1].path_lu)->en - neighbors[k+1][d2-1].en;
              neighbors[k+1][d2-1].barrier_lu = neighbors[k][d2].en - neighbors[k+1][d2-1].en;
            }
            float b_n = MAX2(b_min[k+1][d2-1], neighbors[k+1][d2-1].en - mfe_s1 + neighbors[k+1][d2-1].barrier_lu);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k+1][d2-1];
              neighbors[k][d2].dir = DIR_RD;
            }
          }
        }
      }
      /* is the ld neighbor inside definition space */
      if(k - 1 >= 0 && d2 - 1 >= 0){
        /* is this neighbor really existent */
        if(neighbors[k-1][d2-1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k-1][d2-1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k-1][d2-1].path_ru){
              neighbors[k-1][d2-1].path_ru = get_path(seq, neighbors[k-1][d2-1].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k-1][d2-1].barrier_ru = getSaddlePoint(neighbors[k-1][d2-1].path_ru)->en - neighbors[k-1][d2-1].en;
              neighbors[k-1][d2-1].barrier_ru = neighbors[k][d2].en - neighbors[k-1][d2-1].en;
            }
            float b_n = MAX2(b_min[k-1][d2-1], neighbors[k-1][d2-1].en - mfe_s1 + neighbors[k-1][d2-1].barrier_ru);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k-1][d2-1];
              neighbors[k][d2].dir = DIR_LD;
            }
          }
        }
      }
      /* is the lu neighbor inside definition space */
      if(k - 1 >= 0 && d2 + 1 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[k-1][d2+1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[k-1][d2+1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[k-1][d2+1].path_rd){
              neighbors[k-1][d2+1].path_rd = get_path(seq, neighbors[k-1][d2+1].s, neighbors[k][d2].s, maxKeep);
              //neighbors[k-1][d2+1].barrier_rd = getSaddlePoint(neighbors[k-1][d2+1].path_rd)->en - neighbors[k-1][d2+1].en;
              neighbors[k-1][d2+1].barrier_rd = neighbors[k][d2].en - neighbors[k-1][d2+1].en;
            }
            float b_n = MAX2(b_min[k-1][d2+1], neighbors[k-1][d2+1].en - mfe_s1 + neighbors[k-1][d2+1].barrier_rd);
            if(b_min[k][d2] > b_n){
              b_min[k][d2] = b_n;
              neighbors[k][d2].prev = &neighbors[k-1][d2+1];
              neighbors[k][d2].dir = DIR_LU;
            }
          }
        }
      }
    }
    for(j = bp_dist+i; j>=a-i-1; j--){
      d1 = i+1;
      if(j < 0) break;
      if(j > maxD2) continue;
      if(d1 > maxD1) break;
      if(!neighbors[d1][j].s) continue;

      /*
      * for all neighbors of x[d1,j]
      * if there is a non-infinite barrier between the fixpoint and our neighbor,
      * we calculate the direct path between the neighbor and our current position
      * to obtain the barrier from the fixpoint to the current position ...
      */
      
#if 0
      /* is the r neighbor inside definition space */
      if(d1 + 2 <= maxD1 && j <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[d1+2][j].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1+2][j] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            /* as we investigate our right neighbor, the current state is its left neighbor */
            if(!neighbors[d1+2][j].path_l){
              neighbors[d1+2][j].path_l = get_path(seq, neighbors[d1+2][j].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1+2][j].barrier_l = getSaddlePoint(neighbors[d1+2][j].path_l)->en - neighbors[d1+2][j].en;
              neighbors[d1+2][j].barrier_l = neighbors[d1][j].en - neighbors[d1+2][j].en;
            }
            /* if it seems to be that we have a better barrier if the path is along our right neighbor, */
            /* we just set the pointer prev to it */
            float b_n = MAX2(b_min[d1+2][j], neighbors[d1+2][j].en - mfe_s1 + neighbors[d1+2][j].barrier_l);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1+2][j];
              neighbors[d1][j].dir = DIR_R;
            }
          }
        }
      }
      /* is the l neighbor inside definition space */
      if(d1 - 2 >= 0 && j <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[d1-2][j].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1-2][j] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1-2][j].path_r){
              neighbors[d1-2][j].path_r = get_path(seq, neighbors[d1-2][j].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1-2][j].barrier_r = getSaddlePoint(neighbors[d1-2][j].path_r)->en - neighbors[d1-2][j].en;
              neighbors[d1-2][j].barrier_r = neighbors[d1][j].en - neighbors[d1-2][j].en;
            }
            float b_n = MAX2(b_min[d1-2][j], neighbors[d1-2][j].en - mfe_s1 + neighbors[d1-2][j].barrier_r);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1-2][j];
              neighbors[d1][j].dir = DIR_L;
            }
          }
        }
      }
      /* is the u neighbor inside definition space */
      if(d1 <= maxD1 && j + 2<= maxD2){
        /* is this neighbor really existent */
        if(neighbors[d1][j+2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1][j+2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1][j+2].path_d){
              neighbors[d1][j+2].path_d = get_path(seq, neighbors[d1][j+2].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1][j+2].barrier_d = getSaddlePoint(neighbors[d1][j+2].path_d)->en - neighbors[d1][j+2].en;
              neighbors[d1][j+2].barrier_d = neighbors[d1][j].en - neighbors[d1][j+2].en;
            }
            float b_n = MAX2(b_min[d1][j+2], neighbors[d1][j+2].en - mfe_s1 + neighbors[d1][j+2].barrier_d);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1][j+2];
              neighbors[d1][j].dir = DIR_U;
            }
          }
        }
      }
      /* is the d neighbor inside definition space */
      if(d1 <= maxD1 && j -2 >= 0){
        /* is this neighbor really existent */
        if(neighbors[d1][j-2].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1][j-2] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1][j-2].path_u){
              neighbors[d1][j-2].path_u = get_path(seq, neighbors[d1][j-2].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1][j-2].barrier_u = getSaddlePoint(neighbors[d1][j-2].path_u)->en - neighbors[d1][j-2].en;
              neighbors[d1][j-2].barrier_u = neighbors[d1][j].en - neighbors[d1][j-2].en;
            }
            float b_n = MAX2(b_min[d1][j-2], neighbors[d1][j-2].en - mfe_s1 + neighbors[d1][j-2].barrier_u);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1][j-2];
              neighbors[d1][j].dir = DIR_D;
            }
          }
        }
      }
#endif
      /* is the ru neighbor inside definition space */
      if(d1 + 1 <= maxD1 && j + 1 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[d1+1][j+1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1+1][j+1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1+1][j+1].path_ld){
              neighbors[d1+1][j+1].path_ld = get_path(seq, neighbors[d1+1][j+1].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1+1][j+1].barrier_ld = getSaddlePoint(neighbors[d1+1][j+1].path_ld)->en - neighbors[d1+1][j+1].en;
              neighbors[d1+1][j+1].barrier_ld = neighbors[d1][j].en - neighbors[d1+1][j+1].en;
            }
            float b_n = MAX2(b_min[d1+1][j+1], neighbors[d1+1][j+1].en - mfe_s1 + neighbors[d1+1][j+1].barrier_ld);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1+1][j+1];
              neighbors[d1][j].dir = DIR_RU;
            }
          }
        }
      }
      /* is the rd neighbor inside definition space */
      if(d1 + 1 <= maxD1 && j - 1 >= 0){
        /* is this neighbor really existent */
        if(neighbors[d1+1][j-1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1+1][j-1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1+1][j-1].path_lu){
              neighbors[d1+1][j-1].path_lu = get_path(seq, neighbors[d1+1][j-1].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1+1][j-1].barrier_lu = getSaddlePoint(neighbors[d1+1][j-1].path_lu)->en - neighbors[d1+1][j-1].en;
              neighbors[d1+1][j-1].barrier_lu = neighbors[d1][j].en - neighbors[d1+1][j-1].en;
            }
            float b_n = MAX2(b_min[d1+1][j-1], neighbors[d1+1][j-1].en - mfe_s1 + neighbors[d1+1][j-1].barrier_lu);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1+1][j-1];
              neighbors[d1][j].dir = DIR_RD;
            }
          }
        }
      }
      /* is the ld neighbor inside definition space */
      if(d1 - 1 >= 0 && j - 1 >= 0){
        /* is this neighbor really existent */
        if(neighbors[d1-1][j-1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1-1][j-1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1-1][j-1].path_ru){
              neighbors[d1-1][j-1].path_ru = get_path(seq, neighbors[d1-1][j-1].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1-1][j-1].barrier_ru = getSaddlePoint(neighbors[d1-1][j-1].path_ru)->en - neighbors[d1-1][j-1].en;
              neighbors[d1-1][j-1].barrier_ru = neighbors[d1][j].en - neighbors[d1-1][j-1].en;
            }
            float b_n = MAX2(b_min[d1-1][j-1], neighbors[d1-1][j-1].en - mfe_s1 + neighbors[d1-1][j-1].barrier_ru);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1-1][j-1];
              neighbors[d1][j].dir = DIR_LD;
            }
          }
        }
      }
      /* is the lu neighbor inside definition space */
      if(d1 - 1 >= 0 && j + 1 <= maxD2){
        /* is this neighbor really existent */
        if(neighbors[d1-1][j+1].s){
          /* do we already know the barrier between this neighbor and the fixpoint ? */
          if(b_min[d1-1][j+1] != (float)INF/100.){
            /* so we calculate the path and the barrier between this neighbor */
            if(!neighbors[d1-1][j+1].path_rd){
              neighbors[d1-1][j+1].path_rd = get_path(seq, neighbors[d1-1][j+1].s, neighbors[d1][j].s, maxKeep);
              //neighbors[d1-1][j+1].barrier_rd = getSaddlePoint(neighbors[d1-1][j+1].path_rd)->en - neighbors[d1-1][j+1].en;
              neighbors[d1-1][j+1].barrier_rd = neighbors[d1][j].en - neighbors[d1-1][j+1].en;
            }
            float b_n = MAX2(b_min[d1-1][j+1], neighbors[d1-1][j+1].en - mfe_s1 + neighbors[d1-1][j+1].barrier_rd);
            if(b_min[d1][j] > b_n){
              b_min[d1][j] = b_n;
              neighbors[d1][j].prev = &neighbors[d1-1][j+1];
              neighbors[d1][j].dir = DIR_LU;
            }
          }
        }
      }





    }
  }
  
  printf("best barrier along all investigated paths might be: %6.2f (%6.2f)\n", b_min[bp_dist][0], mfe_s1 + b_min[bp_dist][0]); 
  backtrack(&neighbors[bp_dist][0], -1);
  
  destroy_TwoDfold_variables(twoD_vars);
  destroy_TwoDpfold_variables(q_vars);


}

static void backtrack(nb_t *final, int direction){
  if(!final) return;
  if(final->prev)
    backtrack(final->prev, final->dir);

  path_t *r = NULL;
  switch (direction){
#if 0
    case DIR_R:   r = final->path_l;
                  break;
    case DIR_L:   r = final->path_r;
                  break;
    case DIR_D:   r = final->path_u;
                  break;
    case DIR_U:   r = final->path_d;
                  break;
#endif
    case DIR_RU:  r = final->path_ld;
                  break;
    case DIR_RD:  r = final->path_lu;
                  break;
    case DIR_LU:  r = final->path_rd;
                  break;
    case DIR_LD:  r = final->path_ru;
                  break;
  }
  if(r)
    printf("\t%s %6.2f\n", r->s, r->en);
//    for(;r->s;r++)
//      printf("\t%s %6.2f\n", r->s, r->en);

}

void testParallelFold(void){
  temperature = 24.0;
  char *sequence = "AAGCGGCGGCUUUCCGGCUUA";
  char *structure;
  structure = (char *) space(strlen(sequence) + 1);
  float mfe = fold(sequence, structure);
  printf("%s [%6.5f]\n", structure, mfe);

}


void GetBasinStructure(void){
  char *s1;
  seq = get_line(stdin);
  s1 = get_line(stdin);
  initRNAWalk(seq, circ);
  char *basinStructure = structureWalk(seq, s1, GRADIENT_WALK, circ);
  fprintf(stdout, "%s\n", basinStructure);



}

void PathFinder(){
  char *s1 = NULL, *s2=NULL, *line, *start_struct=NULL, *target_struct=NULL;
  int i, istty, n;

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

    seq = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",seq);
    free(line);
    n = (int) strlen(seq);

    s1 = (char *) space((unsigned) n+1);
    s2 = (char *) space((unsigned) n+1);
    if ((start_struct = get_line(stdin))==NULL)
      nrerror("1st structure missing\n");
    else if(strlen(start_struct) != n)
      nrerror("sequence and 1st structure have unequal length");
    strncpy(s1, start_struct, n);

    if ((target_struct = get_line(stdin))==NULL)
      nrerror("2nd structure missing\n");
    else if(strlen(target_struct) != n)
      nrerror("sequence and 2nd et structure have unequal length");
      strncpy(s2, target_struct, n);
    if (istty)
      printf("length = %d\n", n);

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
    do_backtrack = 1;
    dangles = 2;
    double mfe, kT;
    char *ss;
    /* simulate stochastic backtracking in case we need the qm1 array from pf_fold */
    st_back=1;
    ss = (char *) space(strlen(seq)+1);
    mfe = (circ) ? circfold(seq,ss) : fold(seq, ss);
    kT = (temperature+K0)*GASCONST/1000.; /* in Kcal */
    pf_scale = exp(-(1.03*mfe)/kT/n);
    (circ) ? (void) pf_circ_fold(seq,ss) : (void) pf_fold(seq, ss);
    free(ss);
    if(circ){
      if(!get_circ_pf_arrays(&S_p, &S1_p, &ptype_p, &qb_p, &qm_p, &q1k_p, &qln_p, &qm1_p, &qm2_p, &qo_p, &qho_p, &qio_p, &qmo_p)){
        nrerror("wtf circ");
      }
    }
    else
      if(!get_pf_arrays(&S_p, &S1_p, &ptype_p, &qb_p, &qm_p, &q1k_p, &qln_p)){
        nrerror("wtf");
      }

    fprintf(stdout, "%s\n", seq);
    if(whatToDo != TRANSITION_RATES){
      fprintf(stdout, "%s %6.2f\n", s1, (circ) ? energy_of_circ_struct(seq,s1) : energy_of_struct(seq, s1));
      fprintf(stdout, "%s %6.2f\n", s2, (circ) ? energy_of_circ_struct(seq,s2) : energy_of_struct(seq, s2));
    }
    int numSteps;
    path_t *foldingPath, *Saddle, *r;
    curr_iteration = 0;
    switch(whatToDo){
      case KLKIN:                         {
                                            klkin(seq, s1, s2, maxKeep);
                                          }
                                          break;

      case FIND_DISTANCE_BASED_MFE_PATH:  {
                                            int steps, i, j;
                                            foldingPath = get_path(seq, s1, s2, maxKeep/*, &numSteps, circ*/);
                                            Saddle = getSaddlePoint(foldingPath/*, numSteps*/);

                                            fprintf(stdout, "# direct Path:\n# barrier: %6.2f\n\n", Saddle->en);
                                            for(r = foldingPath; r->s; r++){
                                              fprintf(stdout, "%s %6.2f\n", r->s, r->en);
                                            }
                                            free(foldingPath);
                                          
                                            /* this was the old way to compute the folding path... now the new one is following */
                                            fprintf(stdout, "\n# searching for alternative paths\n# ");
                                            fflush(stdout);
                                            foldingPath = levelSaddlePoint2(s1, s2/*, &numSteps*/, 0);
                                            fprintf(stdout, "\n# done\n\n# Path with detours:\n# barrier: %6.2f\n\n", getSaddlePoint(foldingPath/*, numSteps*/)->en);
                                            for(r = foldingPath; r->s; r++){
                                              fprintf(stdout, "%s %6.2f\n", r->s, r->en);
                                            }
                                          }
                                          break;
      case TRANSITION_RATES:              {
                                            transition_rates((const char *)seq, (const char *)s1, (const char *)s2);
                                          }
                                          break;
      default:                            init_rand();
                                          levelSaddlePoint(s1, s2);
                                          break;
    }
  
    free(seq); free(s1); free(s2); free(S); free(S1);free(start_struct);free(target_struct);
  }while(1);
}

void levelSaddlePoint(char *s1, char *s2){

  int iterator = maxIterations;
  int steps;
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
  int steps1, steps2;
  float newSaddleEn = oldSaddleEn;
  char *intermediates[10];
  for(d=0; d<10;d++)
    intermediates[d] = (char *)space((strlen(s1)+1)*sizeof(char));
  float *intermediateEn = (float *)space(sizeof(float));
  int insertedIntermediates = 0;
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
    free(path_left);
    free(path_right);
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


  int iterator = maxIterations;
  int steps, i, j;
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
  pt1 = make_pair_table(s1);
  pt2 = make_pair_table(s2);
  int n = pt1[0];
  
  short *intersect = (short *) space(sizeof(short)*(n+2));
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
    TwoDfold_vars *twoD_vars = get_TwoDfold_variables(seq, s1, s2, circ);
    TwoDfold_solution **dfold_structs = (circ) ? TwoDfold_circ_bound(twoD_vars, a+maximum_distance1, b+maximum_distance2) : TwoDfold_bound(twoD_vars, a+maximum_distance1, b+maximum_distance2);
    destroy_TwoDfold_variables(twoD_vars);

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
          free(path_left);
          free(path_right);
        }
      }
      free(dfold_structs[i]);
    }
  }
  if(bestMeshPoints.count > 0){
    /*
    now as we know n better SaddlePoints, we can iterate deeper to maybe obtain an even better saddle
    */
    
    if(iteration < maxIterations){
      meshpoint *cur, *next;
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
          if(path_left != NULL) free(path_left);
          path_left = t_path_left;
          if(path_right != NULL) free(path_right);
          path_right = t_path_right;
          steps1 = t_steps1;
          steps2 = t_steps2;
          oldSaddleEn = newSaddleEn;
        }
        else{
          free(t_path_left);
          free(t_path_right);
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
    newPath = (path_t *)space((steps1+steps2) * sizeof(path_t));
    memcpy((path_t *)newPath, (path_t *)path_left, steps1*sizeof(path_t));
    memcpy(((path_t *)newPath)+(steps1), ((path_t *)path_right) + 1, (steps2-1)*sizeof(path_t));
  }
  else{
    free(pt1); free(pt2); free(intersect);
    return foldingPath;
  }
  
  free(pt1); free(pt2); free(intersect);
  free(foldingPath);
  if(path_left != NULL) free(path_left);
  if(path_right != NULL) free(path_right);

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


char *pair_table_to_dotbracket(short *pt){
  char *dotbracket = (char *)space((pt[0]+1)*sizeof(char));
  int i;
  for(i=1; i<=pt[0]; i++)
    dotbracket[i-1] = (pt[i]) ? ((i < pt[i]) ? '(' : ')') : '.';
  dotbracket[i-1] = '\0';
  return dotbracket;
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

typedef struct transition_rates_kl{
  double p_nw;
  double p_sw;
  double p_ne;
  double p_se;
  double p_self;
} transition_rates_kl;


typedef struct kl_basin{
  int     k;
  int     l;
  char    *s;
  float   en;
  double  pf;
} kl_basin;

int sort_neighborhood_by_energy_asc(const void *p1, const void *p2){
  if(((kl_basin *)p1)->en > ((kl_basin *)p2)->en) return 1;
  else if(((kl_basin *)p1)->en < ((kl_basin *)p2)->en) return -1;
  return 0;
}


#define KLINDX(a,b,maxa)   (((b) * (maxa+1)) + a)

void transition_rates(const char *s, const char *s1, const char *s2){
  TwoDfold_vars *twoD_vars    = get_TwoDfold_variables(s, s1, s2, 0);
  TwoDfold_solution **kl_mfes = TwoDfold_bound(twoD_vars, -1, -1);
  TwoDpfold_vars  *kl_pf_vars = get_TwoDpfold_variables_from_MFE(twoD_vars);
  FLT_OR_DBL      **kl_pfs    = TwoDpfold_bound(kl_pf_vars, -1, -1);
  init_rand();
  
  int k, l, i, j, iteration;
  short         *pt_ref1      = twoD_vars->reference_pt1;
  short         *pt_ref2      = twoD_vars->reference_pt2;
  int           *loopidx_ref1 = pair_table_to_loop_index(pt_ref1);
  int           *loopidx_ref2 = pair_table_to_loop_index(pt_ref2);
  int           length        = twoD_vars->seq_length;
  char          *ptype        = twoD_vars->ptype;
  unsigned int  *my_iindx     = twoD_vars->my_iindx;
  double kT                   = (temperature+K0)*GASCONST/1000.;   /* kT in kcal/mol  */
  unsigned int  *assignment   = (unsigned int *)space(sizeof(unsigned int) * ((twoD_vars->maxD1 + 1)*(twoD_vars->maxD2 + 1)));
  unsigned int  macro_states  = 0;
  FLT_OR_DBL    pf_scale      = kl_pf_vars->pf_scale;
  int           *k_min_values   = kl_pf_vars->k_min_values;
  int           *k_max_values   = kl_pf_vars->k_max_values;
  int           **l_min_values  = kl_pf_vars->l_min_values;
  int           **l_max_values  = kl_pf_vars->l_max_values;
  int           idx_1n          = my_iindx[1] - length;
  
  /* convert mfe and partition function output to a datastructure that we may easily use for all further analysis */
  kl_basin *basin_data = (kl_basin *)space(sizeof(kl_basin) * ((twoD_vars->maxD1 + 1)*(twoD_vars->maxD2 + 1)));
  double Q = 0;
  for(k = 0; k<=twoD_vars->maxD1; k++){
    for(l = 0; l<=twoD_vars->maxD2; l++){
      if(kl_mfes[k][l].s){
        basin_data[macro_states].k  = k;
        basin_data[macro_states].l  = l;
        basin_data[macro_states].s  = strdup(kl_mfes[k][l].s);
        basin_data[macro_states].en = kl_mfes[k][l].en;
        basin_data[macro_states].pf = (double)kl_pfs[k][l];
        Q                           += (double)kl_pfs[k][l];
        macro_states++;
      }
    }
  }
  /* shorten allocated memory to fit actual data */
  basin_data = (kl_basin *)realloc(basin_data, macro_states*sizeof(kl_basin));
  /* sort neighborhoods according their mfe */
  qsort(basin_data, macro_states, sizeof(kl_basin), sort_neighborhood_by_energy_asc);

  /* go through the macrostate list and extract 'k,l to index' information */
  for(k = 0; k<macro_states;k++){
    assignment[KLINDX(basin_data[k].k, basin_data[k].l, twoD_vars->maxD1)] = k;
  }

  /* get memory for rate matrix */
  double **rate_matrix2 = (double **)space(sizeof(double *) * macro_states);
  for(i = 0; i < macro_states; i++){
    rate_matrix2[i] = (double *)space(sizeof(double) * macro_states);
  }
  
  
  
  for(k = 0; k<macro_states;k++){
    //printf("macrostate %d [%d,%d] %d\n", k, basin_data[k].k, basin_data[k].l, assignment[KLINDX(basin_data[k].k, basin_data[k].l, twoD_vars->maxD1)]);
    
    int idx_self, idx_ne, idx_se, idx_sw, idx_nw;
    int avail_ne, avail_se, avail_sw, avail_nw;
    avail_ne = avail_se = avail_sw = avail_nw = 0;
    
    if((basin_data[k].k+1 >= k_min_values[idx_1n]) && (basin_data[k].k+1 <= k_max_values[idx_1n])){
      if((basin_data[k].l+1 >= l_min_values[idx_1n][basin_data[k].k+1]) && (basin_data[k].l+1 <= l_max_values[idx_1n][basin_data[k].k+1]))
        avail_ne = 1;
      if((basin_data[k].l-1 >= l_min_values[idx_1n][basin_data[k].k+1]) && (basin_data[k].l-1 <= l_max_values[idx_1n][basin_data[k].k+1]))
        avail_se = 1;
    }
    if((basin_data[k].k-1 >= k_min_values[idx_1n]) && (basin_data[k].k-1 <= k_max_values[idx_1n])){
      if((basin_data[k].l+1 >= l_min_values[idx_1n][basin_data[k].k-1]) && (basin_data[k].l+1 <= l_max_values[idx_1n][basin_data[k].k-1]))
        avail_nw = 1;
      if((basin_data[k].l-1 >= l_min_values[idx_1n][basin_data[k].k-1]) && (basin_data[k].l-1 <= l_max_values[idx_1n][basin_data[k].k-1]))
        avail_sw = 1;
    }
    
    idx_self  = k;
    idx_ne    = (avail_ne) ? assignment[KLINDX(basin_data[k].k+1, basin_data[k].l+1, twoD_vars->maxD1)] : -1;
    idx_se    = (avail_se) ? assignment[KLINDX(basin_data[k].k+1, basin_data[k].l-1, twoD_vars->maxD1)] : -1;
    idx_sw    = (avail_sw) ? assignment[KLINDX(basin_data[k].k-1, basin_data[k].l-1, twoD_vars->maxD1)] : -1;
    idx_nw    = (avail_nw) ? assignment[KLINDX(basin_data[k].k-1, basin_data[k].l+1, twoD_vars->maxD1)] : -1;

    double    dG_alpha, dG_ne, dG_se, dG_sw, dG_nw;
    dG_alpha  = (-log(basin_data[k].pf)-length*log(pf_scale))*kT;
    dG_ne     = (avail_ne) ? (-log(basin_data[idx_ne].pf)-length*log(pf_scale))*kT : (double)INF/100.;
    dG_se     = (avail_se) ? (-log(basin_data[idx_se].pf)-length*log(pf_scale))*kT : (double)INF/100.;
    dG_sw     = (avail_sw) ? (-log(basin_data[idx_sw].pf)-length*log(pf_scale))*kT : (double)INF/100.;
    dG_nw     = (avail_nw) ? (-log(basin_data[idx_nw].pf)-length*log(pf_scale))*kT : (double)INF/100.;
    /**
    *** rates are calculated as follows
    *** /f[
    *** k_{\alpha \rightarrow \beta} = \frac{1}{N}\sum_{i \in \alpha}^{N} \sum_{j \in \beta \cap \mathcal{N}(i)}\\
    *** /f]
    *** as we only sample the microstates we have to take care about the detailed balanced condition:
    *** /f[
    *** k_{\alpha \rightarrow \beta} \cdot \pi_{\alpha} = k_{\beta \rightarrow \alpha} \cdot \pi_{\beta}
    *** /f]
    *** we do this by updating k_{\beta \rightarrow \alpha} in each step we generate the neighboring structure
    *** /f$j \in \beta /f$ with /f$ j \in \mathcal{N}(i) /f$
    *** then for each microrate /f$k_{i \rightarrow j} /f$ the reverse microrate /f$k_{j \rightarrow i} \cdot \frac{\pi_{j | \beta}{i | \alpha} /f$
    *** is added to /f$k_{j \rightarrow i} /f$
    *** at the end of all rate calculations, the rate /f$q_{i,i} = -\sum_{i \neq j} q_{i,j} /f$ is calculated
    **/
    for(iteration = 1; iteration <= maxIterations; iteration++){
      //printf("iteration %d\n", iteration);
      unsigned int  n_nw, n_ne, n_se, n_sw;
      double        p_nw, p_ne, p_se, p_sw;
      double        pr_nw, pr_se, pr_sw, pr_ne;
      p_nw = p_ne = p_se = p_sw = pr_nw = pr_ne = pr_se = pr_sw = 0.;
      n_nw = n_ne = n_se = n_sw = 0;

      /* sample a structure */
      char  *s_kl       = TwoDpfold_pbacktrack(kl_pf_vars, basin_data[k].k, basin_data[k].l);
      short *pt_kl      = make_pair_table(s_kl);
      int   *loopidx_kl = pair_table_to_loop_index(pt_kl);
      float en_kl       = energy_of_struct(s, s_kl);
      
      double pr_alpha   = exp((dG_alpha-en_kl)/kT);

      /* generate 'k +- 1', 'l +- 1' - neighbors */
      for(i = 1; i<length; i++){
        /**
        **********************************
        **  remove base pair
        **********************************
        **/
        if(i<pt_kl[i]){
          int delta1, delta2;
          delta1 = delta2 = 0;
          char  *s_neighbor = strdup(s_kl);
          s_neighbor[i-1]   = s_neighbor[pt_kl[i]-1] = '.';
          float en_neighbor = energy_of_struct(s, s_neighbor);
          /* find out which neighbor we've just generated */
          delta1 = (pt_ref1[i] != 0 && (pt_ref1[i] == pt_kl[i])) ? 1 : -1;
          delta2 = (pt_ref2[i] != 0 && (pt_ref2[i] == pt_kl[i])) ? 1 : -1;
          double p_transition     = MIN2(1., exp(-(en_neighbor-en_kl)/kT));
          double p_rev_transition = MIN2(1., exp(-(en_kl-en_neighbor)/kT));
          if((delta1 > 0) && (delta2 > 0))      { p_ne += p_transition; n_ne++; pr_ne += p_rev_transition * exp((dG_ne-en_neighbor)/kT) / pr_alpha;}
          else if((delta1 > 0) && (delta2 < 0)) { p_se += p_transition; n_se++; pr_se += p_rev_transition * exp((dG_se-en_neighbor)/kT) / pr_alpha;}
          else if((delta1 < 0) && (delta2 > 0)) { p_nw += p_transition; n_nw++; pr_nw += p_rev_transition * exp((dG_nw-en_neighbor)/kT) / pr_alpha;}
          else                                  { p_sw += p_transition; n_sw++; pr_sw += p_rev_transition * exp((dG_sw-en_neighbor)/kT) / pr_alpha;}
          free(s_neighbor);
        }
        /**
        **********************************
        **  insert base pair
        **********************************
        **/
        else if(pt_kl[i] == 0){
          for(j = i+TURN+1; j<=length; j++){
            if(pt_kl[j] != 0) continue;
            if(loopidx_kl[i] != loopidx_kl[j]) continue;
            if(ptype[my_iindx[i]-j]){
              int delta1, delta2;
              char *s_neighbor = strdup(s_kl);
              s_neighbor[i-1] = '(';
              s_neighbor[j-1] = ')';
              float en_neighbor = energy_of_struct(s, s_neighbor);
              /* find out which neighbor we've just generated */
              delta1 = (pt_ref1[i] == j) ? -1 : 1;
              delta2 = (pt_ref2[i] == j) ? -1 : 1;
              double p_transition = MIN2(1., exp(-(en_neighbor-en_kl)/kT));
              double p_rev_transition = MIN2(1., exp(-(en_kl-en_neighbor)/kT));
              if((delta1 > 0) && (delta2 > 0))      { p_ne += p_transition; n_ne++; pr_ne += p_rev_transition * exp((dG_ne-en_neighbor)/kT) / pr_alpha;}
              else if((delta1 > 0) && (delta2 < 0)) { p_se += p_transition; n_se++; pr_se += p_rev_transition * exp((dG_se-en_neighbor)/kT) / pr_alpha;}
              else if((delta1 < 0) && (delta2 > 0)) { p_nw += p_transition; n_nw++; pr_nw += p_rev_transition * exp((dG_nw-en_neighbor)/kT) / pr_alpha;}
              else                                  { p_sw += p_transition; n_sw++; pr_sw += p_rev_transition * exp((dG_sw-en_neighbor)/kT) / pr_alpha;}
              free(s_neighbor);
            }
          }
        }
      } /* end generating neighbors */

      if(avail_ne) rate_matrix2[idx_self][idx_ne] += p_ne/maxIterations;
      if(avail_se) rate_matrix2[idx_self][idx_se] += p_se/maxIterations;
      if(avail_nw) rate_matrix2[idx_self][idx_nw] += p_nw/maxIterations;
      if(avail_sw) rate_matrix2[idx_self][idx_sw] += p_sw/maxIterations;
      
      /* add rates for reverse fluxes from neighbors */
      if(avail_ne) rate_matrix2[idx_ne][idx_self] += pr_ne/maxIterations;
      if(avail_se) rate_matrix2[idx_se][idx_self] += pr_se/maxIterations;
      if(avail_nw) rate_matrix2[idx_nw][idx_self] += pr_nw/maxIterations;
      if(avail_sw) rate_matrix2[idx_sw][idx_self] += pr_sw/maxIterations;
      
      /* free allocated memory */
      free(s_kl);
      free(pt_kl);
      free(loopidx_kl);
    }/* end of iteration */
  } /* end of macrostates */

  /* calculate the rate /f$q_{i,i} = -\sum_{i \neq j} q_{i,j} /f$ */
  for(i = 0; i < macro_states; i++){
    double qii = 0;
    for(j = 0; j < macro_states; j++){
      qii -= rate_matrix2[i][j];
    }
    rate_matrix2[i][i] = qii;
  }

  FILE *outfp = fopen("rates.out", "w");
  if (!outfp){
    fprintf(stderr, "can't open file %s\n", "rates.out");
    exit(1);
  }
  for(i = 0; i < macro_states; i++){
    
    printf("%4d %s %6.2f 0 0 0 0 %6.2f 0 %6.2f [%3d,%3d]\n",
            i+1,
            basin_data[i].s,
            basin_data[i].en,
            (-log(basin_data[i].pf)-length*log(kl_pf_vars->pf_scale)) * kT,
            (-log(basin_data[i].pf)-length*log(kl_pf_vars->pf_scale)) * kT,
            basin_data[i].k,
            basin_data[i].l
          );
    for(j = 0; j < macro_states; j++){
      fprintf(outfp, " %10.4g ", rate_matrix2[i][j]);
    }
    fprintf(outfp, "\n");
  }
  fclose(outfp);

#if 0
  /* check detailed balance */
  for(i = 0; i < macro_states; i++){
    for(j = 0; j < macro_states; j++){
      double bla1 = rate_matrix2[i][j] * basin_data[i].pf/Q;
      double bla2 = rate_matrix2[j][i] * basin_data[j].pf/Q;
      if(abs(bla1 - bla2) > 0){
        printf("[%d,%d]\n", i, j);
        printf("k_a->b * pi_a = %6.8f * %6.8f = %6.8f\n", rate_matrix2[i][j], basin_data[i].pf/Q, bla1);
        printf("k_b->a * pi_b = %6.8f * %6.8f = %6.8f\n\n", rate_matrix2[j][i], basin_data[j].pf/Q, bla2);
      }
    }
  }
#endif
 
  for(k = 0; k<=twoD_vars->maxD1; k++){
    for(l = 0; l<=twoD_vars->maxD2; l++) if(kl_mfes[k][l].s) free(kl_mfes[k][l].s);
    free(kl_mfes[k]);
    free(kl_pfs[k]);
  }
  free(kl_mfes);
  free(kl_pfs);

  for(k = 0; k<macro_states;k++){
    free(basin_data[k].s);
  }
  free(basin_data);
  destroy_TwoDfold_variables(twoD_vars);
  destroy_TwoDpfold_variables(kl_pf_vars);
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


void fake_barriers_output(TwoDpfold_vars *vars){

}
