/*
 * paths.cpp
 *
 *  Created on: Mar 7, 2016
 *      Author: entzian
 */

#include "paths.h"
#include <ViennaRNA/2Dfold.h>


/**
 * Computes the direct folding path between the given structures. This path is used as template for an alternative path search.
 * The alternative path consists of meshpoints. The first meshpoint is the saddlepoint on the direct path. All other meshpoints are
 * computed depending on the number of maxIterations.
 * @param seq - the RNA sequence (ACGU)
 * @param s1 - first structure in dot-bracket format.
 * @param s2 - second structure in dot-bracket format.
 * @param maxIterations - maximal iterations of structurewalks.
 * @param maxKeep - how many structures are being kept (get_path in findpath.c)
 * @param method - gradient or random walk (GRADIENT_WALK, MC_METROPOLIS)
 * @param maxStorage - maximal number of meshpoints (see insert_meshpoint in meshpoint.c)
 */
void levelSaddlePoint(const char *seq, const char *s1, const char *s2, int maxIterations, int maxKeep, int method, int maxStorage){

  int iterator = maxIterations;
  vrna_path_t *foldingPath = get_path(seq, s1, s2, maxKeep/*, &steps, circ*/);
  vrna_path_t *Saddle = getSaddlePoint(foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  fprintf(stdout, "old Path:\nbarrier: %6.2f\n\n", Saddle->en);
  int d;
  for(d=0; foldingPath[d].s; d++){
    fprintf(stdout, "%s %6.2f\n", foldingPath[d].s, foldingPath[d].en);
  }
  vrna_path_t *newLeftSaddle, *newRightSaddle;
  vrna_path_t *path_left, *path_right;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list(&bestMeshPoints);
  fprintf(stdout, "\nsearching for alternative paths...");
  fflush(stdout);


  vrna_md_t md;
  vrna_md_set_default(&md);
  /* set user-defined model details */
  md.circ     = circ;
  md.uniq_ML  = 1;
  initRNAWalk(seq, &md);
  while(iterator > 0){
    char *newSaddle = structureWalk(seq, Saddle->s, method);
    path_left = get_path(seq, s1, newSaddle, maxKeep/*, &steps1, circ*/);
    path_right = get_path(seq, newSaddle, s2, maxKeep/*, &steps2, circ*/);

    newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
    newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
    newSaddleEn = MAX2(newLeftSaddle->en,  newRightSaddle->en);
    iterator--;

    if(newSaddleEn < oldSaddleEn){
    	//TODO: -GE: is newSaddle the correct meshpoint? Or should it be newLeftSaddle or newRightSaddle?
      insert_meshpoint(newSaddle, newSaddleEn, &bestMeshPoints, maxStorage);
    }
    else{
    	free(newSaddle);
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

  free_path(foldingPath);
  freeRNAWalkArrays();
}

/**
 * Compute a path between structure s1 and s2.
 * @param seq - the RNA sequence (ACGU)
 * @param s1 - first structure in dot-bracket format.
 * @param s2 - second structure in dot-bracket format.
 * @param iteration - the current iteration in a series of recursive calls. It should be initialized with 0.
 * @param maxIterations - maximal iterations of structurewalks.
 * @param maxKeep - how many structures are being kept (get_path in findpath.c)
 * @param maxStorage - maximal number of meshpoints (see insert_meshpoint in meshpoint.c)
 */
vrna_path_t *levelSaddlePoint2(const char *seq, const char *s1, const char *s2/*, int *num_entry*/, int iteration,
		int maxIterations, int maxKeep, int maxStorage, int maximum_distance1, int maximum_distance2){

  int i;
  vrna_path_t *foldingPath = get_path(seq, s1, s2, maxKeep/*, &steps, circ*/);
  vrna_path_t *Saddle = getSaddlePoint(foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  vrna_path_t *newLeftSaddle, *newRightSaddle;
  vrna_path_t *path_left = NULL;
  vrna_path_t *path_right = NULL;
  vrna_path_t *newPath = NULL;
  int steps1, steps2;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list(&bestMeshPoints);

  /* begin the nice pathfinding routine */
  short *pt1, *pt2;
  pt1 = vrna_ptable(s1);
  pt2 = vrna_ptable(s2);
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

    vrna_fold_compound_t *vc = vrna_fold_compound_TwoD(seq, s1, s2, &md, VRNA_OPTION_MFE);
    vrna_sol_TwoD_t *mfe_s = vrna_mfe_TwoD(vc, a+maximum_distance1, b+maximum_distance2);

    vrna_fold_compound_free(vc);

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
      vrna_path_t *t_path_left = NULL, *t_path_right =  NULL;
      path_left = NULL;
      path_right = NULL;
      for(cur = bestMeshPoints.first; cur != NULL; cur = cur->next){
        t_path_left = levelSaddlePoint2(seq, s1, cur->s, /*&t_steps1,*/ iteration+1,
        		maxIterations, maxKeep, maxStorage, maximum_distance1, maximum_distance2);
        newLeftSaddle = getSaddlePoint(t_path_left/*, t_steps1*/);

        t_path_right = levelSaddlePoint2(seq, cur->s, s2, /*&t_steps2,*/ iteration+1,
        		maxIterations, maxKeep, maxStorage, maximum_distance1, maximum_distance2);
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
    newPath = (vrna_path_t *)vrna_alloc((steps1+steps2) * sizeof(vrna_path_t));
    memcpy((vrna_path_t *)newPath, (vrna_path_t *)path_left, steps1*sizeof(vrna_path_t));
    memcpy(((vrna_path_t *)newPath)+(steps1), ((vrna_path_t *)path_right) + 1, (steps2-1)*sizeof(vrna_path_t));
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

  //fprintf(stdout, "\rsearching for alternative paths...(%2.1f %% done %d %d)", (float)(curr_iteration)/(float)pow(2*maxStorage, maxIterations) * 100.,curr_iteration, (int)pow(2*maxStorage, maxIterations)  );
  fprintf(stdout, ".");
  fflush(stdout);


  return newPath;
}

vrna_path_t *getSaddlePoint(vrna_path_t *foldingPath){
  vrna_path_t *r, *saddle;
  saddle = foldingPath;
  for(r = foldingPath; r->s; r++)
    if(saddle->en < r->en)
      saddle = r;
  //if(saddle->en == foldingPath->en) return --r;
  return saddle;
}
