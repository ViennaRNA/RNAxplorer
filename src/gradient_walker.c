/*
 ============================================================================
 Name        : GradientWalker.c
 Author      : GE
 Version     :
 Copyright   : Your copyright notice
 Description : Compute recursively the steepest descent structure.
 ============================================================================
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/io/utils.h>
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/walk.h"
#include "ViennaRNA/neighbor.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/read_epars.h"


PRIVATE void
printStructure_pt(vrna_fold_compound_t *vc, short * pt)
{
  char * structure = vrna_db_from_ptable (pt);
  float energy = vrna_eval_structure_pt(vc,pt) / 100.0f;
  printf ("%s %9.4g\n", structure,energy);
  free (structure);
}

int
gradient_walker(double temperature_celsius, int shift_moves, char *parameter_file, const char *sequence, char **structures)
{
  double temperature = temperature_celsius;
  int shifts = shift_moves;
  if(parameter_file != NULL){
    read_parameter_file(parameter_file);
  }

  int moveset = VRNA_MOVESET_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT;
  if(shifts)
    moveset = VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT | VRNA_PATH_NO_TRANSITION_OUTPUT;


  vrna_md_t md;
  vrna_md_set_default (&md);
  md.dangles = 2;
  md.noLP = 0;
  md.noGU = 0;
  md.temperature = temperature;
  vrna_fold_compound_t *vc = vrna_fold_compound (sequence, &md, VRNA_OPTION_EVAL_ONLY);
  char *structure;
  short *pt;
  int i;
  for(i = 0;structures[i] !=NULL; i++){
    structure = structures[i];
    pt = vrna_ptable (structure);
    vrna_move_t *moves = vrna_path_gradient(vc,pt, moveset);
    free(moves);
    printStructure_pt(vc, pt);
    free (pt);
  }

  //printf ("finish!\n");
  vrna_fold_compound_free (vc);

  return EXIT_SUCCESS;
}
