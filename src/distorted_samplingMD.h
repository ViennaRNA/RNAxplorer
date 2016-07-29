#ifndef   _RNAXPLORER_DISTORTED_SAMPLINGMD_H_
#define   _RNAXPLORER_DISTORTED_SAMPLINGMD_H_

//#include <ViennaRNA/data_structures.h>
#include "distorted_sampling.h"

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

typedef struct {
  double kT;
  int *idx;
  int numberOfReferences;
  const char **references;
  double *distortions;
  unsigned int ** referencesBPs;
  short ** referencesAsPairtaibles;
} kl_soft_constraints_MD;

/**
 * Computes the distortions for a given set of structures.
 * The structures should be unique!
 * @param sequence - the rna sequence
 * @param structures - the dot-bracket structures
 * @param numberOfStructures - length of structures = length of the output
 * @return - float array with distortions (one for each structure)
 */
double * rxp_computeDistortions(vrna_fold_compound_t* fc,const char **structures, size_t numberOfStructures);

void print_matrix(char* desc, int m, int n, double* a, int lda);

kl_soft_constraints_MD *kl_init_datastructures_MD(vrna_fold_compound_t *vc, const char **referenceStructures,
    int numberOfReferenceStructures, double *distortions);

void free_kl_soft_constraints_MD(void *data);

FLT_OR_DBL kl_pseudo_energy_MD(int i, int j, int k, int l, char decomp, void *data);

FLT_OR_DBL kl_exp_pseudo_energy_MD(int i, int j, int k, int l, char decomp, void *data);

void fillGridStepwiseBothRef_MD(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int shift,
    int shift_to_first, int verbose, int maxIterations, int maxSteps);

void fillGridStepwiseFirstRef_MD(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax, int verbose,
    int maxIterations, int maxSteps);

void fillGridStepwiseSecondRef_MD(vrna_fold_compound_t *vc, gridLandscapeT *grid, float relaxFactor, int relax,
    int verbose, int maxIterations, int maxSteps);

gridLandscapeT*
estimate_landscapeMD(vrna_fold_compound_t *vc, const char ** refStructures, size_t numberOfReferences,
    int maxIterations, char *extended_options);

/**
 * add the distortion softconstraints (data and functionpointer)to the foldcompound.
 * @param vc - the foldcompound which contains energy parameters (model details).
 * @param structures - the reference structures in dot-bracket notation
 * @param numberOfReferences - size of structures / distortions.
 * @param distortion - pointer to one value which receives the output
 */
void addSoftconstraintsMD(vrna_fold_compound_t *vc, const char ** structures, int numberOfReferences, double * distortions);

#endif
