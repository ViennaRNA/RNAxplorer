#ifndef   _RNAXPLORER_DISTORTED_SAMPLING_H_
#define   _RNAXPLORER_DISTORTED_SAMPLING_H_

#include <ViennaRNA/data_structures.h>


typedef struct {
  int k;
  int l;
  int num_structs;
  int max_structs;
  char **structures;
  double pf;
  double mfe;
} gridpointT;

gridpointT** estimate_landscape( vrna_fold_compound_t *vc,
                    const char *s1,
                    const char *s2,
                    int maxIterations,
                    char *extended_options);

#endif
