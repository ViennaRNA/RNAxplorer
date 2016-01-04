#ifndef   _RNAXPLORER_DISTORTED_SAMPLING_H_
#define   _RNAXPLORER_DISTORTED_SAMPLING_H_

#include <ViennaRNA/data_structures.h>

void estimate_landscape( vrna_fold_compound_t *vc,
                    const char *s1,
                    const char *s2,
                    int maxIterations,
                    char *extended_options);

#endif
