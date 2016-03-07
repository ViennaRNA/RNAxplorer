#ifndef   _RNAXPLORER_BARRIER_LOWER_BOUND_H_
#define   _RNAXPLORER_BARRIER_LOWER_BOUND_H_

#include <ViennaRNA/utils.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/2Dfold.h>

float barrier_estimate_2D(char *seq, vrna_md_t *md_p, char *s1, char *s2, int maximum_distance1, int maximum_distance2);


#endif
