/*
 * paths.h
 *
 *  Created on: Mar 7, 2016
 *      Author: entzian
 */

#ifndef _PATHS_H_
#define _PATHS_H_

#include "RNAwalk.h"
#include <ViennaRNA/findpath.h>


/**
 * !
 * @param seq - the RNA sequence (AGCU).
 * @param s1 - the first secondary structure in dot-bracket format.
 * @param s2 - the second secondary structure in dot-bracket format.
 * @param maxIterations -
 * @param maxKeep - maximum structures that will be kept (see findpath.h)
 * @param method - method for structurewalk (MC_METROPOLIS or GRADIENT_WALK) (see RNAwalk.h)
 * @param maxStorage - for insert_meshpoint
 */
void levelSaddlePoint(const char *seq, const char *s1, const char *s2, int maxIterations, int maxKeep, int method, int maxStorage);

vrna_path_t *levelSaddlePoint2(const char *seq, const char *s1, const char *s2/*, int *num_entry*/, int iteration,
		int maxIterations, int maxKeep, int maxStorage, int maximum_distance1, int maximum_distance2);
vrna_path_t *getSaddlePoint(vrna_path_t *foldingPath);

#endif /* _PATHS_H_ */
