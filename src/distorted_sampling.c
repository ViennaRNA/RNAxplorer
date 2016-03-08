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
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/eval.h>
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

kl_soft_constraints *kl_init_datastructures(vrna_fold_compound_t *vc, const char *s1, const char *s2, double x,
		double y) {
	kl_soft_constraints *data;
	unsigned int n;
	char *s = vc->sequence;

	n = strlen(s);

	/* alloc all memory */
	data = (kl_soft_constraints *) vrna_alloc(sizeof(kl_soft_constraints));
	data->kT = vc->exp_params->kT;
	data->idx = vrna_idx_row_wise(n);
	data->ref1 = strdup(s1);
	data->ref2 = strdup(s2);
	data->reference_pt1 = vrna_ptable(data->ref1);
	data->reference_pt2 = vrna_ptable(data->ref2);
	data->referenceBPs1 = vrna_refBPcnt_matrix(data->reference_pt1, TURN); /* matrix containing number of basepairs of reference structure1 in interval [i,j] */
	data->referenceBPs2 = vrna_refBPcnt_matrix(data->reference_pt2, TURN); /* matrix containing number of basepairs of reference structure2 in interval [i,j] */
	data->x = x;
	data->y = y;

	return data;
}

static void free_kl_soft_constraints(void *data) {
	kl_soft_constraints *dat = (kl_soft_constraints *) data;

	free(dat->idx);
	free(dat->ref1);
	free(dat->ref2);
	free(dat->reference_pt1);
	free(dat->reference_pt2);
	free(dat->referenceBPs1);
	free(dat->referenceBPs2);
	free(dat);
}

FLT_OR_DBL kl_pseudo_energy(int i, int j, int k, int l, char decomp, void *data) {

	int d1, d2, ij, kl;
	int *idx = ((kl_soft_constraints *) data)->idx;
	short *reference_pt1 = ((kl_soft_constraints *) data)->reference_pt1;
	short *reference_pt2 = ((kl_soft_constraints *) data)->reference_pt2;
	unsigned int *referenceBPs1 = ((kl_soft_constraints *) data)->referenceBPs1;
	unsigned int *referenceBPs2 = ((kl_soft_constraints *) data)->referenceBPs2;
	double x = ((kl_soft_constraints *) data)->x;
	double y = ((kl_soft_constraints *) data)->y;

	int base_da = (reference_pt1[i] != (unsigned int) j) ? 1 : -1;
	int base_db = (reference_pt2[i] != (unsigned int) j) ? 1 : -1;
	ij = idx[i] - j;
	kl = idx[k] - l;
	d1 = d2 = 0;

	switch (decomp) {
	/* cases where we actually introduce a base pair */

	case VRNA_DECOMP_PAIR_HP:
		d1 = base_da + referenceBPs1[ij];
		d2 = base_db + referenceBPs2[ij];
		break;

	case VRNA_DECOMP_PAIR_IL:
		d1 = base_da + referenceBPs1[ij] - referenceBPs1[kl];
		d2 = base_db + referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_PAIR_ML:
		d1 = base_da + referenceBPs1[ij] - referenceBPs1[kl];
		d2 = base_db + referenceBPs2[ij] - referenceBPs2[kl];
		break;

		/* cases where we split a segment into one or two subsegments */

	case VRNA_DECOMP_ML_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_ML_ML:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_ML_ML_ML:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_ML_UP:
		d1 = referenceBPs1[ij];
		d2 = referenceBPs2[ij];
		break;

	case VRNA_DECOMP_ML_ML_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_ML_COAXIAL: /* (i,j) stacks onto (k,l), lets ignore this case for now */
		d1 = 0;
		d2 = 0;
		break;

	case VRNA_DECOMP_EXT_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_EXT_UP:
		d1 = referenceBPs1[ij];
		d2 = referenceBPs2[ij];
		break;

	case VRNA_DECOMP_EXT_EXT_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_EXT_STEM:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		break;

	case VRNA_DECOMP_EXT_EXT_STEM: /* fall through */
	case VRNA_DECOMP_EXT_STEM_EXT:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - j];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - j];
		break;

	case VRNA_DECOMP_EXT_EXT_STEM1:
		d1 = referenceBPs1[ij] - referenceBPs1[idx[i] - k] - referenceBPs1[idx[l] - (j - 1)];
		d2 = referenceBPs2[ij] - referenceBPs2[idx[i] - k] - referenceBPs2[idx[l] - (j - 1)];
		break;

	case VRNA_DECOMP_EXT_STEM_OUTSIDE:
		d1 = referenceBPs1[ij] - referenceBPs1[kl];
		d2 = referenceBPs2[ij] - referenceBPs2[kl];
		if (k > i) {
			d1 -= referenceBPs1[idx[i] - (k - 1)];
			d2 -= referenceBPs2[idx[i] - (k - 1)];
		}
		if (l < j) {
			d1 -= referenceBPs1[idx[l + 1] - j];
			d2 -= referenceBPs2[idx[l + 1] - j];
		}
		break;

	default:
		d1 = d2 = 0;
		break;
	}

	return (x * d1 + y * d2) * 100;
}

FLT_OR_DBL kl_exp_pseudo_energy(int i, int j, int k, int l, char decomp, void *data) {

	double kT = ((kl_soft_constraints *) data)->kT;
	return exp((-10. * (double) kl_pseudo_energy(i, j, k, l, decomp, data)) / kT);
}

gridLandscapeT*
estimate_landscape(vrna_fold_compound_t *vc, const char *s1, const char *s2, int maxIterations, char *extended_options) {

	unsigned int *mm1, *mm2;
	int i, j, MAX_k, MAX_l, length, *my_iindx;
	int idx_1n;
	int bp_ref1, bp_ref2, bp_mfe, bp_dist, bp_dist_mfe_ref1, bp_dist_mfe_ref2;
	float e_ref1, e_ref2;
	double mmfe, rescale;
	short *pt_ref1, *pt_ref2, *pt_mfe;
	char *s, *mfe_struct;
	double distortion_x, distortion_y;
	kl_soft_constraints *data;
	/* parse extended options string */
	int plain, do_more, both_at_once, relax, verbose, shift, shift_to_first;

	do_more = 0;
	both_at_once = 0;
	relax = 0;
	plain = 0;
	verbose = 0;
	shift = 0;
	shift_to_first = 0;

	if (extended_options) {
		if (strchr(extended_options, 'M')) /* More sampling required */
			do_more = 1;
		if (strchr(extended_options, 'B')) /* alter both potentials at once */
			both_at_once = 1;
		if (strchr(extended_options, 'R')) /* relax potential instead of increasing it */
			relax = 1;
		if (strchr(extended_options, 'S')) /* shift potential to other structure */
			shift = 1;
		if (strchr(extended_options, 'F')) /* shift to first structure */
			shift_to_first = 1;
		if (strchr(extended_options, 'V')) /* verbose */
			verbose = 1;
		if (strchr(extended_options, 'P')) /* plain sampling, no distortion */
			plain = 1;
	}

	s = vc->sequence;
	length = vc->length;
	my_iindx = vc->iindx;
	idx_1n = my_iindx[1] - length;

	/* get mfe for this sequence */
	mfe_struct = (char *) vrna_alloc(sizeof(char) * (length + 1));
	mmfe = (double) vrna_mfe(vc, mfe_struct);

	/* get free energies of the reference structures */
	e_ref1 = vrna_eval_structure(vc, s1);
	e_ref2 = vrna_eval_structure(vc, s2);

	vrna_init_rand();

	/* ######################################################################### */
	/* #### generate a pseudo 2D map via distortion of the energy landscape #### */
	/* ######################################################################### */

	/* first get pair table and loop index of the reference structures */
	pt_ref1 = vrna_ptable(s1);
	pt_ref2 = vrna_ptable(s2);
	pt_mfe = vrna_ptable(mfe_struct);

	/* then compute maximum matching with disallowed reference structures */
	mm1 = maximumMatchingConstraint(s, pt_ref1);
	mm2 = maximumMatchingConstraint(s, pt_ref2);

	/* get number of bp in reference structures */
	for (bp_ref1 = bp_ref2 = bp_mfe = 0, i = 1; i < length; i++) {
		if (pt_ref1[i] > i)
			bp_ref1++;
		if (pt_ref2[i] > i)
			bp_ref2++;
		if (pt_mfe[i] > i)
			bp_mfe++;
	}

	/* compute maximum d1 and maximum d2 */
	MAX_k = bp_ref1 + mm1[idx_1n];
	MAX_l = bp_ref2 + mm2[idx_1n];

	/* get base pair distance between both references */
	bp_dist = vrna_bp_distance(s1, s2);
	bp_dist_mfe_ref1 = vrna_bp_distance(mfe_struct, s1);
	bp_dist_mfe_ref2 = vrna_bp_distance(mfe_struct, s2);

	free(mfe_struct);
	free(mm1);
	free(mm2);
	free(pt_mfe);
	free(pt_ref1);
	free(pt_ref2);

	if (!plain) {
		//if(mmfe == e_ref1)
		if (bp_dist_mfe_ref1 == 0) {
			/*
			 distortion_x = 0;
			 distortion_y = distortion_x - (e_ref1 - e_ref2) / bp_dist;
			 we use the Wolfram alpha solution below ;)
			 */
			distortion_x = 0;
			if (bp_dist_mfe_ref2 != 0) {
				distortion_y = (distortion_x * bp_dist_mfe_ref2 - mmfe + e_ref2) / bp_dist_mfe_ref2;
			}
		}

		//if(mmfe == e_ref2)
		if (bp_dist_mfe_ref2 == 0) {
			/*
			 distortion_y = 0;
			 distortion_x = (e_ref1 - e_ref2) / bp_dist + distortion_y;
			 we use the Wolfram alpha solution below ;)
			 */
			distortion_y = 0;
			if (bp_dist_mfe_ref1 != 0) {
				distortion_x = (distortion_y * bp_dist_mfe_ref1 + e_ref1 - mmfe) / bp_dist_mfe_ref1;
			}
		}

		//GE: what is if s1=s2 and s1!=mfe and s2!=mfe ?

		//else
		if (bp_dist != 0 && (bp_dist_mfe_ref1 + bp_dist_mfe_ref2 != bp_dist)) {
			/*
			 distortion_x = ((e_ref1 * bp_dist_mfe_ref2) / (bp_dist * bp_dist_mfe_ref1)) - (mmfe / bp_dist_mfe_ref1);
			 distortion_y = ((e_ref2 * bp_dist_mfe_ref1) / (bp_dist * bp_dist_mfe_ref2)) - (mmfe / bp_dist_mfe_ref2);
			 we use the Wolfram alpha solution below ;)
			 */
			double nenner = bp_dist * (bp_dist_mfe_ref1 + bp_dist_mfe_ref2 - bp_dist);
			distortion_x = (bp_dist_mfe_ref2 * e_ref1 - bp_dist_mfe_ref2 * e_ref2 - bp_dist * mmfe + bp_dist * e_ref2)
					/ nenner;
			distortion_y = (-bp_dist_mfe_ref1 * e_ref1 + bp_dist_mfe_ref1 * e_ref2 - bp_dist * mmfe + bp_dist * e_ref1)
					/ nenner;
		}

		printf("d_x = %1.10f, d_y = %1.10f\n", distortion_x, distortion_y);

	}

	/* create the 2D landscape data structure */
	gridpointT **landscape;
	landscape = (gridpointT **) vrna_alloc(sizeof(gridpointT) * (MAX_k + 1));
	for (i = 0; i <= MAX_k; i++)
		landscape[i] = (gridpointT *) vrna_alloc(sizeof(gridpointT) * (MAX_l + 1));

	/* alloc memory for 1000 structures per partition in the landscape */
	for (i = 0; i <= MAX_k; i++)
		for (j = 0; j <= MAX_l; j++) {
			landscape[i][j].max_structs = 1000;
			landscape[i][j].structures = (char **) vrna_alloc(sizeof(char *) * landscape[i][j].max_structs);
		}

	/* prepare pf fold */
	if (!plain) {
		rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y);
	} else {
		rescale = mmfe;
	}
	vrna_exp_params_rescale(vc, &rescale);

	if (!plain) {
		/* apply distortion soft constraints */
		vrna_sc_init(vc);
		data = kl_init_datastructures(vc, (const char *) s1, (const char *) s2, distortion_x, distortion_y);

		vrna_sc_add_data(vc, (void *) data, &free_kl_soft_constraints);
		vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy);
	}

	/* compute (distorted) partition function */

	(void) vrna_pf(vc, NULL);

	/* prefill the landscape with 1000 sample structures according to current distortion values */
	for (i = 0; i < maxIterations; i++) {
		char *sample = vrna_pbacktrack(vc);
		/* get k,l coords and free energy */
		int k = vrna_bp_distance(s1, sample);
		int l = vrna_bp_distance(s2, sample);
		double fe = (double) vrna_eval_structure(vc, sample);

		/* check if we have sufficient memory allocated and alloc more if necessary */
		if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
			landscape[k][l].max_structs *= 2;
			landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
					sizeof(char *) * landscape[k][l].max_structs);
		}

		/* insert structure */
		landscape[k][l].structures[landscape[k][l].num_structs] = sample;
		landscape[k][l].num_structs++;
		if (landscape[k][l].mfe > fe)
			landscape[k][l].mfe = fe;

	}

	double bla_x = data->x;
	double bla_y = data->y;

#define RELAX_FACTOR   1.5         /* if set to 1. final relaxation sets x and/or y to 0. */

	if (do_more && (!plain)) {

		if (!both_at_once) {
			/* change potential of first reference */
			if (bla_x > 0)
				for (j = 0; j < bp_dist; j++) {
					if (relax) {
						data->x -= RELAX_FACTOR * bla_x / bp_dist;
					} else {
						data->x += RELAX_FACTOR * bla_x / bp_dist;
					}

					if (verbose)
						fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
					(void) vrna_pf(vc, NULL);
					for (i = 0; i < maxIterations; i++) {
						char *sample = vrna_pbacktrack(vc);
						/* get k,l coords and free energy */
						int k = vrna_bp_distance(s1, sample);
						int l = vrna_bp_distance(s2, sample);
						double fe = (double) vrna_eval_structure(vc, sample);

						/* check if we have sufficient memory allocated and alloc more if necessary */
						if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
							landscape[k][l].max_structs *= 2;
							landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
									sizeof(char *) * landscape[k][l].max_structs);
						}

						/* insert structure */
						landscape[k][l].structures[landscape[k][l].num_structs] = sample;
						landscape[k][l].num_structs++;
						if (landscape[k][l].mfe > fe)
							landscape[k][l].mfe = fe;
					}
				}

			/* change potential of second reference */
			data->x = bla_x;
			if (bla_y > 0.)
				for (j = 0; j < bp_dist; j++) {
					if (relax) {
						data->y -= RELAX_FACTOR * bla_y / bp_dist;
					} else {
						data->y += RELAX_FACTOR * bla_y / bp_dist;
					}
					if (verbose)
						fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
					(void) vrna_pf(vc, NULL);
					for (i = 0; i < maxIterations; i++) {
						char *sample = vrna_pbacktrack(vc);
						/* get k,l coords and free energy */
						int k = vrna_bp_distance(s1, sample);
						int l = vrna_bp_distance(s2, sample);
						double fe = (double) vrna_eval_structure(vc, sample);

						/* check if we have sufficient memory allocated and alloc more if necessary */
						if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
							landscape[k][l].max_structs *= 2;
							landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
									sizeof(char *) * landscape[k][l].max_structs);
						}

						/* insert structure */
						landscape[k][l].structures[landscape[k][l].num_structs] = sample;
						landscape[k][l].num_structs++;
						if (landscape[k][l].mfe > fe)
							landscape[k][l].mfe = fe;
					}
				}
		} else { /* change potential of both references at the same time */
			data->x = bla_x;
			data->y = bla_y;
			for (j = 0; j < bp_dist; j++) {
				if (shift) {
					if (shift_to_first) {
						data->x += RELAX_FACTOR * bla_x / bp_dist;
						data->y -= RELAX_FACTOR * bla_y / bp_dist;
					} else {
						data->x -= RELAX_FACTOR * bla_x / bp_dist;
						data->y += RELAX_FACTOR * bla_y / bp_dist;
					}
				} else {
					if (relax) {
						data->x -= RELAX_FACTOR * bla_x / bp_dist;
						data->y -= RELAX_FACTOR * bla_y / bp_dist;
					} else {
						data->x += RELAX_FACTOR * bla_x / bp_dist;
						data->y += RELAX_FACTOR * bla_y / bp_dist;
					}
				}

				if (verbose)
					fprintf(stderr, "d_x = %1.10f, d_y = %1.10f\n", data->x, data->y);
				(void) vrna_pf(vc, NULL);
				for (i = 0; i < maxIterations; i++) {
					char *sample = vrna_pbacktrack(vc);
					/* get k,l coords and free energy */
					int k = vrna_bp_distance(s1, sample);
					int l = vrna_bp_distance(s2, sample);
					double fe = (double) vrna_eval_structure(vc, sample);

					/* check if we have sufficient memory allocated and alloc more if necessary */
					if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
						landscape[k][l].max_structs *= 2;
						landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
								sizeof(char *) * landscape[k][l].max_structs);
					}

					/* insert structure */
					landscape[k][l].structures[landscape[k][l].num_structs] = sample;
					landscape[k][l].num_structs++;
					if (landscape[k][l].mfe > fe)
						landscape[k][l].mfe = fe;
				}
			}
		}
	}

#if 1
	/* now we try to fill the landscape with more values */
	if (0)
		for (j = 1; j <= MAX_k; j++) {
			distortion_x -= .2;
			distortion_y -= .2;
			printf("d_x = %1.10f, d_y = %1.10f\n", distortion_x, distortion_y);

			/* prepare pf fold */
			rescale = mmfe + (bp_dist_mfe_ref1 * distortion_x) + (bp_dist_mfe_ref2 * distortion_y);
			vrna_exp_params_rescale(vc, &rescale);

			printf("pf_scale = %6.10f (%6.10f)\n", vc->exp_params->pf_scale, rescale);

			/* apply distortion soft constraints */
			vrna_sc_remove(vc);
			kl_soft_constraints *data = kl_init_datastructures(vc, (const char *) s1, (const char *) s2, distortion_x,
					distortion_y);
			vrna_sc_add_data(vc, (void *) data, &free_kl_soft_constraints);
			vrna_sc_add_exp_f(vc, &kl_exp_pseudo_energy);

			(void) vrna_pf(vc, NULL);

			for (i = 0; i < maxIterations; i++) {
				char *sample = vrna_pbacktrack(vc);
				/* get k,l coords and free energy */
				int k = vrna_bp_distance(s1, sample);
				int l = vrna_bp_distance(s2, sample);
				double fe = (double) vrna_eval_structure(vc, sample);

				/* check if we have sufficient memory allocated and alloc more if necessary */
				if (landscape[k][l].num_structs + 2 >= landscape[k][l].max_structs) {
					landscape[k][l].max_structs *= 2;
					landscape[k][l].structures = (char **) vrna_realloc(landscape[k][l].structures,
							sizeof(char *) * landscape[k][l].max_structs);
				}

				/* insert structure */
				landscape[k][l].structures[landscape[k][l].num_structs] = sample;
				landscape[k][l].num_structs++;
				if (landscape[k][l].mfe > fe)
					landscape[k][l].mfe = fe;

			}

		}
#else
	/* go along the direct path diagonal */
	int p;
	int k_new = bp_dist;
	int l_new = 0;
	for(i = 1; i < bp_dist; i++) {
		double fe_diag = 0;
		k_new--;
		l_new++;
		if(landscape[k_new][l_new].num_structs > 0) {
			fe_diag = landscape[k_new][l_new].mfe;
		} else {
			/* if we haven't seen any structure here before, we just generate one out of the structures we know */
			for(j = 0; j < landscape[k_new+1][l_new-1].num_structs; j++) {
				char *tmp_struct = landscape[k_new+1][l_new-1].structures[j];
				short *pt_tmp = vrna_ptable(tmp_struct);
				int *loopidx_tmp = vrna_loopidx_from_ptable(pt_tmp);
				/* remove a base pair in the current structure such that we get closer to ref1 and further away from ref2 */
				for(p = 1; p < length; p++) {
					if(pt_ref2[p] < p) continue;
					if(pt_ref2[p] == pt_tmp[p]) {
						char *tmp2_struct = strdup(tmp_struct);
						tmp2_struct[p-1] = tmp2_struct[pt_tmp[p]-1] = '.';
						double tmp_fe = (double)vrna_eval_structure(vc, tmp2_struct);
						if(tmp_fe < fe_diag) fe_diag = tmp_fe;
						free(tmp2_struct);
					}
				}

				/* insert a base pair into the current structure such that we get closer to ref1 and further away from ref2 */
				for(p = 1; p < length; p++) {
					if(pt_ref1[p] < p) continue;
					if((pt_tmp[p] == 0) && (pt_tmp[pt_ref1[p]] == 0)) {
						if(loopidx_tmp[p] == loopidx_tmp[pt_ref1[p]]) {
							char *tmp_struct = strdup(tmp_struct);
							tmp2_struct[p-1] = '('; tmp2_struct[pt_ref1[p]-1] = ')';
							double tmp_fe = (double)vrna_eval_structure(vc, tmp2_struct);
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

	// evaluate structures, set properties.
	for (i = 0; i <= MAX_k; i++) {
		for (j = 0; j <= MAX_l; j++) {
			if (landscape[i][j].num_structs > 0) {
				int k;
				for (k = 0; k < landscape[i][j].num_structs; k++)
					landscape[i][j].mfe = vrna_eval_structure(vc, landscape[i][j].structures[k]);
				landscape[i][j].k = i;
				landscape[i][j].l = j;
			}
		}
	}
	gridLandscapeT* grid = malloc(sizeof(gridLandscapeT));
	grid->size1 = MAX_k + 1;
	grid->size2 = MAX_l + 1;
	grid->landscape = landscape;
	printf("%d %d", grid->size1, grid->size2);
	return grid;
}

void printLandscape(gridLandscapeT *grid, vrna_fold_compound_t *vc) {
	for (int i = 0; i < grid->size1; i++) {
		for (int j = 0; j < grid->size2; j++) {
			if (grid->landscape[i][j].num_structs > 0) {
				int k;
				for (k = 0; k < grid->landscape[i][j].num_structs; k++)
					printf("%d\t%d\t%6.2f\t(%d) %s\n", i, j, vrna_eval_structure(vc, grid->landscape[i][j].structures[k]),
							grid->landscape[i][j].num_structs, grid->landscape[i][j].structures[k]);
			}
		}
	}
}

void gridLandscape_free(gridLandscapeT *grid) {
	for (int i = 0; i < grid->size1; i++) {
		for (int j = 0; j < grid->size2; j++) {

			for (int k = 0; k < grid->landscape[i][j].num_structs; k++) {
				free(grid->landscape[i][j].structures[k]);
			}
			free(grid->landscape[i][j].structures);
		}
		free(grid->landscape[i]);
	}
	free(grid->landscape);
	free(grid);
}

