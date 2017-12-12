#ifndef   _RNAXPLORER_REPELLANT_SAMPLING_H_
#define   _RNAXPLORER_REPELLANT_SAMPLING_H_

void
repellant_sampling(vrna_fold_compound_t *fc);

void
rnax_add_repulsion(vrna_fold_compound_t *fc,
                   const char *structure,
                   double     strength);

#endif
