#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/io/file_formats.h>

#include "RNAwalk.h"
#include "meshpoint.h"
#include "barrier_lower_bound.h"
#include "distorted_sampling.h"
#include "distorted_samplingMD.h"
#include "repellant_sampling.h"
#include "paths.h"

#include "RNAxplorer_cmdl.h"

enum strategies_avail_e {
  MOVE_GRADIENT_WALK,
  PATHFINDER_SADDLE_GRADIENT_WALK,  /* perform gradient walks from saddle point to find meshpoint(s) */
  PATHFINDER_SADDLE_MONTE_CARLO,    /* perform rejection-less monte carlo (Gillespie) simulation away from saddle point */
  PATHFINDER_SADDLE_MONTE_CARLO_SA, /* perform rejection-less monte carlo (Gillespie) simulation away from saddle point while cooling down the system (simulated annealing) */
  PATHFINDER_TWO_D_REPRESENTATIVES, /* use 2D representatives as meshpoints */
  FINDPATH,                         /* use findpath heuristic to obtain optimal refolding path */
  TWO_D_LOWER_BOUND,                /* compute optimal folding path based on 2D representatives using Dijkstra algorithm on andscape projection */
  REPELLENT_SAMPLING,
  ATTRACTION_SAMPLING,
  TEMPERATURE_SCALING_SAMPLING
};

struct options_s {
  enum strategies_avail_e strategy;
  vrna_md_t               md;

  /* findpath option(s) */
  int                     max_keep;

  /* common options */
  int                     iterations;
  int                     samples;

  /* PathFinder options */
  int                     max_storage;

  /* 2D fold options */
  int                     max_d1;
  int                     max_d2;

  /* simulated annealing options */
  int                     simulated_annealing;
  float                   t_start;
  float                   t_end;
  float                   cooling_rate;
};

typedef int (xplorer_func)(const char       *rec_id,
                           const char       *orig_sequence,
                           char             **structures,
                           struct options_s *opt);

typedef struct {
  enum strategies_avail_e strategy;
  xplorer_func            *f;
  const char              *name;
} strategies;


struct options_s *
process_arguments(int   argc,
                  char  *argv[]);


char **
extract_structures(unsigned int n,
                   int          maybe_multiline,
                   char         **rec_rest);


int
moves_gradient_descent(const char       *rec_id,
                       const char       *orig_sequence,
                       char             **structures,
                       struct options_s *opt);


int
paths_findpath(const char       *rec_id,
               const char       *orig_sequence,
               char             **structures,
               struct options_s *opt);


int
paths_pathfinder_gd(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_mc(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_mcsa(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt);


int
paths_pathfinder_db(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
paths_pathfinder_dbba(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt);


int
sampling_repulsion(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt);


int
sampling_attraction(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt);


int
sampling_temperature(const char       *rec_id,
                     const char       *orig_sequence,
                     char             **structures,
                     struct options_s *opt);


/**
*** \file RNAxplorer.c
**/


#define NUM_STRATEGIES    10

static strategies known_strategies[NUM_STRATEGIES] = {
  /* code, function, name */
  {
    MOVE_GRADIENT_WALK,
    &moves_gradient_descent,
    "Gradient Descent Moves"
  },
  {
    PATHFINDER_SADDLE_GRADIENT_WALK,
    &paths_pathfinder_gd,
    "PathFinder - Gradient Descent from Saddle"
  },
  {
    PATHFINDER_SADDLE_MONTE_CARLO,
    &paths_pathfinder_mc,
    "PathFinder - MCMC from Saddle"
  },
  {
    PATHFINDER_SADDLE_MONTE_CARLO_SA,
    &paths_pathfinder_mcsa,
    "PathFinder - MCMC from Saddle (Simulated Annealing)"
  },
  {
    PATHFINDER_TWO_D_REPRESENTATIVES,
    &paths_pathfinder_db,
    "PathFinder - 2D Representatives"
  },
  {
    FINDPATH,
    &paths_findpath,
    "Findpath"
  },
  {
    TWO_D_LOWER_BOUND,
    &paths_pathfinder_dbba,
    "Energy Barrier Lower Bound from 2D Representation"
  },
  {
    REPELLENT_SAMPLING,
    &sampling_repulsion,
    "Repellent Sampling Scheme"
  },
  {
    ATTRACTION_SAMPLING,
    &sampling_attraction,
    "Directed Sampling Scheme"
  },
  {
    TEMPERATURE_SCALING_SAMPLING,
    &sampling_temperature,
    "Temperature Scaling Sampling Scheme"
  }
};

static char       *extended_options = NULL;

//percentage of distortion per reference.
size_t            length_indicesAndPercentages  = 0;
double            *indicesAndPercentages        = NULL;


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
int
main(int  argc,
     char *argv[])
{
  FILE              *input_stream     = stdin;
  xplorer_func      *processing_func  = NULL;

  struct options_s  *options = process_arguments(argc, argv);

  for (int i = 0; i < NUM_STRATEGIES; i++)
    if (options->strategy == known_strategies[i].strategy) {
      processing_func = known_strategies[i].f;
      break;
    }

  int           istty_in  = isatty(fileno(input_stream));
  int           istty_out = isatty(fileno(stdout));

  unsigned int  read_opt = 0;

  if (istty_in)
    vrna_message_input_seq("Input sequence (upper or lower case) followed by structures");

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty_in)
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;

  /* initialize random number generator from libRNA */
  vrna_init_rand();


  /* main loop that processes each record obtained from input stream */
  do {
    char          *rec_sequence, *rec_id, **rec_rest;
    unsigned int  rec_type;
    int           maybe_multiline;

    rec_id          = NULL;
    rec_rest        = NULL;
    maybe_multiline = 0;

    rec_type = vrna_file_fasta_read_record(&rec_id,
                                           &rec_sequence,
                                           &rec_rest,
                                           input_stream,
                                           read_opt);

    if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))
      break;

    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    unsigned int  n = strlen(rec_sequence);

    char          **structures = extract_structures(n, maybe_multiline, rec_rest);

    processing_func(rec_id, rec_sequence, structures, options);

    free(rec_sequence);
    free(rec_id);

    /* free the rest of current dataset */
    if (structures) {
      for (int i = 0; structures[i]; i++)
        free(structures[i]);
      free(structures);
    }

    if (istty_in)
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structures");
  } while (1);

  return EXIT_SUCCESS;
}


struct options_s *
default_options(void)
{
  struct options_s *options = (struct options_s *)vrna_alloc(sizeof(struct options_s));

  /* default strategy */
  options->strategy = FINDPATH;

  /* default energy model settings */
  vrna_md_set_default(&(options->md));
  options->md.uniq_ML = 1; /* we certainly require unique multibranch loop decomposition in any case */

  /* findpath option(s) */
  options->max_keep = 10;

  /* common options to many methods */
  options->samples    = 1000;
  options->iterations = 1;

  /* PathFinder options */
  options->max_storage = 10;

  /* 2D fold options */
  options->max_d1 = 5;
  options->max_d2 = 5;

  /* simulated annealing options */
  options->simulated_annealing  = 0;
  options->t_start              = 37.0 + K0;
  options->t_end                = 0. + K0;
  options->cooling_rate         = 0.9998;

  return options;
}


struct options_s *
process_arguments(int   argc,
                  char  *argv[])
{
  struct RNAxplorer_args_info args_info;

  struct options_s            *options = default_options();

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAxplorer_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* temperature */
  if (args_info.temp_given)
    temperature = args_info.temp_arg;

  /* method */
  if (args_info.method_given) {
    char *m = args_info.method_arg;
    if (!strcmp(m, "MC")) {
      options->strategy = PATHFINDER_SADDLE_MONTE_CARLO;
    } else if (!strcmp(m, "MC-SA")) {
      options->strategy   = PATHFINDER_SADDLE_MONTE_CARLO_SA;
      simulatedAnnealing  = 1;
    } else if (!strcmp(m, "GW")) {
      options->strategy = PATHFINDER_SADDLE_GRADIENT_WALK;
    } else if (!strcmp(m, "DB-MFE")) {
      options->strategy = PATHFINDER_TWO_D_REPRESENTATIVES;
    } else if (!strcmp(m, "BLUBB")) {
      options->strategy = TWO_D_LOWER_BOUND;
    } else if (!strcmp(m, "SM")) {
      options->strategy = ATTRACTION_SAMPLING;
    } else if (!strcmp(m, "RS")) {
      /* Repellant Sampling */
      options->strategy = REPELLENT_SAMPLING;
    }
  }

  /* maximum number of simulations / iterations */
  if (args_info.extended_opt_given)
    extended_options = strdup(args_info.extended_opt_arg);

  /* maximum number of simulations / iterations */
  if (args_info.iterations_given)
    options->iterations = args_info.iterations_arg;

  /* maxkeep for Flamm et al. direct path heuristics */
  if (args_info.maxKeep_given)
    options->max_keep = args_info.maxKeep_arg;

  /* Amount of best solutions to store per iteration */
  if (args_info.maxStore_given)
    options->max_storage = args_info.maxStore_arg;

  if (args_info.circ_given)
    options->md.circ = 1;

  if (args_info.cooling_rate_given)
    treduction = args_info.cooling_rate_arg;

  if (args_info.tstart_given)
    tstart = args_info.tstart_arg + K0;

  if (args_info.tstop_given)
    tstop = args_info.tstop_arg + K0;

  if (args_info.penalizeBackWalks_given)
    backWalkPenalty = 1;

  if (args_info.basinStructure_given)
    options->strategy = MOVE_GRADIENT_WALK;

  if (args_info.maxD_given)
    options->max_d1 = options->max_d2 = args_info.maxD_arg;

  if (args_info.maxD1_given)
    options->max_d1 = args_info.maxD1_arg;

  if (args_info.maxD2_given)
    options->max_d2 = args_info.maxD2_arg;

  if (args_info.betaScale_given)
    options->md.betaScale = args_info.betaScale_arg;

  if (args_info.p0_given) {
    int     i, j = 1, lmintmp;
    double  poptmp = 0.;

    length_indicesAndPercentages  = args_info.p0_given;
    indicesAndPercentages         = (double *)calloc(2 * args_info.p0_given + 1, sizeof(double));
    *indicesAndPercentages        = 1;
    for (i = 0; i < args_info.p0_given; i++) {
      if (sscanf(args_info.p0_arg[i], "%d=%lg", &lmintmp, &poptmp) == 0)
        exit(EXIT_FAILURE);

      if (lmintmp < 1) {
        fprintf(stderr, "States in --p0 must be >=1\n");
        exit(EXIT_FAILURE);
      } else {
        *(indicesAndPercentages + j)      = (double)lmintmp;
        *(indicesAndPercentages + j + 1)  = poptmp;
        *indicesAndPercentages            += 2;
        j                                 += 2;
      }
    }
  }

  /* free allocated memory of command line data structure */
  RNAxplorer_cmdline_parser_free(&args_info);

  return options;
}


char **
extract_structures(unsigned int n,
                   int          maybe_multiline,
                   char         **rec_rest)
{
  char    **structures = NULL;
  size_t  num_structures = 0;
  size_t  size_structures = 10;
  size_t  l, l_prev;
  int     i, read_on;

  if ((rec_rest) && (rec_rest[0])) {
    structures = (char **)vrna_alloc(sizeof(char *) * size_structures);

    read_on = 0;
    l_prev  = 0;

    for (i = 0; rec_rest[i]; i++) {
      switch (rec_rest[i][0]) {
        case '\0':  /* fall-through */
        case '#':   /* fall-through */
        case '%':   /* fall-through */
        case ';':   /* fall-through */
        case '/':   /* fall-through */
        case '*':   /* fall-through */
        case ' ':   /* fall-through */
        case '\t':
          break;

        case '(':
        case ')':
        case '.':
        case '+':
          l = strlen(rec_rest[i]);
          if (l + l_prev < n) {
            if (maybe_multiline) {
              read_on = 1;
            } else {
              vrna_message_error(
                "sequence and structure (at line %d) have unequal lengths (%u vs. %u)",
                i,
                n,
                l + l_prev);
            }
          } else if (l + l_prev > n) {
            vrna_message_error(
              "sequence and structure (at line %d) have unequal lengths (%u vs. %u)",
              i,
              n,
              l + l_prev);
          }

          if (l_prev == 0)
            structures[num_structures] = (char *)vrna_alloc(sizeof(char) * (n + 1));

          memcpy(structures[num_structures] + l_prev, &(rec_rest[i][0]), sizeof(char) * l);
          structures[num_structures][l_prev + l] = '\0';

          if (l + l_prev == n) {
            num_structures++;
            l_prev  = 0;
            read_on = 0;
            if (num_structures == size_structures - 1) {
              size_structures *= 1.4;
              structures      = (char **)vrna_realloc(structures, sizeof(char *) * size_structures);
            }
          } else if (read_on) {
            l_prev = l;
          }

          break;
      }

      free(rec_rest[i]);
    }

    free(rec_rest);

    if (num_structures > 0) {
      structures =
        (char **)vrna_realloc(structures, sizeof(char *) * (num_structures + 1));
      structures[num_structures] = NULL;
    } else {
      free(structures);
      structures = NULL;
    }
  }

  return structures;
}


/*
 *  #############################################################
 *  # Below are the wrappers for the different operation modes  #
 *  #############################################################
 */

/*
 *  *******************
 *  * 1. Move Modes   *
 *  *******************
 */
int
moves_gradient_descent(const char       *rec_id,
                       const char       *orig_sequence,
                       char             **structures,
                       struct options_s *opt)
{
  char *rec_sequence = strdup(orig_sequence);

  vrna_seq_toRNA(rec_sequence);
  vrna_seq_toupper(rec_sequence);

  initRNAWalk(rec_sequence, &(opt->md));

  for (int i = 0; structures[i]; i++) {
    char *basinStructure = structureWalk(rec_sequence, structures[i], GRADIENT_WALK);
    fprintf(stdout, "%4d\t%s\n", i, basinStructure);
    free(basinStructure);
  }

  freeRNAWalkArrays();

  return 1; /* success */
}


/*
 *  ****************************
 *  * 2. Path / Barrier Modes  *
 *  ****************************
 */
int
paths_findpath(const char       *rec_id,
               const char       *orig_sequence,
               char             **structures,
               struct options_s *opt)
{
  char        *rec_sequence = strdup(orig_sequence);
  vrna_path_t *foldingPath, *Saddle, *r;

  vrna_seq_toRNA(rec_sequence);
  vrna_seq_toupper(rec_sequence);

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute direct folding path");
    free(rec_sequence);
    return 0; /* failure */
  }

  foldingPath = get_path(rec_sequence,
                         structures[0],
                         structures[1],
                         opt->max_keep);

  Saddle = getSaddlePoint(foldingPath);

  fprintf(stdout, "# direct Path:\n# barrier: %6.2f\n\n", Saddle->en);
  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  free(rec_sequence);

  return 1; /* success */
}


int
paths_pathfinder_gd(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   GRADIENT_WALK,
                   opt->max_storage);

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_mc(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   MC_METROPOLIS,
                   opt->max_storage);

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_mcsa(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   MC_METROPOLIS,
                   opt->max_storage);

  free(sequence);

  return 1; /* success */
}


int
paths_pathfinder_db(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char        *sequence;
  vrna_path_t *foldingPath, *Saddle, *r;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute folding path");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  foldingPath = levelSaddlePoint2(sequence,
                                  structures[0],
                                  structures[1],
                                  0,
                                  opt->iterations,
                                  opt->max_keep,
                                  opt->max_storage,
                                  opt->max_d1,
                                  opt->max_d2);

  fprintf(stdout,
          "\n# done\n\n# Path with detours:\n# barrier: %6.2f\n\n",
          getSaddlePoint(foldingPath)->en);

  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  free(sequence);
  fflush(stdout);

  return 1; /* success */
}


int
paths_pathfinder_dbba(const char        *rec_id,
                      const char        *orig_sequence,
                      char              **structures,
                      struct options_s  *opt)
{
  char *sequence;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to compute energy barrier approximation");
    return 0; /* failure */
  }

  sequence = strdup(orig_sequence);
  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  barrier_estimate_2D(sequence,
                      &(opt->md),
                      structures[0],
                      structures[1],
                      opt->max_d1,
                      opt->max_d2);

  free(sequence);

  return 1; /* success */
}


/*
 *  *********************
 *  * 3. Sampling Modes *
 *  *********************
 */
int
sampling_repulsion(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt)
{
  char                  *sequence = strdup(orig_sequence);

  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  vrna_fold_compound_t  *fc = vrna_fold_compound(sequence,
                                                 &(opt->md),
                                                 VRNA_OPTION_DEFAULT);

  repellant_sampling(fc);

  vrna_fold_compound_free(fc);
  free(sequence);

  return 1; /* success */
}


int
sampling_attraction(const char        *rec_id,
                    const char        *orig_sequence,
                    char              **structures,
                    struct options_s  *opt)
{
  char                  *sequence;
  int                   num_ref = 0;
  vrna_fold_compound_t  *fc;
  gridLandscapeT        *grid;

  if ((!structures) || (!structures[0]) || (!structures[1])) {
    vrna_message_warning("Too few structures to perform directed sampling");
    return 0; /* failure */
  }

  /* count number of reference structures */
  for (num_ref = 0; structures[num_ref]; num_ref++);

  sequence = strdup(orig_sequence);

  vrna_seq_toRNA(sequence);
  vrna_seq_toupper(sequence);

  fc = vrna_fold_compound(sequence,
                          &(opt->md),
                          VRNA_OPTION_DEFAULT);

#if 0
  grid = estimate_landscape(fc,
                            structures,
                            num_ref,
                            opt->iterations,
                            extended_options);
#endif
  grid = estimate_landscapeMD(fc,
                              (const char **)structures,
                              num_ref,
                              opt->iterations,
                              extended_options,
                              indicesAndPercentages,
                              length_indicesAndPercentages);

  printLandscape(grid, fc);

  free_gridLandscape(grid);

  vrna_fold_compound_free(fc);
  free(sequence);

  return 1; /* success */
}


int
sampling_temperature(const char       *rec_id,
                     const char       *orig_sequence,
                     char             **structures,
                     struct options_s *opt)
{
  vrna_message_warning("Not implemented yet!");
  return 1; /* success */
}
