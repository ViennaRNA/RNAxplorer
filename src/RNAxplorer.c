#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <ViennaRNA/model.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/io/file_formats.h>
#include <ViennaRNA/datastructures/hash_tables.h>
#include <ViennaRNA/landscape/move.h>
#include <ViennaRNA/landscape/walk.h>
#include <ViennaRNA/landscape/neighbor.h>

#include "RNAwalk.h"
#include "meshpoint.h"
#include "barrier_lower_bound.h"
#include "distorted_sampling.h"
#include "distorted_samplingMD.h"
#include "repellant_sampling.h"
#include "PathFinder.h"
#include "gradient_walker.h"

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
  TEMPERATURE_SCALING_SAMPLING,
  REPELLENT_SAMPLING_HEURISTIC,
  RETRIEVE_LOCAL_MINIMA             /* perform gradient walks for and arbitrary set of structures */
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

  rnax_path_finder_opt_t  pathfinder;

  /* retrieve local minima options */
  double temperature_celsius;
  int shift_moves;
  char *parameter_file;

  /* repulsive sampling options*/
  char *sequence;
  char *struc1;
  char *struc2;
  int granularity;
  int num_samples;
  float exploration_factor;
  float min_exploration_percent;
  int cluster; //flag
  char *lmin_file;
  char *TwoD_file;
  int non_red; //flag
  char *non_red_file;
  int explore_two_neighborhood; //flag
  int post_filter_two; //flag
  int ediff_penalty; //flag
  float mu;
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

int
sampling_repellent_heuristic(const char       *rec_id,
                   const char       *orig_sequence,
                   char             **structures,
                   struct options_s *opt);

int
retrieve_local_minima(const char       *rec_id,
                      const char       *orig_sequence,
                      char             **structures,
                      struct options_s *opt);


/**
*** \file RNAxplorer.c
**/


#define NUM_STRATEGIES    12

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
  },
  {
    REPELLENT_SAMPLING_HEURISTIC,
    &sampling_repellent_heuristic,
    "Repellent Sampling Heuristic"
  },
  {
      RETRIEVE_LOCAL_MINIMA,
      &retrieve_local_minima,
      "Retrieve local minima for an arbitrary set of secondary structures"
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
  int read_fasta_records = 0;
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

    if (rec_type & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT)){
        if(read_fasta_records == 0){
            vrna_message_input_seq("Input at least one valid fasta record! (a header line '> header' and a sequence line (upper or lower case) followed by structures)");
        }
        break;
    }

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
    read_fasta_records++;

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

  /* free options */
  free(options->TwoD_file);
  free(options->non_red_file);
  free(options->sequence);
  free(options->struc1);
  free(options->struc2);

  return EXIT_SUCCESS;
}


struct options_s *
default_options(void)
{
  struct options_s *options = (struct options_s *)vrna_alloc(sizeof(struct options_s));

  /* default strategy */
  options->strategy = PATHFINDER_SADDLE_GRADIENT_WALK;

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
  options->temperature_celsius = 37;
  if (args_info.temp_given){
    temperature = args_info.temp_arg;
    options->temperature_celsius = args_info.temp_arg;
  }
  options->shift_moves = 0;
  if (args_info.shift_moves_flag){
      options->shift_moves = 1;
  }
  options->parameter_file = NULL;
  if (args_info.parameter_file_given){
      options->parameter_file = args_info.parameter_file_arg;
  }

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
    } else if (!strcmp(m, "RSH")) {
      /* Repellant Sampling Heuristic*/
      options->strategy = REPELLENT_SAMPLING_HEURISTIC;
    } else if (!strcmp(m, "RL")) {
      options->strategy = RETRIEVE_LOCAL_MINIMA;
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

  /* repulsive sampling args */
  if(args_info.sequence_given){
      options->sequence = vrna_alloc(strlen(args_info.sequence_arg)+1);
      strcpy(options->sequence, args_info.sequence_arg);
  }
  if(args_info.struc1_given){
      options->struc1 = vrna_alloc(strlen(args_info.struc1_arg)+1);
      strcpy(options->struc1, args_info.struc1_arg);
  }
  if(args_info.struc2_given){
      options->struc2 = vrna_alloc(strlen(args_info.struc2_arg)+1);
      strcpy(options->struc2, args_info.struc2_arg);
  }

  if(args_info.granularity_given){
      options->granularity = args_info.granularity_arg;
  }
  if(args_info.num_samples_given){
      options->num_samples = args_info.num_samples_arg;
  }
  if(args_info.exploration_factor_given){
      options->exploration_factor = args_info.exploration_factor_arg;
  }
  if(args_info.min_exploration_percent_given){
      options->min_exploration_percent = args_info.min_exploration_percent_arg;
  }
  if(args_info.cluster_given){
      options->cluster = args_info.cluster_flag;
  }
  if(args_info.lmin_file_given){
      options->lmin_file = vrna_alloc(strlen(args_info.lmin_file_arg)+1);
      strcpy(options->lmin_file, args_info.lmin_file_arg);
  }
  if(args_info.TwoD_file_given){
      options->TwoD_file = vrna_alloc(strlen(args_info.TwoD_file_arg)+1);
      strcpy(options->TwoD_file, args_info.TwoD_file_arg);
  }
  if(args_info.nonred_given){
      options->non_red = args_info.nonred_flag;
  }
  if(args_info.nonred_file_given){
      options->non_red_file = vrna_alloc(strlen(args_info.nonred_file_arg)+1);
      strcpy(options->non_red_file, args_info.nonred_file_arg);
  }
  if(args_info.explore_two_neighborhood_given){
      options->explore_two_neighborhood = args_info.explore_two_neighborhood_flag;
  }
  if(args_info.post_filter_two_given){
      options->post_filter_two = args_info.post_filter_two_flag;
  }
  if(args_info.ediff_penalty_given){
      options->ediff_penalty = args_info.ediff_penalty_flag;
  }
  if(args_info.mu_given){
      options->mu = args_info.mu_arg;
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

  vrna_fold_compound_t *fc = vrna_fold_compound(rec_sequence,
                                                &(opt->md),
                                                VRNA_OPTION_EVAL_ONLY);

  for (int i = 0; structures[i]; i++) {
    short *pt = vrna_ptable(structures[i]);

    (void)vrna_path(fc, pt, 0, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);

    char *loc_min = vrna_db_from_ptable(pt);
    fprintf(stdout, "%4d\t%s\n", i, loc_min);

    free(loc_min);
    free(pt);
  }

  vrna_fold_compound_free(fc);
  free(rec_sequence);

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

  vrna_fold_compound_t *fc = vrna_fold_compound(rec_sequence,
                                                &(opt->md),
                                                VRNA_OPTION_EVAL_ONLY);
  
  foldingPath = vrna_path_findpath(fc,
                                   structures[0],
                                   structures[1],
                                   opt->max_keep);

  Saddle = getSaddlePoint(foldingPath);

  fprintf(stdout, "# direct Path:\n# barrier: %6.2f\n\n", Saddle->en - foldingPath->en);
  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  vrna_fold_compound_free(fc);
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

#if 1
  vrna_path_t *foldingPath, *Saddle, *r;
  rnax_path_finder_opt_t *options = rnax_path_finder_options();
  options->md          = opt->md;
  options->iterations  = opt->iterations;
  options->max_keep    = opt->max_keep;

  foldingPath = rnax_path_finder(sequence,
                                 structures[0],
                                 structures[1],
                                 options);

  fprintf(stdout,
          "# barrier: %6.2f\n\n",
          getSaddlePoint(foldingPath)->en);

  for (r = foldingPath; r->s; r++)
    fprintf(stdout, "%s %6.2f\n", r->s, r->en);

  free_path(foldingPath);
  free(options);
  fflush(stdout);
#else
  levelSaddlePoint(sequence,
                   structures[0],
                   structures[1],
                   opt->iterations,
                   opt->max_keep,
                   GRADIENT_WALK,
                   opt->max_storage);
#endif

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


struct sc_data{
    vrna_hash_table_t base_pairs;
    vrna_hash_table_t weights;
};

typedef struct key_value_ {
  vrna_move_t *key;
  int value;
} key_value;

static unsigned
hash_function_base_pairs(void           *hash_entry,
                   unsigned long  hashtable_size)
{
  key_value *kv = (key_value *)hash_entry;
  unsigned long  hash_value = ((unsigned long)kv->key->pos_5) << 32;
  hash_value = hash_value & kv->key->pos_3;
  return hash_value;
}


static int
hash_comparison_base_pairs(void *x,
                     void *y)
{
  key_value *hem_x  = ((key_value *)x);
  key_value *hem_y  = ((key_value *)y);

  if ((x == NULL) ^ (y == NULL))
    return 1;

  return !(hem_x->key->pos_5 == hem_y->key->pos_5 && hem_x->key->pos_3 == hem_y->key->pos_3);
}


static int
free_base_pairs(void *x)
{
    //do nothing (free whole array)
  return 0;
}

static vrna_hash_table_t
create_hashtable(int hashbits)
{
  vrna_callback_ht_free_entry       *my_free          = free_base_pairs;
  vrna_callback_ht_compare_entries  *my_comparison    = hash_comparison_base_pairs;
  vrna_callback_ht_hash_function    *my_hash_function = hash_function_base_pairs;
  vrna_hash_table_t                 ht                = vrna_ht_init(hashbits,
                                                                     my_comparison,
                                                                     my_hash_function,
                                                                     my_free);

  return ht;
}

typedef struct hashtable_list_ {
  unsigned long     length;
  unsigned long     allocated_size;
  float      *list_weights;
  int        *list_counts;
  key_value         **list_key_value_pairs;
  vrna_hash_table_t ht_pairs; // lookup table;
} hashtable_list;

static hashtable_list
create_hashtable_list(int hashbits)
{
  hashtable_list ht_list;

  ht_list.allocated_size          = 10;
  ht_list.length                  = 0;
  ht_list.list_weights = vrna_alloc(sizeof(float) * ht_list.allocated_size);
  ht_list.list_counts    = vrna_alloc(sizeof(int) * ht_list.allocated_size);
  ht_list.list_key_value_pairs    = vrna_alloc(sizeof(key_value *) * ht_list.allocated_size);
  ht_list.ht_pairs         = create_hashtable(hashbits);
  return ht_list;
}


static
void
free_hashtable_list(hashtable_list *ht_list)
{
  vrna_ht_free(ht_list->ht_pairs);
  free(ht_list->list_weights);
  free(ht_list->list_counts);

  //int i = 0;
  //for (; i < ht_list->length; i++)
  //  free(ht_list->list_key_value_pairs[i]);

  free(ht_list->list_key_value_pairs);
}

static
void
hashtable_list_add_weight_and_count(hashtable_list *htl, vrna_move_t key, float weight)
{
  if (htl->ht_pairs != NULL) {
    key_value to_check;
    to_check.key = &key;
    //to_check->value = 0; //not checked anyways --> not set

    //to_check.key = energy;
    //to_check.value = count;
    key_value *lookup_result = NULL;
    lookup_result = vrna_ht_get(htl->ht_pairs, (void *)&to_check);
    if (lookup_result == NULL) {
      //value is not in list.
      if (htl->length >= htl->allocated_size) {
        htl->allocated_size           += 10;
        htl->list_weights  =
          vrna_realloc(htl->list_weights, sizeof(int) * htl->allocated_size);
        htl->list_counts  =
          vrna_realloc(htl->list_counts, sizeof(float) * htl->allocated_size);
        htl->list_key_value_pairs = vrna_realloc(htl->list_key_value_pairs,
                                                 sizeof(key_value *) * htl->allocated_size);
      }

      float weight;
      int           list_index = htl->length;
      htl->list_weights[list_index]  = weight;
      htl->list_counts[list_index]  = 1;
      key_value *to_insert = vrna_alloc(sizeof(key_value));
      to_insert->key = vrna_alloc(sizeof(vrna_move_t));
      to_insert->key->pos_5 = key.pos_5;
      to_insert->key->pos_3 = key.pos_3;
      to_insert->value = list_index;
      htl->list_key_value_pairs[list_index] = to_insert;
      htl->length++;
      int           res         = vrna_ht_insert(htl->ht_pairs, (void *)to_insert);
      if (res != 0)
        fprintf(stderr, "dos.c: hash table insert failed!");
    } else {
      // the energy-index pair is already in the list.
      int list_index = lookup_result->value;
      htl->list_counts[list_index] += 1;
      htl->list_weights[list_index] += weight;
    }
  }
}

void store_basepair_sc(vrna_fold_compound_t *fc, hashtable_list *data, char *structure, float weight, int distance_based /*= False */){
    if(distance_based == 1){
        return;
    }
    else{
        short *pt = vrna_ptable(structure);
        // count number of pairs in structure to repell
        int cnt = 0;
        for(int i = 1; i < pt[0]; i++){
            if(pt[i] > i){
                cnt = cnt + 1;
            }
        }
        if(cnt > 0){
            weight = weight / (float)cnt;
        }
        // add repulsion
        int i;
        for(i = 1; i < pt[0]; i++){
            if(pt[i] > i){
                vrna_move_t key;
                key.pos_5 = i;
                key.pos_3 = pt[i];
                hashtable_list_add_weight_and_count(data, key, weight);
                /*
                if key not in data['base_pairs']:
                    data['base_pairs'][key] = 1
                    data['weights'][key] = weight
                else:
                    data['base_pairs'][key] = data['base_pairs'][key] + 1
                    data['weights'][key] = data['weights'][key] + weight
                 */
            }
        }
        // remove previous soft constraints
        //fc->sc_remove();
        vrna_sc_init(fc);

        // add latest penalties for base pairs
        int j;
        int list_index;
        for(int k = 0; k < data->length; k++){ // in data['weights'].keys():
            i = data->list_key_value_pairs[k]->key->pos_5;
            j = data->list_key_value_pairs[k]->key->pos_3;
            list_index = data->list_key_value_pairs[k]->value;
            int pair_weight = data->list_weights[list_index];
            //fc->sc_add_bp(i, j, pair_weight);
            vrna_sc_add_bp(fc, i, j, pair_weight, VRNA_OPTION_DEFAULT);
        }
    }
}

short * detect_local_minimum(vrna_fold_compound_t *fc, short *structure_pt){
    // perform gradient walk from sample to determine direct local minimum
    //pt = RNA.IntVector(RNA.ptable(structure))
    //fc->path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    vrna_path_gradient(fc, structure_pt, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);
    return structure_pt; //RNA.db_from_ptable(list(pt))
}

/**
 * Take a local minimum and detect nearby local minimum via extended gradient walks.
 * If a lower energy structure within radius 2 is detected (lower than the local minimum
 * of a normal gradient walk), then another gradient walk is applied for the lower structure.
 * @brief detect_local_minimum_two
 * @param fc
 * @param structure_pt
 * @return
 */
short * detect_local_minimum_two(vrna_fold_compound_t *fc, short *structure_pt){

    short *deepest_neighbor = vrna_ptable_copy(structure_pt);
    int deepest_neighbor_ddG  = 1;
    //pt                    = RNA.IntVector(RNA.ptable(structure))
    vrna_move_t *neigh                 = vrna_neighbors(fc, structure_pt, VRNA_MOVESET_DELETION | VRNA_MOVESET_INSERTION);

    vrna_move_t *nb;
    for(nb = neigh; nb->pos_5 != 0 && nb->pos_3 != 0; nb++){
        int dG_nb   = vrna_eval_move_pt(fc, structure_pt, nb->pos_5, nb->pos_3);
        short *pt_nb = vrna_ptable_copy(structure_pt);
        vrna_move_apply(pt_nb, nb);
        vrna_move_t *neigh_2                = vrna_neighbors(fc, pt_nb, VRNA_MOVESET_DELETION | VRNA_MOVESET_INSERTION);

         vrna_move_t *nb_2;
        for(nb_2 = neigh_2; nb_2->pos_5 != 0 && nb_2->pos_3 != 0; nb_2++){
            int dG_nb2   = vrna_eval_move_pt(fc, pt_nb, nb->pos_5, nb->pos_3);

            int ddG     = dG_nb + dG_nb2;
            if(ddG < 0 && ddG < deepest_neighbor_ddG){
                deepest_neighbor_ddG  = ddG;
                memcpy(deepest_neighbor, pt_nb, sizeof(short)*(pt_nb[0]+1));
                vrna_move_apply(deepest_neighbor, nb_2);

            }
        }
        free(neigh_2);
    }
    free(neigh);
    if(memcmp(deepest_neighbor, structure_pt, sizeof(short)*(structure_pt[0]+1)) != 0){
        free(deepest_neighbor);
        deepest_neighbor = detect_local_minimum(fc, deepest_neighbor);
        return detect_local_minimum_two(fc, deepest_neighbor);
    }
    else{
        free(deepest_neighbor);
        return structure_pt;
    }
}


void reduce_lm_two_neighborhood(vrna_fold_compound_t *fc, short **lm, int verbose /*= False*/){
//TODO: translate this
    /*
    cnt       = 1;
    cnt_max   = len(lm)
    lm_remove = list()
    lm_novel  = dict()

    if verbose:
        sys.stderr.write("Applying 2-Neighborhood Filter...")
        sys.stderr.flush()

    for s in lm:
        if verbose:
            sys.stderr.write("\rApplying 2-Neighborhood Filter...%6d / %6d" % (cnt, cnt_max))
            sys.stderr.flush()

        ss = detect_local_minimum_two(fc, s)
        if s != ss:
            # store structure for removal
            lm_remove.append(s)
            if ss not in lm:
                # store structure for novel insertion
                if ss not in lm_novel:
                    lm_novel[ss] = { 'count' : lm[s]['count'], 'energy' : fc.eval_structure(ss) }
                else:
                    lm_novel[ss]['count'] = lm_novel[ss]['count'] + lm[s]['count']
            else:
                # update current minima list
                lm[ss]['count'] = lm[ss]['count'] + lm[s]['count']

        cnt = cnt + 1

    # remove obsolete local minima
    for a in lm_remove:
        del lm[a]

    # add newly detected local minima
    lm.update(lm_novel)

    if verbose:
        sys.stderr.write("\rApplying 2-Neighborhood Filter...done             \n")
  */
}

short ** generate_samples(vrna_fold_compound_t *fc, int number, int non_redundant /*=False */){
    //samples = list()
    short **samples;

    if(non_redundant){
        //TODO: samples = fc.pbacktrack_nr(number)
    }
    else{
        int i;
        for(i =0; i < number; i++){
            //TODO samples.append(fc.pbacktrack())
        }
    }
    return samples;
}

int
sampling_repellent_heuristic(const char       *rec_id,
                             const char       *orig_sequence,
                             char             **structures,
                             struct options_s *opt){
    if(opt->sequence && strlen(opt->sequence) == 0){
        fprintf(stderr, "Error: the input sequence is not given!");
        exit(1);
    }
    //if(strlen(opt->sequence) != strlen(opt->struc1) || strlen(opt->sequence) != strlen(opt->struc2)){
    //    fprintf(stderr, "Error: the input sequence has to have the same length as structure 1 and structure 2!");
    //    exit(1);
    //}
    if(rec_id)
        printf("%s\n",rec_id);


    char *sequence = NULL;
    char *structure1 = NULL;
    char *structure2 = NULL;

    if(structures){
        int i;
        for(i = 0; structures[i]; i++){
            // count only
        }
        if(i>=1 && strlen(orig_sequence) == strlen(structures[0]) && strlen(orig_sequence) == strlen(structures[1])){
            sequence = vrna_alloc(strlen(orig_sequence));
            strcpy(sequence, orig_sequence);
            structure1 = structures[0];
            structure2 = structures[1];
        }
    }
    if(sequence == NULL && opt->sequence != NULL && opt->struc1 != NULL && opt->struc1 != NULL &&
            strlen(opt->sequence) == strlen(opt->struc1) && strlen(orig_sequence) == strlen(opt->struc2)){
        sequence = vrna_alloc(strlen(opt->sequence));
        strcpy(sequence, opt->sequence);
        structure1 = opt->struc1;
        structure2 = opt->struc2;
    }
    if(sequence == NULL){
        fprintf(stderr, "Error: provide at least a fasta file with structures as standard input or use the command line parameters for sequence, structure 1 and structure 2!\n");
        exit(EXIT_FAILURE);
    }
    printf("%s\n%s\n%s\n",opt->sequence,opt->struc1,opt->struc2);


    //struct sc_data *my_sc_sdata;
    //my_sc_sdata.base_pairs = create_hashtable(27);
    hashtable_list sc_data = create_hashtable_list(27);

/*

sc_data = {
  'base_pairs': {},
  'weights': {},
}


def store_basepair_sc(fc, data, structure, weight, distance_based = False):
    if distance_based:
        return
    else:
        pt = RNA.ptable(structure)
        # count number of pairs in structure to repell
        cnt = 0
        for i in range(1, len(pt)):
            if pt[i] > i:
                cnt = cnt + 1

        if cnt > 0:
            weight = weight / cnt

        # add repulsion
        for i in range(1, len(pt)):
            if pt[i] > i:
                key = (i, pt[i])
                if key not in data['base_pairs']:
                    data['base_pairs'][key] = 1
                    data['weights'][key] = weight
                else:
                    data['base_pairs'][key] = data['base_pairs'][key] + 1
                    data['weights'][key] = data['weights'][key] + weight

        # remove previous soft constraints
        fc.sc_remove()

        # add latest penalties for base pairs
        for k in data['weights'].keys():
            i = k[0]
            j = k[1]
            fc.sc_add_bp(i, j, data['weights'][k])


def move_apply(structure_in, move):
    #s_length = structure_in[0]
    output_structure = list(structure_in)
    if(move.pos_5 > 0 and move.pos_3 > 0):
        output_structure[move.pos_5] = move.pos_3
        output_structure[move.pos_3] = move.pos_5
    elif(move.pos_5 < 0 and move.pos_3 < 0):
        output_structure[-move.pos_5] = 0
        output_structure[-move.pos_3] = 0
    else:
        print("Error: shfit moves are not supported")
    res = output_structure
    return res


def detect_local_minimum(fc, structure):
    # perform gradient walk from sample to determine direct local minimum
    pt = RNA.IntVector(RNA.ptable(structure))
    fc.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
    return RNA.db_from_ptable(list(pt))


def detect_local_minimum_two(fc, structure):
    """
    Take a local minimum and detect nearby local minimum via extended gradient walks.
    If a lower energy structure within radius 2 is detected (lower than the local minimum
    of a normal gradient walk), then another gradient walk is applied for the lower structure.
    """
    deepest_neighbor      = None
    deepest_neighbor_ddG  = 1
    pt                    = RNA.IntVector(RNA.ptable(structure))
    neigh                 = fc.neighbors(pt, RNA.MOVESET_DELETION | RNA.MOVESET_INSERTION)

    for nb in neigh:
        dG_nb   = fc.eval_move_pt(pt, nb.pos_5, nb.pos_3)
        pt_nb   = RNA.IntVector(move_apply(pt, nb))
        neigh_2 = fc.neighbors(pt_nb, RNA.MOVESET_DELETION | RNA.MOVESET_INSERTION)

        for nb_2 in neigh_2:
            dG_nb2  = fc.eval_move_pt(pt_nb, nb_2.pos_5, nb_2.pos_3)
            ddG     = dG_nb + dG_nb2
            if ddG < 0 and ddG < deepest_neighbor_ddG:
                deepest_neighbor_ddG  = ddG
                deepest_neighbor      = move_apply(pt_nb, nb_2)

    if deepest_neighbor:
        deepest_neighbor = detect_local_minimum(fc, RNA.db_from_ptable(deepest_neighbor))
        return detect_local_minimum_two(fc, deepest_neighbor)
    else:
        return structure


def reduce_lm_two_neighborhood(fc, lm, verbose = False):
    cnt       = 1;
    cnt_max   = len(lm)
    lm_remove = list()
    lm_novel  = dict()

    if verbose:
        sys.stderr.write("Applying 2-Neighborhood Filter...")
        sys.stderr.flush()

    for s in lm:
        if verbose:
            sys.stderr.write("\rApplying 2-Neighborhood Filter...%6d / %6d" % (cnt, cnt_max))
            sys.stderr.flush()

        ss = detect_local_minimum_two(fc, s)
        if s != ss:
            # store structure for removal
            lm_remove.append(s)
            if ss not in lm:
                # store structure for novel insertion
                if ss not in lm_novel:
                    lm_novel[ss] = { 'count' : lm[s]['count'], 'energy' : fc.eval_structure(ss) }
                else:
                    lm_novel[ss]['count'] = lm_novel[ss]['count'] + lm[s]['count']
            else:
                # update current minima list
                lm[ss]['count'] = lm[ss]['count'] + lm[s]['count']

        cnt = cnt + 1

    # remove obsolete local minima
    for a in lm_remove:
        del lm[a]

    # add newly detected local minima
    lm.update(lm_novel)

    if verbose:
        sys.stderr.write("\rApplying 2-Neighborhood Filter...done             \n")


def generate_samples(fc, number, non_redundant=False):
    samples = list()

    if non_redundant:
        samples = fc.pbacktrack_nr(number)
    else:
        for i in range(0, number):
            samples.append(fc.pbacktrack())

    return samples


"""
Do main stuff
"""
# init random number generator in RNAlib
RNA.init_rand()

# prepare RNAlib fold_compound
md = RNA.md()
md.uniq_ML = 1
md.compute_bpp = 0

kT = RNA.exp_param(md).kT

fc = RNA.fold_compound(sequence, md)
fc_base = RNA.fold_compound(sequence, md)

# compute MFE structure
(ss, mfe) = fc.mfe()
fc.pf()


minima = dict()
sample_list = []

pending_lm = dict()
num_sc = 1

num_iter      = int(math.ceil(float(num_samples) / float(granularity)))
samples_left  = num_samples

for it in range(0, num_iter):
    # determine number of samples for this round
    # usually, this is 'granularity'
    if samples_left < granularity:
        current_num_samples = samples_left
    else:
        current_num_samples = granularity

    samples_left = samples_left - current_num_samples

    # generate samples through stocastic backtracing
    sample_set = generate_samples(fc, current_num_samples, nonredundant)

    sys.stderr.write("\rsamples so far: %6d / %6d" % (num_samples - samples_left, num_samples))
    sys.stderr.flush()

    # store samples of this round to global list of samples
    sample_list = sample_list + sample_set

    current_lm = dict()

    # go through list of sampled structures and determine corresponding local minima
    for s in sample_set:
        ss = detect_local_minimum(fc_base, s)

        if ss not in current_lm:
            current_lm[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
        else:
            current_lm[ss]['count'] = current_lm[ss]['count'] + 1

    # explore the 2-neighborhood of current local minima to reduce total number of local minima
    if explore_two_neighborhood:
        reduce_lm_two_neighborhood(fc_base, current_lm, verbose)

    # transfer local minima obtained in this iteration to list of pending local minima
    for ss in current_lm:
        if ss not in pending_lm:
            pending_lm[ss] = current_lm[ss]
        else:
            pending_lm[ss]['count'] = pending_lm[ss]['count'] + current_lm[ss]['count']

    del current_lm

    if it < num_iter - 1:
        # find out which local minima we've seen the most in this sampling round
        struct_cnt_max = max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
        struct_cnt_min = min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['count']))
        struct_en_max = max(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))
        struct_en_min = min(pending_lm.iterkeys(), key=(lambda a: pending_lm[a]['energy']))

        cnt_once  = 0
        cnt_other = 0

        for key, value in sorted(pending_lm.iteritems(), key=lambda (k, v): v['energy']):
            if value['count'] == 1:
                cnt_once = cnt_once + 1
            else:
                cnt_other = cnt_other + 1

            # check whether we've seen other local minim way more often than those we've seen just once
            if cnt_other > 0 and cnt_once > 0:
                if cnt_once <= (mu * cnt_other):
                    if ediff_penalty:
                        repell_en = kt_fact * (value['energy'] - pending_lm[struct_en_min]['energy'])
                    else:
                        repell_en = kt_fact * kT / 1000.

                    store_basepair_sc(fc, sc_data, struct_cnt_max, repell_en)

                    fc.pf()

                    for cmk in pending_lm.keys():
                        if cmk not in minima:
                            minima[cmk] = pending_lm[cmk]
                        else:
                            minima[cmk]['count'] = minima[cmk]['count'] + pending_lm[cmk]['count']

                    pending_lm = dict()

                    break
    else:
        for cmk in pending_lm.keys():
            if cmk not in minima:
                minima[cmk] = pending_lm[cmk]
            else:
                minima[cmk]['count'] = minima[cmk]['count'] + pending_lm[cmk]['count']

eprint(" ... done")

if post_filter_two:
    reduce_lm_two_neighborhood(fc_base, minima, verbose)


RNAlocmin_output(sequence, minima, lmin_file)


sample_file=""
rind = lmin_file.rfind(".")
if rind >= 0 :
    sample_file = lmin_file[:rind] + ".samples"
else:
    sample_file = lmin_file[:rind] + ".samples"
f = open(sample_file, 'w')
f.write("     %s\n" % sequence)
for s in sample_list:
    f.write(s+"\n")
f.close()


if fake_2D_file:
    (ss, mfe) = fc_base.mfe()
    RNA2Dfold_output(sequence, ss, mfe, structure1, structure2, minima, TwoD_file)


# read a list of sample structures and produce list of local minima for it
if nonredundant_sample_file:
    lmin_nonred_file = "local_minima_nonred.txt"
    nonredundant_samples = []
    with open(nonredundant_sample_file) as f:
        nonredundant_samples = f.readlines()

    nonredundant_samples = [x.strip() for x in nonredundant_samples]

    nonredundant_samples.pop(0)

    nonredundant_minima = dict()

    for s in nonredundant_samples:
        pt = RNA.IntVector(RNA.ptable(s))
        fc_base.path(pt, 0, RNA.PATH_DEFAULT | RNA.PATH_NO_TRANSITION_OUTPUT)
        ss = RNA.db_from_ptable(list(pt))
        if ss not in nonredundant_minima:
             nonredundant_minima[ss] = { 'count' : 1, 'energy' : fc_base.eval_structure(ss) }
             new_minima = new_minima + 1
        else:
             nonredundant_minima[ss]['count'] = nonredundant_minima[ss]['count'] + 1

    f = open(lmin_nonred_file, 'w')
    f.write("     %s\n" % sequence)
    for i,s in enumerate(sorted(nonredundant_minima.keys(), key=lambda x: nonredundant_minima[x]['energy'])):
        f.write("%4d %s %6.2f %6d\n" % (i, s, nonredundant_minima[s]['energy'], nonredundant_minima[s]['count']))
    f.close()

    if fake_2D_file:
        distances = [ [ None for j in range(0, 200) ] for i in range(0, 200) ];

        for s in nonredundant_minima.keys():
            d1 = RNA.bp_distance(structure1, s)
            d2 = RNA.bp_distance(structure2, s)
            if not distances[d1][d2] or nonredundant_minima[s]['energy'] < distances[d1][d2]:
                distances[d1][d2] = nonredundant_minima[s]['energy']

        f = open("sv11_fake_nonred.2D.out", 'w')
        f.write("%s\n%s (%6.2f)\n%s\n%s\n\n\n" % (sequence, structure1, mfe, structure1, structure2))
        for i in range(0, 200):
            for j in range(0, 200):
                if distances[i][j] != None:
                    f.write("%d\t%d\ta\tb\tc\t%6.2f\n" % (i, j, distances[i][j]))
        f.close()

*/


    /* free */
    free(sequence);

    return EXIT_SUCCESS;
}

int
retrieve_local_minima(const char       *rec_id,
                      const char       *orig_sequence,
                      char             **structures,
                      struct options_s *opt){
  double temperature_celsius = opt->temperature_celsius;
  int shift_moves = opt->shift_moves;
  char *parameter_file = opt->parameter_file;
  gradient_walker(temperature_celsius, shift_moves, parameter_file, orig_sequence, structures);

}


