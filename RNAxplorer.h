#ifndef __PATHFINDER__
#define __PATHFINDER__

#define MIN2(A, B)        ((A) < (B) ? (A) : (B))
#define MAX2(A, B)        ((A) > (B) ? (A) : (B))

#define FIND_BEST_FOLDINGPATH           0
#define FIND_BASIN_STRUCTURE            1
#define FIND_DISTANCE_BASED_MFE_PATH    2
#define KLKIN                           3
#define TRANSITION_RATES                4

void    RNAxplorer(void);
void    levelSaddlePoint(char *s1, char *s2);
path_t  *levelSaddlePoint2(char *s1, char *s2/*, int *num_entry*/, int iteration);
path_t  *getSaddlePoint(path_t *foldingPath/*, int steps*/);
void    GetBasinStructure(void);

char    *pair_table_to_dotbracket(short *pt);

void    print_structure(short* pt, int E);

/**
*** Calculate the transition rates between k,l-macrostates
*** This function creates a transition rate matrix and a fake
*** barrier output to use with treekin
***
*** rates are calculated as follows
*** /f[
*** k_{\alpha \rightarrow \beta} = \frac{1}{N}\sum_{i \in \alpha}^{N} \sum_{j \in \beta \cap \mathcal{N}(i)}\\
*** /f]
*** as we only sample the microstates we have to take care about the detailed balanced condition:
*** /f[
*** k_{\alpha \rightarrow \beta} \cdot \pi_{\alpha} = k_{\beta \rightarrow \alpha} \cdot \pi_{\beta}
*** /f]
*** we do this by updating k_{\beta \rightarrow \alpha} in each step we generate the neighboring structure
*** /f$j \in \beta /f$ with /f$ j \in \mathcal{N}(i) /f$
*** then for each microrate /f$k_{i \rightarrow j} /f$ the reverse microrate /f$k_{j \rightarrow i} \cdot \frac{\pi_{j | \beta}{i | \alpha} /f$
*** is added to /f$k_{j \rightarrow i} /f$
*** at the end of all rate calculations, the rate /f$q_{i,i} = -\sum_{i \neq j} q_{i,j} /f$ is calculated
***
***
**/
void transition_rates(const char *s, const char *s1, const char *s2);

#endif
