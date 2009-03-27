#ifndef __PATHFINDER__
#define __PATHFINDER__

#define MIN2(A, B)        ((A) < (B) ? (A) : (B))
#define MAX2(A, B)        ((A) > (B) ? (A) : (B))

#define FIND_BEST_FOLDINGPATH   0
#define FIND_BASIN_STRUCTURE    1
#define FIND_DISTANCE_BASED_MFE_PATH    2

void    PathFinder(void);
void    levelSaddlePoint(char *s1, char *s2);
path_t  *levelSaddlePoint2(char *s1, char *s2/*, int *num_entry*/, int iteration);
path_t  *getSaddlePoint(path_t *foldingPath/*, int steps*/);
void    GetBasinStructure(void);

short   *copy_pair_table(short *template);
char    *pair_table_to_dotbracket(short *pt);

void    print_structure(short* pt, int E);

#endif
