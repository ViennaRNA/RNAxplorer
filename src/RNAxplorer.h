#ifndef __RNAXPLORER__
#define __RNAXPLORER__

#define FIND_BEST_FOLDINGPATH           0
#define FIND_BASIN_STRUCTURE            1
#define FIND_DISTANCE_BASED_MFE_PATH    2
#define FIND_2D_BARRIER_ESTIMATE        3
#define FIND_2D_LANDSCAPE_ESTIMATE      4
#define DISTORTED_SAMPLING              5

void    RNAxplorer(void);

void    GetBasinStructure(void);

//char    *pair_table_to_dotbracket(short *pt);

void    print_structure(short* pt, int E);

#endif
