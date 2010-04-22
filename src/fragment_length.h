#include <stdio.h>

#define LOWER_FL_CUTOFF_QUANTILE 0.001
#define UPPER_FL_CUTOFF_QUANTILE 0.99
#define MIN_EMPIRICAL_COVERAGE 50
#define TEMP_FL_ARRAY_GF 100

struct fragment_length_dist_t {
    int min_fl;
    int max_fl;
    float* density;
};

void
init_fl_dist( struct fragment_length_dist_t** fl_dist, int min_fl, int max_fl );

void
free_fl_dist( struct fragment_length_dist_t* fl_dist );

void
init_fl_dist_from_file( struct fragment_length_dist_t** fl_dist, FILE* fp );

float
get_fl_prb( struct fragment_length_dist_t* fl_dist, int fl );

void
print_fl_dist( struct fragment_length_dist_t* fl_dist );
