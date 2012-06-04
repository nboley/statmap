#include <stdio.h>

#define LOWER_FL_CUTOFF_QUANTILE 0.001
#define UPPER_FL_CUTOFF_QUANTILE 0.99
#define MIN_EMPIRICAL_COVERAGE 50
#define TEMP_FL_ARRAY_GF 10
#define MIN_NUM_READS 10

struct fragment_length_dist_t {
    int min_fl;
    int max_fl;
    float* density;
    float* chipseq_bs_density;
    float* rev_chipseq_bs_density;

    /* to user the vector operations, we need to 
       make sure that the chipseq bs densities
       are correctly aligned at 16 bytes. However, 
       we still need to free the memory so, below,
       we store the start of the allcd addresses */
    float* chipseq_bs_density_start;
    float* rev_chipseq_bs_density_start;
};

void
init_fl_dist( struct fragment_length_dist_t** fl_dist, int min_fl, int max_fl );

void
free_fl_dist( struct fragment_length_dist_t** fl_dist );

void
init_fl_dist_from_file( struct fragment_length_dist_t** fl_dist, FILE* fp );

/* FWD declaration */
struct mapped_reads_db;

void
estimate_fl_dist_from_mapped_reads( struct mapped_reads_db* rdb );

void
build_chipseq_bs_density( struct fragment_length_dist_t* fl_dist );

float
get_fl_prb( struct fragment_length_dist_t* fl_dist, int fl );

void
fprint_fl_dist( FILE* fp, struct fragment_length_dist_t* fl_dist );
