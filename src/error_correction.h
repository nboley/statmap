#include <pthread.h>

#define ERROR_WEIGHT    0.5
#define ERROR_STATS_LOG "error_stats.log"

#define max_num_qual_scores 256

enum error_model_type_t {
    MISMATCH    = 1,
    FASTQ_MODEL = 2,
    ESTIMATED   = 3
};

struct error_data_t;

struct error_model_t {
    enum error_model_type_t error_model_type;
    void* data;
};

void
init_error_model( 
    struct error_model_t** error_model,
    enum error_model_type_t error_model_type
 );

void
update_error_model_from_error_data( 
    struct error_model_t* error_model,
    struct error_data_t* data
);

void
free_error_model( struct error_model_t* error_model );

/*
 *  Functions for saving error information
 *
 */

struct error_data_t {
    /* number of reads processed */
    int num_unique_reads;

    /* maximum read length processed */
    int max_read_length;
    
    /* array of counters of mismatches at each index in a read */
    int* position_cnts;
    int* position_mismatch_cnts;
    
    /* quality score match and mismatch counts */
    int qual_score_cnts[max_num_qual_scores];
    int qual_score_mismatch_cnts[max_num_qual_scores];

    /* mutex used by the global error_data_t struct for thread safety */
    pthread_mutex_t* mutex;
};

void
init_error_data( struct error_data_t** data );

void 
free_error_data( struct error_data_t* data );

void
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int length
);

void
add_error_data(
    struct error_data_t* dest,
    struct error_data_t* src
);

void
fprintf_error_data( FILE* stream, struct error_data_t* data );

void 
log_error_data( struct error_data_t* ed );
