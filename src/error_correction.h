#ifndef ERROR_CORRECTION_H
#define ERROR_CORRECTION_H

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rembedded.h>
#include <Rmath.h>

#include <pthread.h>
#include <limits.h>
#include <assert.h>

#include "read.h"
#include "genome.h"
#include "fragment_length.h"

#define NUM_BASE_SAMPLES_FOR_MIN_PENALTY_COMP 5000
#define NUM_READ_SAMPLES_FOR_MIN_PENALTY_COMP 10
#define ROUND_ERROR 1e-6

struct indexable_subtemplates; // fwd declaration (?)
struct read_subtemplate;
struct penalty_t;
struct read;

struct freqs_array {
    int max_qual_score;
    int max_position;
    double** freqs;
};

struct error_data_t;

void
predict_freqs( 
    struct error_data_t* data, 
    int record_index, 
    struct freqs_array* predicted_freqs 
);

/*
 * Error model functions
 *
 */

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

/* Mapping metaparameters are set by the user and fixed. They are used by the
 * error model code as guides to set more fine grained parameters. */
struct mapping_metaparams {
    enum error_model_type_t error_model_type;
    /* Static struct (for now) */
    float error_model_params[10];
};

struct index_search_params
{
    float min_match_penalty;
    float max_penalty_spread;
};

/* 3 penalties: read penalty, rst penalty, (index probe penalty) */
struct sampled_penalties_t
{
    float read_penalty;
    float read_subtemplate_penalty;
};

struct mapping_params {
    struct mapping_metaparams* metaparams;
    
    int num_penalty_arrays;
    struct penalty_array_t** fwd_penalty_arrays;
    struct penalty_array_t** rev_penalty_arrays;
    
    int total_read_length;
    
    float recheck_min_match_penalty;
    float recheck_max_penalty_spread;
};

float
expected_value_of_rst_subsequence(
        struct penalty_array_t* rst_pens,
        int subseq_start,
        int subseq_length
    );

float
expected_value_of_rst(
        struct read_subtemplate* rst,
        struct penalty_array_t* rst_pens
    );

int
compute_sampled_penalties_for_reads(
        struct rawread_db_t* rdb,
        struct error_model_t* error_model,
        int num_reads,
        float quantile
    );

struct mapping_params*
init_mapping_params_for_read(
        struct read* r,        
        struct mapping_metaparams* metaparams,
        struct error_model_t* error_model
    );

struct index_search_params*
init_index_search_params(
        struct indexable_subtemplates* ists,
        struct mapping_params* mapping_params
    );

void
free_mapping_params( struct mapping_params* p );

/*
 *  Functions for determining read mappability.
 *
 */


int
calc_effective_sequence_length( struct penalty_array_t* penalties );

/* determine whether reads are mappable */
enum bool
filter_read(
        struct read* r,
        struct mapping_params* mapping_params,
        struct genome_data* genome
    );

enum bool
filter_indexable_subtemplates(
    struct indexable_subtemplates* ists,
    struct mapping_params* mapping_params,
    struct genome_data* genome
    );


/*
 *  Functions for saving error information
 *
 */

struct error_data_record_t {
    /* number of reads processed */
    int num_unique_reads;

    /* maximum read length processed */
    int max_read_length;
    int max_qual_score;

    int read_subtemplate_index;
    int strand;
    
    /* the readkey range that this record covers */
    int min_readkey;
    int max_readkey;
    
    /* 
       2D array of counts per base. The first dimension
       stores the qual scores, the second positions. So
       base_type_cnts[qual_score][pos] returns the cnts 
       for qual_score, pos 
    */
    int** base_type_cnts;
    int** base_type_mismatch_cnts;
};

struct error_data_t {
    int num_records;
    int max_read_length;
    int max_qual_score;
    
    struct error_data_record_t** records;
    
    /* mutex used by the global error_data_t struct for thread safety */
    pthread_mutex_t* mutex;
};


void
init_error_data( struct error_data_t** data );

void 
free_error_data( struct error_data_t* data );

void
add_new_error_data_record( 
    struct error_data_t* data, int min_readkey, int max_readkey );

/*
 * Merge record into the error_data_record number i in data->records.
 * If record_index == -1, use the last record.
 *
 */
void
merge_in_error_data( struct error_data_t* data_to_update, 
                     struct error_data_t* data );

void
find_length_and_qual_score_limits( struct error_data_t* data,
                                   int* min_qual_score, int* max_qual_score,
                                   int* max_read_length );

void log_error_data( FILE* ofp, struct error_data_t* data );


/*******************************************************************************
 *
 *
 * Error record data
 *
 *
 ******************************************************************************/


void
init_error_data_record( struct error_data_record_t** data, 
                        int read_subtemplate_index, enum STRAND strand,
                        int max_read_len, int max_qual_score );

void
free_error_data_record( struct error_data_record_t* data );

void
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read_seq,
    char* error_str,
    int read_length,
    int subtemplate_index,
    enum STRAND strand,
    int location_offset
);

void
sum_error_data_records(
    struct error_data_record_t* dest,
    struct error_data_record_t* src
);

void
fprintf_error_data_record( 
    FILE* stream, struct error_data_record_t* data,
    int min_qual_score, int max_qual_score, int max_read_length );

/*******************************************************************************
 *
 *
 * Error data collection
 *
 *
 ******************************************************************************/

void*
update_error_data_from_index_search_results(
    struct read_subtemplate* rst,
    mapped_locations** search_results, 
    struct genome_data* genome, 
    struct error_data_t* error_data);


#endif
