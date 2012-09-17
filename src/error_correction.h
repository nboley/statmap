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

/* This code is preparation for future merge of error_correction branch */
/* Mapping metaparameters are set by the user and fixed. They are used by the
 * error model code as guides to set more fine grained parameters. */
struct mapping_metaparams {
    enum error_model_type_t error_model_type;
    /* Static struct (for now) */
    float error_model_params[10];
};

/*
 * Store separate mapping parameters For
 *  a) index search
 *  b) the final recheck.
 *
 * For now, use the command line parameters for the recheck penalties and
 * generate the index search penalties from the error model.
 */
struct index_search_params
{
    float min_match_penalty;
    float max_penalty_spread;
};

struct mapping_params {
    struct mapping_metaparams* metaparams;

    float recheck_min_match_penalty;
    float recheck_max_penalty_spread;
};

void
init_mapping_params_for_read(
        struct mapping_params** p,
        struct mapping_metaparams* metaparams,
        float recheck_min_match_penalty,
        float recheck_max_penalty_spread
    );

void
free_mapping_params( struct mapping_params* p );

void
find_mapping_params( struct mapping_params* params,
                     struct read_subtemplate* rst,
                     struct genome_data* genome,
                     struct error_model_t* error_model_t
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
merge_in_error_data_record( struct error_data_t* data, int record_index,
                            struct error_data_record_t* record );

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
                        int max_read_len, int max_qual_score );

void
free_error_data_record( struct error_data_record_t* data );

void
update_error_data_record(
    struct error_data_record_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int read_length
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

#endif
