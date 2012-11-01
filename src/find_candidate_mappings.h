/* Copyright (c) 2009-2012, Nathan Boley */

#include <pthread.h>
#include <stdio.h>

#include "quality.h"
#include "genome.h"
#include "rawread.h"
#include "candidate_mapping.h"
#include "mapped_location.h"
#include "error_correction.h"

#define MAX_NUM_INDEX_PROBES 3

#define MAX_NUM_UNTEMPLATED_GS 1
#define UNTEMPLATED_G_MARGINAL_LOG_PRB -1.30103

#define CANT_BUILD_READ_SUBTEMPLATES 1

struct single_map_thread_data {
    int thread_id;
    struct genome_data* genome;

    struct rawread_db_t* rdb;
    
    unsigned int* mapped_cnt;
    pthread_spinlock_t* mapped_cnt_lock;
    readkey_t max_readkey;
    
    struct mapped_reads_db* mpd_rds_db;
    
    struct mapping_metaparams* metaparams;    
    float reads_min_match_penalty;
    struct error_model_t* error_model;
    struct error_data_t* error_data;
    enum bool only_collect_error_data;
};

void*
find_candidate_mappings( void* params );

void
bootstrap_estimated_error_model( 
        struct genome_data* genome,
        struct rawread_db_t* rdb,
        struct mapped_reads_db* mpd_rds_db, // TODO set to NULL for bootstrap?
        struct mapping_metaparams* mapping_metaparams,
        struct error_model_t* error_model
    ); 

void
find_all_candidate_mappings(
        struct genome_data* genome,
        struct rawread_db_t* rdb,
        struct mapped_reads_db* mpd_rds_db,

        struct mapping_metaparams* mapping_metaparams,
        struct error_model_t* error_model
    );

float
subseq_penalty(
        struct read_subtemplate* rst,
        int subseq_offset,
        int subseq_length,
        struct penalty_array_t* penalties
    );