/* Copyright (c) 2009-2012, Nathan Boley */

#include <pthread.h>
#include <stdio.h>

#include "quality.h"
#include "genome.h"
#include "rawread.h"
#include "candidate_mapping.h"
#include "mapped_location.h"

void
search_index(
        struct index_t* index, 
        struct indexable_subtemplate* ist,

        float min_match_penalty,
        float max_penalty_spread,
        mapped_locations** results,

        enum bool only_find_unique_sequence
    );

struct single_map_thread_data {
    int thread_id;
    struct genome_data* genome;

    struct rawread_db_t* rdb;
    
    unsigned int* mapped_cnt;
    pthread_spinlock_t* mapped_cnt_lock;
    readkey_t max_readkey;
    
    struct mapped_reads_db* mpd_rds_db;

    float min_match_penalty;
    float max_penalty_spread;
    
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
    
    struct error_model_t* error_model
); 

void
find_all_candidate_mappings(
        struct genome_data* genome,

        struct rawread_db_t* rdb,
        struct mapped_reads_db* mpd_rds_db,

        struct error_model_t* error_model,

        float min_match_penalty,
        float max_penalty_spread
    );
