/* Copyright (c) 2009-2012, Nathan Boley */

#include <pthread.h>
#include <stdio.h>

#include "quality.h"
#include "genome.h"
#include "rawread.h"
#include "candidate_mapping.h"
#include "mapped_location.h"

struct single_map_thread_data {
    int thread_id;
    struct genome_data* genome;

    struct rawread_db_t* rdb;
    
    unsigned int* mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex;
    readkey_t max_readkey;
    
    candidate_mappings_db* mappings_db;
    pthread_mutex_t* mappings_db_mutex;
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
    candidate_mappings_db* mappings_db,
    
    struct error_model_t* error_model
);

void
find_all_candidate_mappings( 
    struct genome_data* genome,
                             
    struct rawread_db_t* rdb,
    candidate_mappings_db* mappings_db,
                             
    struct error_model_t* error_model,
                             
    float min_match_penalty,
    float max_penalty_spread
);

