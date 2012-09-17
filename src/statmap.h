/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef STATMAP_H
#define STATMAP_H

#include "config.h"
#include "config_parsing.h"
#include "rawread.h"
#include "genome.h"

extern int num_threads;
extern int min_num_hq_bps;
extern int max_reference_insert_len;

extern int num_trace_tracks;
extern char** trace_track_names;

/* Fwd declarations */
struct mapped_reads_db;

void
map_marginal( struct args_t* args, 
              struct genome_data* genome,
              struct rawread_db_t* rdb,
              struct mapped_reads_db** mpd_rds_db,
              enum bool is_nc );

void
build_fl_dist( struct args_t* args, struct mapped_reads_db* mpd_rds_db );

void
iterative_mapping( struct args_t* args, 
                   struct genome_data* genome,
                   struct mapped_reads_db* mpd_rds_db );

int 
main( int argc, char** argv );

#endif /* STATMAP_H */
