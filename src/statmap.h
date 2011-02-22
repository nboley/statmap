/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef STATMAP_H
#define STATMAP_H

#include "config.h"

struct genome_data;

extern int num_threads;
extern int min_num_hq_bps;

extern int num_trace_tracks;
extern char** trace_track_names;


/* fwd declarations */
struct rawread_db_t;
struct args_t;
enum input_file_type_t;

void usage();

/*
 * Try and determine the file type. 
 *
 * In particular, we try and determine what the sequence mutation
 * string types are.
 *
 * The method is to scan the first 10000 reads and record the 
 * max and min untranslated scores. 
 *
 */

enum input_file_type_t
guess_input_file_type( struct args_t* args );

/*
 * Guess the optimal indexed sequence length.
 *
 * To do this, we open up the first read file and then scan for a read. The 
 * seq length of the first read is what we set the index length to.
 *
 */
int
guess_optimal_indexed_seq_len( struct args_t* args);

struct args_t
parse_arguments( int argc, char** argv );

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
