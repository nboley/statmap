#ifndef ITERATIVE_MAPPING_H
#define ITERATIVE_MAPPING_H

#include "config.h"
#include "genome.h"
#include "mapped_location.h"

#define HQ_THRESHOLD 0.95

#define LOCATIONS_GROWTH_FACTOR 1000

#define SKIP_THIS_BP

#define WINDOW_SIZE 20


/******************************************************************************
 * Chr Trace Code
 *
 * Chr traces are just arrays that store cnt expectations at each basepair. 
 * Below is the code for manipulating them.
 *
 *****************************************************************************/

#define TRACE_TYPE float

typedef struct {
    int num_traces;
    unsigned int* trace_lengths;
    TRACE_TYPE** fwd_traces;
    TRACE_TYPE** bkwd_traces;
} stranded_traces_t;

typedef struct {
    /* Usually, this is just the number of chrs */
    int num_traces;
    unsigned int* trace_lengths;
    TRACE_TYPE** traces;
} traces_t;

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size );

void
init_stranded_traces( const genome_data* const genome,
                      stranded_traces_t** traces );

void
close_stranded_traces( stranded_traces_t* traces );


void
init_traces( const genome_data* const genome,
             traces_t** traces );

void
close_traces( traces_t* traces );

/*
double
sum_traces( TRACE_TYPE** chr_traces, 
            unsigned int num_chrs, 
            unsigned int* chr_lens );
*/
/*
void
renormalize_traces( TRACE_TYPE** chr_traces, 
                    int num_chrs, 
                    unsigned int* chr_lens,
                    unsigned long num_reads );
*/

void
write_mapped_reads_to_wiggle( mapped_reads_db* rdb, 
                              genome_data* genome,
                              FILE* wfp );


void
write_wiggle_from_traces( traces_t* traces,
                          char** trace_names,
                          
                          const char* output_fname, 
                          const char* track_name,
                          
                          const double filter_threshold );

/*
 * END Chr Trace Code
 *
 *****************************************************************************/

/*****************************************************************************
 * 
 * Mapped Short Reads Fns
 *
 * Methods for dealing with mapped short reads in the context of iterative 
 * updates.
 *
 *****************************************************************************/


void
mmap_mapped_short_reads_file( char* fname, 
                              char** mapped_short_reads, 
                              size_t* mapped_short_reads_size );

void
build_reads_start_array( char* reads_data, size_t reads_data_size,
                         char*** reads_starts_array, unsigned long* num_reads );

/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/

/* Integrated trace/rawread info */

void
update_traces_from_mapped_chipseq_reads( 
    mapped_reads_db* reads_db,
    traces_t* traces
);


double
update_mapped_chipseq_reads_from_trace( 
    mapped_reads_db* reads_db,
    traces_t* traces
);

int 
update_chipseq_mapping( 
    mapped_reads_db* rdb, 
    genome_data* genome,
    int max_num_iterations 
);

int
update_mapping( mapped_reads_db* rdb, 
                genome_data* genome,
                int max_num_iterations,
                enum assay_type_t       );

#endif
