#ifndef ITERATIVE_MAPPING_H
#define ITERATIVE_MAPPING_H

#include "config.h"
#include "genome.h"
#include "mapped_location.h"

#define HQ_THRESHOLD 0.95

#define LOCATIONS_GROWTH_FACTOR 1000

#define SKIP_THIS_BP

#define WINDOW_SIZE 20


/* Forward declaration for the trace type */
#include "trace.h"

void
naive_update_trace_expectation_from_location( 
    const traces_t* const traces, 
    const mapped_read_location* const loc );


void
update_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    traces_t* traces,
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc)
);


double
update_mapped_reads_from_trace( 
    struct mapped_reads_db* reads_db,
    traces_t* traces,
    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
);

int
update_mapping(
    struct mapped_reads_db* rdb,
    traces_t* starting_trace,
    int max_num_iterations,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),

    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
    );

void
build_random_starting_trace( 
    traces_t* traces, 
    struct mapped_reads_db* rdb,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),

    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
    );

/*****************************************************************************
 * 
 * High Level functions 
 *
 *****************************************************************************/

int
sample_random_traces( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    int trace_dim,
    int num_samples,
    int max_num_iterations,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),
    
    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
                          
    );

/*****************************************************************************
 * 
 * ChIP Seq specific functions 
 *
 *****************************************************************************/

void 
update_chipseq_trace_expectation_from_location(
    const traces_t* const traces, 
    const mapped_read_location* const loc );

double 
update_chipseq_mapped_read_prbs( const traces_t* const traces, 
                                 const mapped_read* const r  );

int
update_chipseq_mapping( struct mapped_reads_db* rdb, 
                        struct genome_data* genome,
                        int max_num_iterations );

/* END ChIP Seq specific functions  ******************************************/

/*****************************************************************************
 * 
 * CAGE specific functions 
 *
 *****************************************************************************/

void update_CAGE_trace_expectation_from_location(
    const traces_t* const traces, 
    const mapped_read_location* const loc );

double update_CAGE_mapped_read_prbs( 
    const traces_t* const traces, 
    const mapped_read* const r  );

int
update_cage_mapping( struct mapped_reads_db* rdb, 
                     struct genome_data* genome,
                     int max_num_iterations );

/* END CAGE specific functions  **********************************************/

void
write_mapped_reads_to_wiggle( struct mapped_reads_db* rdb, 
                              struct genome_data* genome,
                              FILE* wfp );

/*
 * END Chr Trace Code
 *
 *****************************************************************************/

/* fwd declaration for the assay type */
enum assay_type_t;

int
generic_update_mapping( struct mapped_reads_db* rdb, 
                        struct genome_data* genome,
                        enum assay_type_t       );

#endif
