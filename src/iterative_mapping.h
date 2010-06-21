#ifndef ITERATIVE_MAPPING_H
#define ITERATIVE_MAPPING_H

#include "config.h"
struct genome_data;
#include "mapped_read.h"

/*
 * Global arguments that are available to all of the 
 * iterative mapping functions. I should probably 
 * pass these through the calling chain, but that
 * seems like a pain so I just declare them global.
 *
 */
static const void* iterative_mapping_args; 

#define WINDOW_SIZE 20

struct update_mapped_read_rv_t {
    /* Store the probability of observing the 
       returned sequence, conditional on the trace */
    double log_lhd;
    /* Store the max change in conditional mapping 
       probability for the returned read */
    double max_change;
};

/* Forward declaration for the trace type */
struct trace_t;
struct mapped_read_t;

void
naive_update_trace_expectation_from_location( 
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc );


void
update_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc )
);


struct update_mapped_read_rv_t
update_mapped_reads_from_trace( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
);

void
bootstrap_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc)
    );

double
calc_log_lhd ( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
);

int
update_mapping(
    struct mapped_reads_db* rdb,
    struct trace_t* starting_trace,
    int max_num_iterations,
    float max_prb_change_for_convergence,
    
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc),

    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
    );

void
build_random_starting_trace( 
    struct trace_t* traces, 
    struct mapped_reads_db* rdb,
    
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc),

    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
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

    int num_tracks,
    char** track_names,

    int num_samples,
    int max_num_iterations,
    float max_prb_change_for_convergence,
    
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc),
    
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
                          
    );

/*****************************************************************************
 * 
 * ChIP Seq specific functions 
 *
 *****************************************************************************/

void 
update_chipseq_trace_expectation_from_location(
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc );

struct update_mapped_read_rv_t 
update_chipseq_mapped_read_prbs( const struct trace_t* const traces, 
                                 const struct mapped_read_t* const r  );

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
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc );

struct update_mapped_read_rv_t 
update_CAGE_mapped_read_prbs( 
    const struct trace_t* const traces, 
    const struct mapped_read_t* const r  );

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
generic_update_mapping( struct rawread_db_t* rawread_db,
                        struct mapped_reads_db* rdb, 
                        struct genome_data* genome,
                        enum assay_type_t,
                        int num_samples,
                        float max_prb_change_for_convergence);


#endif
