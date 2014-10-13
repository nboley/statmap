#ifndef ITERATIVE_MAPPING_H
#define ITERATIVE_MAPPING_H

#include "config.h"
#include "genome.h"
#include "mapped_read.h"

/*
 * Global arguments that are available to all of the 
 * iterative mapping functions. I should probably 
 * pass these through the calling chain, but that
 * seems like a pain so I just declare them global.
 *
 */
// BUG
// static const void* iterative_mapping_args; 

#define WINDOW_SIZE 100

#define CHIPSEQ_LHD_RATIO_STOP_VALUE 1.05
#define CAGE_LHD_RATIO_STOP_VALUE -1

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
//mapped_read_t;

int
update_mapping(
    struct mapped_reads_db* rdb,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* starting_trace,
    int max_num_iterations,
    float max_prb_change_for_convergence,
    float lhd_ratio_stop_value,
    
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob),

    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           mapped_read_t* r  )
    );

struct cond_prbs_db_t*
build_posterior_db( struct genome_data* genome, 
                    struct mapped_reads_db* mpd_rdb,
                    enum assay_type_t assay_type );

#endif
