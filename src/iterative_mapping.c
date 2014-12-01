/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <sys/time.h> /* gettimeofday() */
#include <signal.h>
#include <sys/wait.h>


#include <fcntl.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>
// #include <omp.h>
#include <errno.h>
#include <float.h>

#include "statmap.h"
#include "config_parsing.h"
#include "iterative_mapping.h"
#include "mapped_read.h"
#include "trace.h"
#include "wiggle.h"
#include "mapped_location.h"
#include "genome.h"
#include "fragment_length.h"
#include "sam.h"
#include "util.h"
#include "log.h"

struct fragment_length_dist_t* global_fl_dist;
struct genome_data* global_genome;
struct trace_t* global_starting_trace;

#define BIN_SIZE 1
#define WINDOW_SIZE 100

/* 
 * This is what update_traces_from_mapped_reads calls - it is a function
 * designed to be called from a thread.
 * 
 */ 
struct update_traces_param {
    /* 
     * we lock the read index, so that threads grab the next available read
     * I like this better than a block of reads approach, because it ensures that 
     * any io is sequential. Of course, this is at the cost of some lock 
     * contention, but that should be minor, espcially if ( the mmapped ) rdb
     * can't fully fit into memory.
     */
    struct mapped_reads_db* reads_db;
    struct cond_prbs_db_t* cond_prbs_db;
    struct trace_t* traces;
    void (* update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob );
};

void*
update_traces_from_mapped_reads_worker( void* params )
{
    struct trace_t* traces = ( (struct update_traces_param*) params)->traces;
    
    struct mapped_reads_db* rdb = 
        ( (struct update_traces_param*) params)->reads_db;

    struct cond_prbs_db_t* cond_prbs_db =
        ( (struct update_traces_param*) params)->cond_prbs_db;

    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob )
        = ( (struct update_traces_param*) 
            params)->update_trace_expectation_from_location;
    
    mapped_read_t* r;

    mapped_read_index* rd_index;
    heap_allocate_mapped_read_index(rd_index);
    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) )     
    {
        init_mapped_read_index( r, rd_index );

        /* Update the trace from this mapping */        
        MPD_RD_ID_T j;
        for( j = 0; j < rd_index->num_mappings; j++ ) {
            float cond_prob = get_cond_prb(
                cond_prbs_db, rd_index->read_id, j );
            
            update_trace_expectation_from_location( 
                traces, rd_index->mappings[j], cond_prob );
        }
    }
    
    free_mapped_read_index( rd_index );
    free( rd_index );
    pthread_exit( NULL );
}

void
update_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob )
)
{        
    zero_traces( traces );
    
    /* Move the database ( this should be a cursor ) back to the first read */
    rewind_mapped_reads_db( reads_db );
    
    /* If the number of threads is one, then just do everything in serial */
    if( num_threads == 1 )
    {
        mapped_read_t* r;
        int read_num = 0;
        
        mapped_read_index* rd_index;
        heap_allocate_mapped_read_index(rd_index);
        while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) )     
        {
            read_num++;
            init_mapped_read_index( r, rd_index );

            /* Update the trace from this mapping */        
            MPD_RD_ID_T j;
            for( j = 0; j < rd_index->num_mappings; j++ ) {
                float cond_prob = get_cond_prb(
                    cond_prbs_db, rd_index->read_id, j );
                
                update_trace_expectation_from_location( 
                    traces, rd_index->mappings[j], cond_prob );
            }
        }
        free_mapped_read_index( rd_index );
        free(rd_index);
    } 
    /* otherwise, if we are expecting more than one thread */
    else {
        /* initialize the thread parameters structure */
        struct update_traces_param params;
        
        params.reads_db = reads_db;
        params.cond_prbs_db = cond_prbs_db;
        params.traces = traces;
        params.update_trace_expectation_from_location 
            = update_trace_expectation_from_location;
        
        /* initialize all of the threads */
        int rc;
        void* status;
        pthread_t thread[num_threads];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        struct update_traces_param* tds;
        tds = malloc(num_threads*sizeof(struct update_traces_param));
        int t;
        for( t = 0; t < num_threads; t++ )
        {  
            memcpy( tds+t,  &params, sizeof(struct update_traces_param) );
            
            rc = pthread_create( &(thread[t]), 
                                 &attr, 
                                 update_traces_from_mapped_reads_worker, 
                                 (void *)(tds + t) 
                ); 

            if (rc) {
                statmap_log( LOG_FATAL, "Return code from pthread_create() is %d", rc );
            }
        }

        /* wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            if (rc) {
                statmap_log( LOG_FATAL, "Return code from pthread_join() is %d", rc );
            }
        }
        pthread_attr_destroy(&attr);
        
        free( tds );
    }
    
    return;
}

/* 
 * This is what update_mapped_reads_from_trace calls - it is a function
 * designed to be called from a thread.
 * 
 */ 
struct update_mapped_reads_param {
    struct update_mapped_read_rv_t rv;
    
    /* 
     * we lock the read index, so that threads grab the next available read
     * I like this better than a block of reads approach, because it ensures 
     * any io is sequential. Of course, this is at the cost of some lock 
     * contention, but that should be minor, espcially if ( the mmapped ) rdb
     * can't fully fit into memory.
     */
    struct mapped_reads_db* reads_db;
    struct cond_prbs_db_t* cond_prbs_db;
    
    struct trace_t* traces;
    struct update_mapped_read_rv_t 
        (* update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                     const struct trace_t* const traces, 
                                     mapped_read_t* r  );
    
};

void
update_mapped_reads_from_trace_worker( void* params )
{
    struct trace_t* traces = ( (struct update_mapped_reads_param*) params)->traces;
    struct mapped_reads_db* reads_db = 
        ( (struct update_mapped_reads_param*) params)->reads_db;

    struct cond_prbs_db_t* cond_prbs_db =
        ( (struct update_mapped_reads_param*) params)->cond_prbs_db;

    struct update_mapped_read_rv_t 
        (* update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                     const struct trace_t* const traces, 
                                     mapped_read_t* r  )
        = ( (struct update_mapped_reads_param*)
            params)->update_mapped_read_prbs;

    struct update_mapped_read_rv_t *rv = 
        &(((struct update_mapped_reads_param*) params)->rv);
    
    mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) ) 
    {
        /* Update the read */
        struct update_mapped_read_rv_t tmp_rv 
            = update_mapped_read_prbs( cond_prbs_db, traces, r );
        
        /* Hand reduce the max */
        if( tmp_rv.max_change > rv->max_change )
            rv->max_change = tmp_rv.max_change;

        
        /* Update the lhd */
        rv->log_lhd += tmp_rv.log_lhd;
    }

    return;
}

void*
update_mapped_reads_from_trace_thread_worker( void* params )
{
    /* call the real worker function, and exit the thread when it returns.
     *
     * Wrapping update_mapped_reads_from_trace_worker allows us to use it
     * unmodified in the single thread case. We must wrap it to make the call
     * to pthread_exit because calling pthread_exit in the only thread in
     * a process is equivalent to calling exit (see `man pthread_exit`) */
    update_mapped_reads_from_trace_worker( params );
    pthread_exit( NULL );
}


struct update_mapped_read_rv_t
update_mapped_reads_from_trace( 
    struct mapped_reads_db* reads_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           mapped_read_t* r  )
    )
{    
    /* reset the cursor */
    rewind_mapped_reads_db( reads_db );

    /* initialize the thread parameters structure */
    struct update_mapped_reads_param params;
        
    params.rv.max_change = 0;
    params.rv.log_lhd = 0;
                
    params.reads_db = reads_db;
    params.cond_prbs_db = cond_prbs_db;
        
    /* traces should be read only, so they are shared */
    params.traces = traces;
    params.update_mapped_read_prbs = update_mapped_read_prbs;
    
    /* If the number of threads is one, then just do everything in serial */
    if( num_threads == 1 )
    {
        update_mapped_reads_from_trace_worker( &params );
    } else {        
        /* initialize all of the threads */
        int rc;
        void* status;
        pthread_t thread[num_threads];

        pthread_attr_t* attrs;
        attrs = malloc(sizeof(pthread_attr_t)*num_threads);
        
        struct update_mapped_reads_param* tds;
        tds = malloc(num_threads*sizeof(struct update_mapped_reads_param));
        int t;
        for( t = 0; t < num_threads; t++ )
        {  
            memcpy( tds+t,  &params, sizeof(struct update_mapped_reads_param) );
            
            pthread_attr_init(attrs+t);
            pthread_attr_setdetachstate(attrs+t, PTHREAD_CREATE_JOINABLE);
            
            rc = pthread_create( thread + t, 
                                 attrs + t, 
                                 update_mapped_reads_from_trace_thread_worker, 
                                 (void *)(tds + t) 
                ); 
            if (rc) {
                statmap_log( LOG_FATAL, 
                             "Return code from pthread_create() is %d", rc );
            }
        }

        /* wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            pthread_attr_destroy( attrs + t );
            
            if (rc) {
                statmap_log( LOG_FATAL, 
                             "Return code from pthread_join() is %d", rc );
            }
        }
        
        
        /* Aggregate over the max in the error difference */
        for( t=0; t < num_threads; t++ )
        {
            params.rv.log_lhd += tds[t].rv.log_lhd;
            params.rv.max_change = MAX( 
                tds[t].rv.max_change, params.rv.max_change );
        }
        
        free( tds );
        free( attrs );
    }
    
    return params.rv;
}

double
calc_log_lhd( 
    struct mapped_reads_db* reads_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           mapped_read_t* r  )
    )
{    
    struct update_mapped_read_rv_t rv;
    
    rv = update_mapped_reads_from_trace( 
        reads_db, cond_prbs_db, traces, update_mapped_read_prbs );
    
    return rv.log_lhd;
}

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
        const float prb ),

    struct update_mapped_read_rv_t 
        (*update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                    const struct trace_t* const traces, 
                                    mapped_read_t* r  )
    )
{
    /* make sure the mapped reads are open for reading */
    assert( rdb->mode == 'r' );

    double prev_log_lhd = 0;
    
    struct update_mapped_read_rv_t rv;
    
    struct timeval last_log_time;
    gettimeofday( &last_log_time, NULL );
    
    int num_iterations = 0;
    for( num_iterations = 0; 
         num_iterations < max_num_iterations; 
         num_iterations++ )
    {
        /* Normalize the trace sum to 1 */
        /* This makes the trace the current marginal read density estimate */
        struct timeval nt_start, nt_stop, umr_start, umr_stop, ut_start, ut_stop;
        
        gettimeofday( &nt_start, NULL );   
        normalize_traces( starting_trace );
        gettimeofday( &nt_stop, NULL );   
        
        gettimeofday( &umr_start, NULL );   
        rv = update_mapped_reads_from_trace(
            rdb, cond_prbs_db,
            starting_trace, 
            update_mapped_read_prbs
        );
        gettimeofday( &umr_stop, NULL );   
        
        gettimeofday( &ut_start, NULL );   
        update_traces_from_mapped_reads( 
            rdb, cond_prbs_db, starting_trace, 
            update_trace_expectation_from_location
        );
        gettimeofday( &ut_stop, NULL );

        if( num_iterations > 0 &&
            ( 
                ( 
                  num_iterations == 1 
                  || num_iterations%100 == 0 
                  || (ut_stop.tv_sec-last_log_time.tv_sec) > 30
                )
                || rv.max_change < max_prb_change_for_convergence
                || ( lhd_ratio_stop_value >= 1.0
                     && pow( 10, rv.log_lhd - prev_log_lhd ) < lhd_ratio_stop_value
                     && pow( 10, rv.log_lhd - prev_log_lhd ) > 0.95
                   )
                )
            )
        {
            gettimeofday( &last_log_time, NULL );  
            
            statmap_log( LOG_NOTICE,
            "Iter %i: \tError: %e \tLog Lhd: %e (ratio %e) \tNorm Trace:  %.5f sec\t Read UT:  %.5f sec\tTrace UT:  %.5f sec", 
            num_iterations, rv.max_change, rv.log_lhd,
            pow( 10, rv.log_lhd - prev_log_lhd ),
            (float)(nt_stop.tv_sec - nt_start.tv_sec) 
                   + ((float)(nt_stop.tv_usec - nt_start.tv_usec))/1000000,
            (float)(umr_stop.tv_sec - umr_start.tv_sec) 
                   + ((float)(umr_stop.tv_usec - umr_start.tv_usec))/1000000,
            (float)(ut_stop.tv_sec - ut_start.tv_sec) 
                   + ((float)(ut_stop.tv_usec - ut_start.tv_usec))/1000000
            );
            
            if( num_iterations > 1 
                && ( rv.max_change < max_prb_change_for_convergence 
                || ( lhd_ratio_stop_value >= 1.0
                     && pow( 10, rv.log_lhd - prev_log_lhd ) < lhd_ratio_stop_value 
                     && pow( 10, rv.log_lhd - prev_log_lhd ) > 0.95
                   )
                )
                ) {
                break;
            }
        }        

        prev_log_lhd = rv.log_lhd;
    }
    
    return 0;
}

/*
 * END Chr Trace Code
 *
 *****************************************************************************/

/*****************************************************************************
 * 
 * ChIP Seq specific functions 
 *
 *****************************************************************************/

void 
update_chipseq_trace_expectation_from_location(
    const struct trace_t* const traces, 
    const mapped_read_location* const loc,
    float cond_prob 
)
{
    const MRL_CHR_TYPE chr_index = get_chr_from_mapped_read_location( loc  );
    const MRL_START_POS_TYPE start = get_start_from_mapped_read_location( loc );
    const MRL_STOP_POS_TYPE stop = get_stop_from_mapped_read_location( loc  );

    assert( mapped_read_location_is_paired( loc ) );
    enum bool first_read_is_rev_comp
        = first_read_in_mapped_read_location_is_rev_comp( loc );
    
    assert( cond_prob >= 0.0 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_chrs );
    assert( traces->chr_lengths[chr_index] >= stop );

    const float scale = (1.0/(stop-start));

    int track_index;
    if( first_read_is_rev_comp )
    {
        track_index = 1;
    } else {
        track_index = 0;
    }

    /* get the set of trace_segments to update */
    struct trace_segments_t* trace_segments
        = &(traces->segments[track_index][chr_index]);
        
    float scale_factor = scale*cond_prob;
    update_trace_segments_from_uniform_kernel(
        trace_segments, scale_factor, start, stop);
    
    return;
}

struct update_mapped_read_rv_t 
update_chipseq_mapped_read_prbs( struct cond_prbs_db_t* cond_prbs_db,
                                 const struct trace_t* const traces, 
                                 mapped_read_t* r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };

    /* FIXME cast away const? */
    mapped_read_index* rd_index;
    stack_allocate_and_init_mapped_read_index( rd_index, r );
    
    /* allocate space to store the temporary values */
    float* new_prbs = malloc( sizeof(float)*(rd_index->num_mappings) );
    
    /* Update the reads from the trace */
    /* set to flt_min to avoid division by 0 */
    double prb_sum = ML_PRB_MIN;
    double error_prb_sum = ML_PRB_MIN;

    /* calculate the mean density */
    double window_density = 0;

    int i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        mapped_read_location* current_loc = rd_index->mappings[i];

        MRL_CHR_TYPE chr_index =
            get_chr_from_mapped_read_location( current_loc );
        MRL_START_POS_TYPE start =
            get_start_from_mapped_read_location( current_loc );
        MRL_STOP_POS_TYPE stop =
            get_stop_from_mapped_read_location( current_loc );

        assert( mapped_read_location_is_paired( current_loc ) );
        
        window_density += accumulate_from_traces(
            traces, chr_index, start, stop);
            
        /* 
         * This is the probability of observing the seqeunce given that it came from 
         * location i, assuming the chip traces have been normalized to 0. 
         */
        window_density *= get_fl_prb( 
            global_fl_dist, 
            get_fl_from_mapped_read_location( current_loc ) );
        
        new_prbs[i] = get_seq_error_from_mapped_read_location(
                current_loc )*window_density;
            
        prb_sum += new_prbs[i];        
        error_prb_sum += get_seq_error_from_mapped_read_location( current_loc );
    }    
    
    /* renormalize the read probabilities */
    if( prb_sum > 0 )
    {
        for( i = 0; i < rd_index->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* absolute value error */
            double normalized_density = new_prbs[i]/prb_sum; 
            
            rv.max_change += fabs( 
                get_cond_prb( cond_prbs_db, rd_index->read_id, i ) -
                normalized_density );
            
            set_cond_prb( cond_prbs_db,
                          rd_index->read_id,
                          i,
                          normalized_density );
        }
        /* BUG - Can we be smarter here? */
        /* if the window density sums to 0, then it is impossible to normalize. 
           so we assume uniformity over the reads. I dont know if this is the right
           thing to do ( in fact, im worried that it may destroy convergence. However,
           it's the only thing to think of except ignoring the read, so I do. Note 
           that this should not happen with any trace that is generated by the same 
           stream. ie, this only happens if we are using a marginal density that 
           did not come from the actual reads. ie, when we map the NC reads in a chipseq
           experiment using the IP marginal density.
        */
    } else {
        for( i = 0; i < rd_index->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* absolute value error */
            mapped_read_location* current_loc = rd_index->mappings[i];

            double normalized_density = 
                get_seq_error_from_mapped_read_location( current_loc ) /
                error_prb_sum; 

            set_cond_prb( cond_prbs_db,
                          rd_index->read_id,
                          i,
                          normalized_density );
        }
    }
    
    rv.log_lhd = log10( prb_sum );
    
    /* Free the tmp array */
    free( new_prbs );

    free_mapped_read_index( rd_index );

    return rv;
}

/*****************************************************************************
 * 
 * CAGE specific functions 
 *
 *****************************************************************************/

void update_CAGE_trace_expectation_from_location(
    const struct trace_t* const traces, 
    const mapped_read_location* const loc,
    float cond_prob
)
{
    MRL_CHR_TYPE chr_index = get_chr_from_mapped_read_location( loc  );
    MRL_START_POS_TYPE start = get_start_from_mapped_read_location( loc  );
    MRL_START_POS_TYPE stop = get_stop_from_mapped_read_location( loc  );
    
    enum bool is_paired 
        = mapped_read_location_is_paired( loc );
    enum bool first_read_is_rev_comp
        = first_read_in_mapped_read_location_is_rev_comp( loc );

    assert( cond_prob >= 0.0 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_chrs );
    assert( traces->chr_lengths[chr_index] >= stop );

    assert( cond_prob >= 0.0 );
    
    /* update the trace */
    /* If the reads are paired */
    if( is_paired ) {
        statmap_log(LOG_FATAL, "Paired CAGE reads are not supported" );
    } 
    /* If the read is *not* paired */
    else {
        int track_index;
        int promoter_pos = -1;
        /* If this is in the reverse ( 3') transcriptome */
        if( first_read_is_rev_comp )
        {
            track_index = 1;
            promoter_pos = stop - 1;
        } 
        /* We are in the 5' ( positive ) transcriptome */
        else {
            track_index = 0;
            promoter_pos = start;
        }

        struct trace_segments_t* trace_segments
            = &(traces->segments[track_index][chr_index]);
        float scale_factor = cond_prob;
        /* Just update the promoter start position (n, n+1) */
        update_trace_segments_from_uniform_kernel(trace_segments, scale_factor,
            promoter_pos, promoter_pos+1);
    }
    
    return;
}

/* Returns the probability of observing the read, conditional on the trace */
struct update_mapped_read_rv_t 
update_CAGE_mapped_read_prbs( 
    struct cond_prbs_db_t* cond_prbs_db,
    const struct trace_t* const traces, 
    mapped_read_t* r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };

    /* FIXME cast away const? */
    mapped_read_index* rd_index;
    stack_allocate_and_init_mapped_read_index( rd_index, r );
    
    /* allocate space to store the temporary values */
    float* new_prbs = malloc( sizeof(float)*(rd_index->num_mappings) );
    
    /* Update the reads from the trace */
    double density_sum = 0;
    int i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        /* calculate the mean density */
        /* We set this to 2*DBL_EPSILON to prevent the division by 0 */
        double window_density = 2*DBL_EPSILON;

        mapped_read_location* current_loc = rd_index->mappings[i];

        MRL_CHR_TYPE chr_index =
            get_chr_from_mapped_read_location( current_loc );
        MRL_START_POS_TYPE start = 
            get_start_from_mapped_read_location( current_loc );
        MRL_START_POS_TYPE stop = 
            get_stop_from_mapped_read_location( current_loc ) ;
        assert( stop > 0 );
        enum bool is_paired 
            = mapped_read_location_is_paired( current_loc );
        enum bool first_read_is_rev_comp
            = first_read_in_mapped_read_location_is_rev_comp( current_loc );
        
        /* If the reads are paired */
        if( is_paired ) {
            statmap_log(LOG_FATAL, "FATAL: paired cage reads are not supported" );
        } 
        /* If the read is *not* paired */
        else {
            int track_index;
            unsigned int win_stop = 0;
            unsigned int win_start = 0;
            
            /* If this is in the rev stranded transcriptome */
            if( first_read_is_rev_comp )
            {
                track_index = 1;
                win_start = stop - MIN( stop, WINDOW_SIZE ) - 1;
                win_stop = MIN( traces->chr_lengths[chr_index], stop + WINDOW_SIZE - 1 );
            } 
            /* We are in the 5' ( positive ) transcriptome */
            else {
                track_index = 0;
                win_start = start - MIN( start, WINDOW_SIZE );
                win_stop = MIN( traces->chr_lengths[chr_index], start + WINDOW_SIZE );
            }

            window_density += accumulate_from_trace(traces, track_index,
                chr_index, win_start, win_stop);
        }
        
        new_prbs[i] = get_seq_error_from_mapped_read_location( current_loc )
                          *window_density;
        
        density_sum += new_prbs[i];
    }
    
    /* renormalize the read probabilities */
    if( rd_index->num_mappings > 0 )
    {
        assert( density_sum > 0 );
        
        for( i = 0; i < rd_index->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* quadratic error */
            // abs_error += pow(new_prbs[i]/density_sum 
            //                 - r.locations[i].cond_prob, 2 ) ;
            /* absolute value error */
            double normalized_density = new_prbs[i]/density_sum; 
            double old_cnd_prb =
                get_cond_prb( cond_prbs_db, rd_index->read_id, i );
            
            // printf("%i\t%e\t%e\t%e\n", i, new_prbs[i], normalized_density, old_cnd_prb );
            
            rv.max_change += fabs( normalized_density - old_cnd_prb );
           
            set_cond_prb( cond_prbs_db,
                          rd_index->read_id,
                          i,
                          normalized_density );
        }

        // printf( "%i\t%e\n", r->num_mappings, rv.max_change );
                    
        rv.log_lhd = log10( density_sum );
    }      

    free_mapped_read_index( rd_index );
    
    /* Free the tmp array */
    free( new_prbs );   

    return rv;
}

/*
 * END Cage Mapping Code
 *
 ******************************************************************************/

void update_ATACSeq_trace_expectation_from_location(
    const struct trace_t* const traces, 
    const mapped_read_location* const loc,
    float cond_prob
)
{
    MRL_CHR_TYPE chrm = 
        get_chr_from_mapped_read_location(loc);

    /* update the left open chromatin region */
    MRL_START_POS_TYPE start = get_start_from_mapped_read_location(loc) + 4;
    update_trace_segments_from_uniform_kernel(
        traces->segments[0] + chrm, 
        cond_prob, 
        start/BIN_SIZE, 
        start/BIN_SIZE+1 );
    
    MRL_STOP_POS_TYPE stop = get_stop_from_mapped_read_location(loc) - 5;
    update_trace_segments_from_uniform_kernel(
        traces->segments[0] + chrm, 
        cond_prob, 
        stop/BIN_SIZE, 
        stop/BIN_SIZE+1 );
    
    return;
}

struct update_mapped_read_rv_t 
update_ATACSeq_mapped_read_prbs(
    struct cond_prbs_db_t* cond_prbs_db, 
    const struct trace_t* const trace,
    mapped_read_t* r )
{
    struct update_mapped_read_rv_t rv = {0.0, 0.0};
    
    mapped_read_index* rd_index;
    stack_allocate_and_init_mapped_read_index(rd_index, r);
    
    /* store the log sequencing errors, and then maximum log sequencing
       error. We need the errors to normalize, and then max log errors
       to center the results so that we can avoid rounding errors */
    double* log_seq_errors = alloca(sizeof(double)*rd_index->num_mappings);
    double max_log_val = -ML_PRB_MAX;
    
    /* store the local read density for each mapping, to properly
       re-weight the sum */
    double* local_read_sum = alloca(sizeof(double)*rd_index->num_mappings);
    
    /* Update the trace from this mapping */        
    MPD_RD_ID_T i;
    for( i = 0; i < rd_index->num_mappings; i++ ) {
        mapped_read_location* loc = rd_index->mappings[i];

        ML_PRB_TYPE log_error =get_log_seq_error_from_mapped_read_location(loc);

        log_seq_errors[i] = log_error;
        max_log_val = MAX(log_error, max_log_val);

        MRL_START_POS_TYPE start = get_start_from_mapped_read_location(
            rd_index->mappings[i]);
        MRL_STOP_POS_TYPE stop = get_stop_from_mapped_read_location(
            rd_index->mappings[i]);

        local_read_sum[i] = accumulate_from_trace(
            trace,
            0,
            get_chr_from_mapped_read_location(rd_index->mappings[i]),
            (start-WINDOW_SIZE)/BIN_SIZE, (start+WINDOW_SIZE)/BIN_SIZE + 1);

        local_read_sum[i] += accumulate_from_trace(
            trace,
            0,
            get_chr_from_mapped_read_location(rd_index->mappings[i]),
            (stop-WINDOW_SIZE)/BIN_SIZE, (stop+WINDOW_SIZE)/BIN_SIZE + 1);
        
        assert(1 == isfinite(local_read_sum[i]));
    }

    double shifted_prb_sum = 0;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        double unnormalized_prb = local_read_sum[i]*pow(
            10, log_seq_errors[i]-max_log_val );
        local_read_sum[i] = unnormalized_prb;
        shifted_prb_sum += unnormalized_prb;
    }
    
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        ML_PRB_TYPE old_cond_prb = get_cond_prb( 
            cond_prbs_db, rd_index->read_id, i);
        ML_PRB_TYPE cond_prb = local_read_sum[i]/shifted_prb_sum;
        rv.max_change = MAX(rv.max_change, fabs(old_cond_prb - cond_prb));
        set_cond_prb( cond_prbs_db, rd_index->read_id, i, cond_prb);
        assert( (cond_prb <= 1 + 1e-6) && (cond_prb) >= 0-1e-6 );
    }
            
    free_mapped_read_index( rd_index );
    return rv;
}


/******************************************************************************
 *
 * Entry point(s) into the iterative mapping code
 *
 */

#if 0
struct trace_t*
build_segmented_trace(
        struct genome_data* genome,
        int num_tracks,
        char** track_names,
        struct mapped_reads_db* mpd_rdb,
        struct cond_prbs_db_t* cond_prbs_db,
        void (* const update_trace_expectation_from_location)(
            const struct trace_t* const traces, 
            const mapped_read_location* const loc,
            const float cond_prob )    )
{
    /* Initialize a full trace that we will use to determine the segments */
    struct trace_t* full_trace = NULL;
    init_full_trace( genome, &full_trace, num_tracks, track_names );
    reset_all_read_cond_probs( mpd_rdb, cond_prbs_db );
    update_traces_from_mapped_reads( mpd_rdb, cond_prbs_db, full_trace,
        update_trace_expectation_from_location );

    /* build list of segments from the full trace */
    struct segments_list *segments_list
        = build_trace_segments_list( full_trace );

    // DEBUG
    log_segments_list( segments_list );

    /* build the segmented trace graph and log it for debugging */
    build_segmented_trace_graph( segments_list, mpd_rdb );

    /* Initialize the segmented trace */
    struct trace_t* segmented_trace = NULL;
    init_trace( genome, &segmented_trace, num_tracks, track_names );

    /* Add the segments as specified by the segment list
       NOTE: this assumes segments in the segment list are sorted by start
       (they are naturally sorted with out current segment finding algo). */
    int i;
    for( i = 0; i < segments_list->length; i++ )
    {
        struct segment *s = segments_list->segments + i;

        /* get the list of segments to add to */
        struct trace_segments_t* update_segments
            = &(segmented_trace->segments[s->track_index][s->chr_index]);

        add_trace_segment_to_trace_segments( update_segments,
            s->track_index, s->chr_index, s->start, (s->stop - s->start) );
    }

    close_traces( full_trace );
    free_segments_list( segments_list );    

    return segmented_trace;
}
#endif

/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/

struct cond_prbs_db_t*
build_posterior_db( struct genome_data* genome, 
                    struct mapped_reads_db* mpd_rdb,
                    enum assay_type_t assay_type )
{   
    void (*update_traces_from_mapped_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float prb ) = NULL;

    struct update_mapped_read_rv_t 
        (*update_reads)( struct cond_prbs_db_t* cond_prbs_db, 
                         const struct trace_t* const traces, 
                         mapped_read_t* r  )
        = NULL;

    switch( assay_type )
    {
    case ATACSeq:
        update_traces_from_mapped_location 
            = update_ATACSeq_trace_expectation_from_location;
        update_reads = update_ATACSeq_mapped_read_prbs;
        break;

    case CAGE:
        update_traces_from_mapped_location 
            = update_CAGE_trace_expectation_from_location;
        update_reads = update_CAGE_mapped_read_prbs;
        break;
    
    case CHIP_SEQ:
        update_traces_from_mapped_location 
            = update_chipseq_trace_expectation_from_location;
        update_reads = update_chipseq_mapped_read_prbs;
        break;    

    // case STRANDED_RNASEQ:

    default:
        statmap_log(LOG_FATAL, "Unrecognized assay type in iterative mapping.");
        abort();
    }

    struct trace_t* traces = NULL;
    char* track_name = "fwd";
    init_binned_trace( genome, &traces, 1, &track_name, 1 );
    set_trace_to_uniform(traces, 1);
    
    struct cond_prbs_db_t* cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, mpd_rdb);
    reset_all_read_cond_probs( mpd_rdb, cond_prbs_db );
    
    update_mapping(
        mpd_rdb, cond_prbs_db, traces, 
        1e5, // max num iterations
        1e-2, // max_prb_change_for_convergence,
        1, // lhd_ratio_stop_value,
        update_traces_from_mapped_location,
        update_reads
        );
    
    write_trace_to_file(traces, "MLE.trace");
    close_traces(traces);
    return cond_prbs_db;
}
