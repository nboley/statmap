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

#include "log.h"

#define DONT_USE_VEC_OPERATIONS
#define LOCK_TRACES


#ifdef USE_VEC_OPERATIONS
typedef float v4sf __attribute__ ((vector_size(16) ));
                  
union f4vector 
{
    v4sf v;
    float f[4];
};
#endif

struct fragment_length_dist_t* global_fl_dist;
struct genome_data* global_genome;
struct trace_t* global_starting_trace;

/* This isn't being used...
static inline void
update_stranded_read_start_density_from_location( 
    const struct trace_t* const traces, 
    const mapped_read_location* const loc,
    const float cond_prob 
)
{
    const MRL_CHR_TYPE chr_index = get_chr_from_mapped_read_location( loc  );
    const MRL_START_POS_TYPE start = get_start_from_mapped_read_location( loc  );
    
    if( first_read_in_mapped_read_location_is_rev_comp( loc ) )
    {
        
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_lock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #endif
        
        
        // #pragma omp atomic
        traces->traces[1][chr_index][start] += cond_prob;
        // traces->traces[1][chr_index][start] += cond_prob;

        
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_unlock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #endif
        
    } else {
        
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_lock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #endif
        

        // #pragma omp atomic
        traces->traces[0][chr_index][start] += cond_prob;        
        // traces->traces[0][chr_index][start] += cond_prob;

        
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_unlock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #endif
        
    }
    
    return;
}
*/

/*
 * This is technically a non-parametric bootstrap
 * but since the reads are conditionally independent 
 * with multinomial reads densities, this is identical
 * to the parametric bootstrap without needing a bisection
 * step. The algorithm is
 * Repeat the following n_read times
 *     1) Choose a random read uniformily
 *     2) Choose a random location, weighted by it's probability
 *     3) update the trace
 */
void
bootstrap_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        float cond_prob
    )
)
{        
    global_fl_dist = reads_db->fl_dist;
    
    /* First, zero the trace to be summed into */
    zero_traces( traces );
    
    if( RAND_MAX < reads_db->num_mapped_reads )
    {
        statmap_log( LOG_ERROR, "cannot bootstrap random reads - number of reads exceeds RAND_MAX. ( This needs to be fixed. PLEASE file a bug report about this ) " );
        return;
    }

    int num_error_reads = 0;
    
    /* Update the trace from reads chosen randomly, with replacement */
    MPD_RD_ID_T i;
    for( i = 0; i < reads_db->num_mapped_reads; i++ )
    {
        MPD_RD_ID_T read_index = rand()%(reads_db->num_mapped_reads);
        assert( read_index < reads_db->num_mapped_reads );
        char* read_start = reads_db->index[ read_index ].ptr;

        mapped_read_t* rd = (mapped_read_t*) read_start;

        mapped_read_index* rd_index = NULL;
        init_mapped_read_index( &rd_index, rd );

        if( i > 0 && i%1000000 == 0 )
        {
            statmap_log( LOG_NOTICE, "Read %i reads in bootstrap",  i  );
        }

        /* If there are no mappings, then this is pointless */
        if( rd_index->num_mappings > 0 )
        {
            /* Choose a random location, proportional to the normalized probabilities */
            MPD_RD_ID_T j = 0;
            /* if there is exactly 1 mapping, then the randomness is pointless */
            if( rd_index->num_mappings > 1 )
            {
                float random_num = (float)rand()/(float)RAND_MAX;
                double cum_dist = 0;
                for( j = 0; j < rd_index->num_mappings; j++ ) 
                {
                    cum_dist += get_cond_prb(
                            cond_prbs_db, rd_index->read_id, j );
                    if( random_num <= cum_dist )
                        break;
                }

                /* 
                 * rounding error makes it possible for j == r.num_mappings. If this is 
                 * the case, then make sure the differnece isn't too big, and then chnage 
                 * it. We could save a few cycles by taking this into account earlier, but
                 * for now we're going to leave the sanity check.
                 */
                assert( j <= rd_index->num_mappings );
                if( j == rd_index->num_mappings )
                {
                    j = rd_index->num_mappings - 1;
                    if( cum_dist < 0.999 )
                    {
                        /* BUG BUG BUG - I dont think that this should ever happen */
                        num_error_reads += 1;
                        // fprintf( stderr, "%e\t%li\t%i\n", cum_dist, r.read_id, get_start_from_mapped_read_location( r.locations + j ) );
                        // fprintf( stderr, "ERROR      : There is a serious error in the bs code. It appears as if the read cond probs dont sum to 1 - contuining but PLEASE file a bug report.\n");
                        continue;
                    }
                }
            }
                        
            /* Add this location to the trace, with prb 1 */
            update_trace_expectation_from_location( traces,
                                                    rd_index->mappings[j],
                                                    1.0 );
        }

        free_mapped_read_index( rd_index );
    }
    
    if( num_error_reads > 0 )
    {
        statmap_log( LOG_ERROR,
                "%i reads ( out of %u ) had errors in which the cum dist didnt add up to 1.",
                num_error_reads,
                reads_db->num_mapped_reads
            );
    }
    
    return;
}

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

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) )     
    {
        mapped_read_index* rd_index = NULL;
        init_mapped_read_index( &rd_index, r );

        /* Update the trace from this mapping */        
        MPD_RD_ID_T j;
        for( j = 0; j < rd_index->num_mappings; j++ ) {
            float cond_prob = get_cond_prb(
                cond_prbs_db, rd_index->read_id, j );
            
            update_trace_expectation_from_location( 
                traces, rd_index->mappings[j], cond_prob );
        }
        
        free_mapped_read_index( rd_index );
    }
    
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
        while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) )     
        {
            read_num++;
            
            mapped_read_index* rd_index = NULL;
            init_mapped_read_index( &rd_index, r );

            /* Update the trace from this mapping */        
            MPD_RD_ID_T j;
            for( j = 0; j < rd_index->num_mappings; j++ ) {
                float cond_prob = get_cond_prb(
                    cond_prbs_db, rd_index->read_id, j );
                
                update_trace_expectation_from_location( 
                    traces, rd_index->mappings[j], cond_prob );
            }
            
            free_mapped_read_index( rd_index );
        }
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
     * I like this better than a block of reads approach, because it ensures that 
     * any io is sequential. Of course, this is at the cost of some lock 
     * contention, but that should be minor, espcially if ( the mmapped ) rdb
     * can't fully fit into memory.
     */
    pthread_spinlock_t* curr_read_index_spinlock;
    MPD_RD_ID_T* curr_read_index;
    
    struct mapped_reads_db* reads_db;
    struct cond_prbs_db_t* cond_prbs_db;
    
    struct trace_t* traces;
    struct update_mapped_read_rv_t 
        (* update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                     const struct trace_t* const traces, 
                                     const mapped_read_t* const r  );
    
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
                                     const mapped_read_t* const r  )
        = ( (struct update_mapped_reads_param*)
            params)->update_mapped_read_prbs;

    mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) ) 
    {
        /* Update the read */
        struct update_mapped_read_rv_t tmp_rv 
            = update_mapped_read_prbs( cond_prbs_db, traces, r );

        /* Hand reduce the max */
        if( tmp_rv.max_change > ( (struct update_mapped_reads_param*) params)->rv.max_change )
            ( (struct update_mapped_reads_param*) params)->rv.max_change = tmp_rv.max_change;
        
        /* Update the lhd */
        ( (struct update_mapped_reads_param*) params)->rv.log_lhd += tmp_rv.log_lhd;
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
                                           const mapped_read_t* const r  )
    )
{    
    /* store the toal accumulated error */
    struct update_mapped_read_rv_t rv = { 0, 0.0 };
    
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
                statmap_log( LOG_FATAL, "Return code from pthread_create() is %d", rc );
            }
        }

        /* wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            pthread_attr_destroy( attrs + t );
            
            if (rc) {
                statmap_log( LOG_FATAL, "Return code from pthread_join() is %d", rc );
            }
        }
        
        
        /* Aggregate over the max in the error difference */
        for( t=0; t < num_threads; t++ )
        {
            rv.log_lhd += tds[t].rv.log_lhd;
            rv.max_change = MAX( tds[t].rv.max_change, rv.max_change );
        }
        
        free( tds );
        free( attrs );
    }
    
    return rv;
}

double
calc_log_lhd( 
    struct mapped_reads_db* reads_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           const mapped_read_t* const r  )
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
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           const mapped_read_t* const r  )
    )
{
    /* make sure the mapped reads are open for reading */
    assert( rdb->mode == 'r' );

    double prev_log_lhd = 0;
    
    struct update_mapped_read_rv_t rv;
    
    int num_iterations = 0;
    for( num_iterations = 0; 
         num_iterations < max_num_iterations; 
         num_iterations++ )
    {
        /* Normalize the trace sum to 1 */
        /* This makes the trace the current marginal read density estimate */
        struct timeval nt_start, nt_stop, umr_start, umr_stop, ut_start, ut_stop;;
        
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
                  || num_iterations%1 == 0 
                  || (ut_stop.tv_sec-nt_start.tv_sec) > 30 
                )
                || rv.max_change < max_prb_change_for_convergence
                || ( lhd_ratio_stop_value >= 1.0
                     && pow( 10, rv.log_lhd - prev_log_lhd ) < lhd_ratio_stop_value
                     && pow( 10, rv.log_lhd - prev_log_lhd ) > 0.95
                   )
                )
            )
        {
            statmap_log( LOG_INFO,
                    "Iter %i: \tError: %e \tLog Lhd: %e (ratio %e) \tNorm Trace:  %.5f sec\t Read UT:  %.5f sec\tTrace UT:  %.5f sec", 
                    num_iterations, rv.max_change, rv.log_lhd,
                    pow( 10, rv.log_lhd - prev_log_lhd ),
                    (float)(nt_stop.tv_sec - nt_start.tv_sec) + ((float)(nt_stop.tv_usec - nt_start.tv_usec))/1000000,
                    (float)(umr_stop.tv_sec - umr_start.tv_sec) + ((float)(umr_stop.tv_usec - umr_start.tv_usec))/1000000,
                    (float)(ut_stop.tv_sec - ut_start.tv_sec) + ((float)(ut_stop.tv_usec - ut_start.tv_usec))/1000000
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

void
build_random_starting_trace( 
    struct trace_t* traces, 
    struct genome_data* genome,

    struct mapped_reads_db* rdb,
    struct cond_prbs_db_t* cond_prbs_db,

    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob),

    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           const mapped_read_t* const r  )
    )
{
    global_fl_dist = rdb->fl_dist;
    global_genome = genome;
    
    /* int traces to a low consant. This should heavily favor prior reads. */
    set_trace_to_uniform( traces, EXPLORATION_PRIOR );

    rewind_mapped_reads_db( rdb );
    mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) ) 
    {        
        /* reset the read cond prbs under a uniform prior */
        reset_read_cond_probs( cond_prbs_db, r, rdb );
        
        /* update the read conditional probabilities from the trace */
        update_mapped_read_prbs( cond_prbs_db, traces, r );        

        mapped_read_index* rd_index = NULL;
        init_mapped_read_index( &rd_index, r );
        
        if( rd_index->num_mappings > 0 )
        {
            /* Choose a random loc, proportional to the normalized probabilities */
            int j = 0;
            /* if there is exactly 1 mapping, then the randomness is pointless */
            if( rd_index->num_mappings > 1 )
            {
                float random_num = (float)rand()/(float)RAND_MAX;
                double cum_dist = 0;
                for( j = 0; j < rd_index->num_mappings; j++ ) 
                {
                    cum_dist += get_cond_prb(
                            cond_prbs_db, rd_index->read_id, j );
                    if( random_num <= cum_dist )
                        break;
                }

                /* 
                 * rounding error makes it possible for j == r.num_mappings. If this is 
                 * the case, then make sure the differnece isn't too big, and then chnage 
                 * it. We could save a few cycles by taking this into account earlier, but
                 * for now we're going to leave the sanity check.
                 */
                assert( j <= rd_index->num_mappings );
                if( j == rd_index->num_mappings )
                {
                    // printf( "WARNING - rounding error! %e\n", cum_dist );
                    j = rd_index->num_mappings - 1;
                    assert( cum_dist > 0.999 );
                }
            }
                        
            /* Add this location to the trace, with prb 1 */
            update_trace_expectation_from_location(
                    traces, rd_index->mappings[j], 1.0 );
        }

        free_mapped_read_index( rd_index );
    }
    
    return;
}

/*
 * END Chr Trace Code
 *
 *****************************************************************************/


/******************************************************************************
 *
 * High-Level Functions 
 *
 *****************************************************************************/

int 
sample_random_trace(
    struct mapped_reads_db* rdb, 
    struct cond_prbs_db_t* cond_prbs_db,
    struct genome_data* genome,

    int num_tracks,
    char** track_names,

    int sample_index,
    
    // (starting) sample meta info
    FILE* ss_mi,
    FILE* s_mi,

    int max_num_iterations,
    float max_prb_change_for_convergence,
    float min_lhd_ratio_change_for_convergence,

    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const mapped_read_location* const loc,
        const float cond_prob
    ),
    
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( struct cond_prbs_db_t* cond_prbs_db,
                                           const struct trace_t* const traces, 
                                           const mapped_read_t* const r  )
    )
{
    statmap_log( LOG_INFO, "Starting sample %i", sample_index+1 );

    struct trace_t* sample_trace;
    init_trace( genome, &sample_trace, num_tracks, track_names );

    build_random_starting_trace( 
        sample_trace, genome, rdb, cond_prbs_db,
        update_trace_expectation_from_location,
        update_mapped_read_prbs
     );

    if( SAVE_STARTING_SAMPLES )
    {
        char buffer[100];
        sprintf( buffer, "%ssample%i.bin.trace", STARTING_SAMPLES_PATH, sample_index+1 );
            
        write_trace_to_file( sample_trace, buffer );

        double log_lhd = calc_log_lhd( 
            rdb, cond_prbs_db, sample_trace, update_mapped_read_prbs );
        fprintf( ss_mi, "%i,%e\n", sample_index+1, log_lhd );
        fflush( ss_mi );
    }
        
    /* maximize the likelihood */
    update_mapping( 
        rdb, cond_prbs_db, sample_trace, max_num_iterations,
        max_prb_change_for_convergence,
        min_lhd_ratio_change_for_convergence, 
        update_trace_expectation_from_location,
        update_mapped_read_prbs
    );

    if( SAVE_SAMPLES )
    {
        char buffer[100];
        sprintf( buffer, "%ssample%i.bin.trace", RELAXED_SAMPLES_PATH, sample_index+1 );

        write_trace_to_file( sample_trace, buffer );
            
        double log_lhd = calc_log_lhd( 
            rdb, cond_prbs_db, sample_trace, update_mapped_read_prbs );
        fprintf( s_mi, "%i,%e\n", sample_index+1, log_lhd );
        fflush( s_mi );
    }
        
    close_traces( sample_trace );
    
    return 0;
}

/*
 *
 * END High-Level Functions 
 *
 ******************************************************************************/


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
    const MRL_START_POS_TYPE start = get_start_from_mapped_read_location( loc  );
    const MRL_STOP_POS_TYPE stop = get_stop_from_mapped_read_location( loc  );

    enum bool is_paired 
        = mapped_read_location_is_paired( loc );
    enum bool first_read_is_rev_comp
        = first_read_in_mapped_read_location_is_rev_comp( loc );
    
    assert( cond_prob >= 0.0 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_chrs );
    assert( traces->chr_lengths[chr_index] >= stop );
    
    /* iteration variable */
    unsigned int k;

    /* update the trace */
    /* If the reads are paired */
    if( is_paired )
    {
        const float scale = (1.0/(stop-start));

        int trace_index;
        if( first_read_is_rev_comp )
        {
            trace_index = 0;
        } else {
            trace_index = 1;
        }
        
        for( k = start/TM_GRAN; k <= (stop-1)/TM_GRAN; k++ )
        {
            /* lock the spinlocks */
            #ifdef LOCK_TRACES
                #ifdef USE_MUTEX
                pthread_mutex_lock( traces->locks[trace_index][ chr_index ] + k );
                #else
                pthread_spin_lock( traces->locks[trace_index][ chr_index ] + k );
                #endif
            #endif
        
            unsigned int LR_start = MAX( start, k*TM_GRAN ); // Locked Region start
            unsigned int LR_stop = MIN( stop, (k+1)*TM_GRAN ); // Locked Region stop
            
            unsigned int j;
            for( j = LR_start; j < LR_stop; j++ )
            {
                assert( chr_index < traces->num_chrs );
                assert( j < traces->chr_lengths[chr_index] );
                
                traces->traces[trace_index][chr_index][j] 
                    += scale*cond_prob;
            }

            /* unlock the spinlocks */
            #ifdef LOCK_TRACES
                #ifdef USE_MUTEX
                pthread_mutex_unlock( traces->locks[trace_index][ chr_index ] + k );
                #else
                pthread_spin_unlock( traces->locks[trace_index][ chr_index ] + k );
                #endif
            #endif
        }
    } 
    /* If the read is *not* paired */
    else {
        /* the binding site must be on the five prime side of the read */
        unsigned int window_start;
        unsigned int window_stop;
        /* determine which trace to lock */
        int trace_index;
        if( first_read_is_rev_comp )
        {
            /* prevent overrunning the trace bounadry - this is fancy to avoid the 
               fact that the indexes are unsigned */
            window_start = stop - MIN( stop, (unsigned int) global_fl_dist->max_fl );
            window_stop = MIN( stop, traces->chr_lengths[chr_index] );
            trace_index = 1;
        } else {
            window_start = start;
            window_stop = MIN( start + global_fl_dist->max_fl, 
                               traces->chr_lengths[chr_index] );
            trace_index = 0;
        }

        for( k = window_start/TM_GRAN; k <= (window_stop-1)/TM_GRAN; k++ )
        {
            int LR_start = MAX( window_start, k*TM_GRAN ); // Locked Region start
            int LR_stop = MIN( window_stop, (k+1)*TM_GRAN ); // Locked Region stop
            
            /* lock the spinlocks */
            #ifdef LOCK_TRACES
                #ifdef USE_MUTEX
                pthread_mutex_lock( traces->locks[trace_index][ chr_index ] + k );
                #else
                pthread_spin_lock( traces->locks[trace_index][ chr_index ] + k );
                #endif
            #endif
            
            float* fl_dist_array;
            int trace_index;
            if( first_read_is_rev_comp )
            {
                /* we use the rev_chipseq_bs_density instead of just
                   moving through the bs density in reverse order so
                   that we can use vector operations */ 
                /* note that, even if the vec operations are disabled,
                   we want this so that gcc can auto-vectorize */
                fl_dist_array = global_fl_dist->rev_chipseq_bs_density;
                trace_index = 1;
            } else {
                fl_dist_array = global_fl_dist->chipseq_bs_density;
                trace_index = 0;
            }            
            
            int j;
            #ifndef USE_VEC_OPERATIONS
            for( j = LR_start; j < LR_stop; j++ )
            {
                float fl_density = fl_dist_array[ j - LR_start ];
                fl_density *= cond_prob;
                traces->traces[trace_index][chr_index][j] += fl_density;
            }
            #else
            
            #define AL_SIZE 16
            union f4vector cond_prob_vec;
            cond_prob_vec.f[0] = cond_prob;
            cond_prob_vec.f[1] = cond_prob;
            cond_prob_vec.f[2] = cond_prob;
            cond_prob_vec.f[3] = cond_prob;

            union f4vector* fl_array = (union f4vector*) fl_dist_array;
            /* make sure this is aligned to the 16 byte boundary */
            assert( ( ( (size_t)fl_dist_array)%AL_SIZE ) == 0 );
            
            union f4vector* trace_array = 
                (union f4vector*) ( traces->traces[trace_index][chr_index] + LR_start );
            /* word align the trace array */
            /* save this so we know what has been skipped */
            size_t old_trace_array_start = (size_t) trace_array;
            trace_array =  (union f4vector*) ( (size_t)trace_array + ( AL_SIZE - ((size_t)trace_array)%AL_SIZE ) );
            
            /* deal with the unaligned floats */
            for( j = 0; j < (int) ( (size_t)trace_array - old_trace_array_start ); j++ )
            {
                float fl_density = fl_dist_array[ j ];
                fl_density *= cond_prob;
                traces->traces[trace_index][chr_index][j+LR_start] += fl_density;
            }
            
            /* we subtract 1 to account for the revised alignment */
            for( j = 0; j < (LR_stop-LR_start)/4 - 1; j++ )
            {
                union f4vector fl_density __attribute__ ((aligned (128)));
                fl_density.v = __builtin_ia32_mulps( fl_array[j].v, cond_prob_vec.v);
                trace_array[j].v += fl_density.v;
            }
            
            /* deal with the excess - its always greater than 1 
               because of the alignment floats*/
            for( j = LR_stop - LR_stop%4 - 4; j < LR_stop; j++ )
            {
                float fl_density = fl_dist_array[ j-LR_start ];
                fl_density *= cond_prob;
                traces->traces[trace_index][chr_index][j] += fl_density;
                assert( fl_density < 1 );
            }
            #endif
        
            /* unlock the spinlocks */
            #ifdef LOCK_TRACES
                #ifdef USE_MUTEX
                pthread_mutex_unlock( traces->locks[trace_index][ chr_index ] + k );
                #else
                pthread_spin_unlock( traces->locks[trace_index][ chr_index ] + k );
                #endif
            #endif
        }        
    }

    return;
}

struct update_mapped_read_rv_t 
update_chipseq_mapped_read_prbs(     struct cond_prbs_db_t* cond_prbs_db,
                                     const struct trace_t* const traces, 
                                     const mapped_read_t* const r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };

    mapped_read_index* rd_index;
    /* FIXME cast away const? */
    init_mapped_read_index( &rd_index, (mapped_read_t*) r );
    
    /* allocate space to store the temporary values */
    float* new_prbs = malloc( sizeof(float)*(rd_index->num_mappings) );
    
    /* Update the reads from the trace */
    /* set to flt_min to avoid division by 0 */
    double prb_sum = ML_PRB_MIN;
    double error_prb_sum = ML_PRB_MIN;

    int i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        /* calculate the mean density */
        double window_density = 0;

        mapped_read_location* current_loc = rd_index->mappings[i];

        MRL_CHR_TYPE chr_index =
            get_chr_from_mapped_read_location( current_loc );
        MRL_START_POS_TYPE start =
            get_start_from_mapped_read_location( current_loc );
        MRL_STOP_POS_TYPE stop =
            get_stop_from_mapped_read_location( current_loc );

        enum bool is_paired 
            = mapped_read_location_is_paired( current_loc );
        enum bool first_read_is_rev_comp
            = first_read_in_mapped_read_location_is_rev_comp( current_loc );

        /* If the reads are paired */
        unsigned int j = 0;
        if( is_paired )
        {
            for( j = start; j <= stop; j++ )
            {
                assert( j < traces->chr_lengths[chr_index] );
                window_density += traces->traces[0][chr_index][j]
                                  + traces->traces[1][chr_index][j];
            }

            /* 
             * This is the probability of observing the seqeunce given that it came from 
             * location i, assuming the chip traces have been normalized to 0. 
             */
            float fl_prb = get_fl_prb( 
                global_fl_dist, 
                get_fl_from_mapped_read_location( current_loc ) 
            );
            
            window_density *= fl_prb;
        } 
        /* If the read is *not* paired */
        else {
            unsigned int k;
            
            /* the binding site must be on the five prime side of the read */
            unsigned int window_start;
            unsigned int window_stop;
            /* determine which trace to use */
            if( first_read_is_rev_comp )
            {
                /* we don't use MAX bc stop-max_fl might be negative */
                window_start = stop - MIN( stop, (unsigned int) global_fl_dist->max_fl );
                window_stop = MIN( stop, global_genome->chr_lens[chr_index] );
            } else {
                window_start = start;
                window_stop = MIN( start + global_fl_dist->max_fl, 
                                   global_genome->chr_lens[chr_index] );
            }

            if( first_read_is_rev_comp )
            {
                for( k = window_start; k < window_stop; k++ )
                {
                    /* 
                       This is a bit confusing. First, at this stage we dont care 
                       whether the inferred binding site came from 5'or 3' reads  
                       - so we add up the binding site density from both strands. 
                       Next, we normalize this value by the precomputed frag length
                       normalization factor.
                    */
                    float value =  traces->traces[0][chr_index][k];
                    value += traces->traces[1][chr_index][k];
                    value *= global_fl_dist->chipseq_bs_density[ 
                                 global_fl_dist->max_fl - 1 - (k-window_start) ];
                    window_density += value;
                }
            } else {
                for( k = window_start; k < window_stop; k++ )
                {
                    /* same comment as directly above */
                    float value = traces->traces[0][chr_index][k];
                    value += traces->traces[1][chr_index][k]; 
                    value *= global_fl_dist->chipseq_bs_density[ k-window_start ];
                    window_density += value;
                }
            }
        }
        
        new_prbs[i] = 
            get_seq_error_from_mapped_read_location( current_loc )*
            window_density;
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
    if( is_paired )
    {
        statmap_log( LOG_FATAL, "paired cage reads are not supported" );
        exit( -1 );                
    } 
    /* If the read is *not* paired */
    else {
        int trace_index;
        int promoter_pos = -1;
        /* If this is in the reverse ( 3') transcriptome */
        if( first_read_is_rev_comp )
        {
            trace_index = 1;
            promoter_pos = stop - 1;
        } 
        /* We are in the 5' ( positive ) transcriptome */
        else {
            trace_index = 0;
            promoter_pos = start;
        }
        
        /* lock the spinlock */
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[ trace_index ][ chr_index ] + (promoter_pos/TM_GRAN) );
        #else
        pthread_spin_lock( traces->locks[ trace_index ][ chr_index ] + (promoter_pos/TM_GRAN) );
        #endif        
        traces->traces[ trace_index ][ chr_index ][ promoter_pos ] += cond_prob; 
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[ trace_index ][ chr_index ] + (promoter_pos/TM_GRAN) );
        #else
        pthread_spin_unlock( traces->locks[ trace_index ][ chr_index ] + (promoter_pos/TM_GRAN) );
        #endif
    }
    
    return;
}

/* Returns the probability of observing the read, conditional on the trace */
inline struct update_mapped_read_rv_t 
update_CAGE_mapped_read_prbs( 
    struct cond_prbs_db_t* cond_prbs_db,
    const struct trace_t* const traces, 
    const mapped_read_t* const r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };

    mapped_read_index* rd_index;
    /* FIXME cast away const? */
    init_mapped_read_index( &rd_index, (mapped_read_t*) r );
    
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
        unsigned int j = 0;
        if( is_paired )
        {
            statmap_log( LOG_FATAL, "paired cage reads are not supported" );
            continue;
        } 
        /* If the read is *not* paired */
        else {
            /* store the trace that we care about */
            TRACE_TYPE** trace;
            unsigned int win_stop = 0;
            unsigned int win_start = 0;
            
            /* If this is in the rev stranded transcriptome */
            if( first_read_is_rev_comp )
            {
                trace = traces->traces[1];
                win_start = stop - MIN( stop, WINDOW_SIZE ) - 1;
                win_stop = MIN( traces->chr_lengths[chr_index], stop + WINDOW_SIZE - 1 );
            } 
            /* We are in the 5' ( positive ) transcriptome */
            else {
                trace = traces->traces[0];
                win_start = start - MIN( start, WINDOW_SIZE );
                win_stop = MIN( traces->chr_lengths[chr_index], start + WINDOW_SIZE );
            }
            
            for( j = win_start; j < win_stop; j++ )
            {
                window_density += trace[chr_index][j];
            }
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

/******************************************************************************
 *
 * Entry point(s) into the iterative mapping code
 *
 */


int
update_chipseq_mapping_wnc(  
    int sample_index,
    
    struct mapped_reads_db* ip_rdb, 
    struct cond_prbs_db_t* ip_cond_prbs_db,
    struct trace_t** ip_trace,
    
    struct mapped_reads_db* nc_rdb,
    struct cond_prbs_db_t* nc_cond_prbs_db,
    struct trace_t** nc_trace,
    
    struct genome_data* genome,

    float max_prb_change_for_convergence,
    /* true if we should use a random start - otherwise, we use uniform */
    enum bool random_start )
{
    struct timeval start, stop;
    
    int error = 0;
    
    /* BUG!!! */
    /* Set the global fl dist - we should probably be passing this in through
       the call chain via a void ptr, but thats a lot of work with very little
       benefit ( except some cleanliness of course )
     */
    global_fl_dist = ip_rdb->fl_dist;
    global_genome = genome;

    /* reset the read cond prbs under a uniform prior */
    reset_all_read_cond_probs( ip_rdb, ip_cond_prbs_db );
    reset_all_read_cond_probs( nc_rdb, nc_cond_prbs_db );
            
    /* iteratively map from a uniform prior */
    gettimeofday( &start, NULL );
    statmap_log( LOG_NOTICE, "Starting iterative mapping." );
    
    /* initialize the trace that we will store the expectation in */
    /* it has dimension 2 - for the negative and positive stranded reads */
    char* ip_track_names[2] = {"IP_fwd_strand", "IP_bkwd_strand"};
    init_trace( genome, ip_trace, 2, ip_track_names );
    char* nc_track_names[2] = {"NC_fwd_strand", "NC_bkwd_strand"};
    init_trace( genome, nc_trace, 2, nc_track_names );    

    if( false == random_start )
    {
        set_trace_to_uniform( *ip_trace, 1 );    
    } else {
        build_random_starting_trace( 
            *ip_trace, genome, 
            ip_rdb, ip_cond_prbs_db,
            update_chipseq_trace_expectation_from_location,
            update_chipseq_mapped_read_prbs
        );

        if( SAVE_STARTING_SAMPLES )
        {
            char buffer[100];
            sprintf( buffer, "%ssample%i.bin.trace", STARTING_SAMPLES_PATH, sample_index+1 );
            
            write_trace_to_file( *ip_trace, buffer );
        }
    }
    
    /* reset the read cond prbs under a uniform prior */
    reset_all_read_cond_probs( ip_rdb, ip_cond_prbs_db );
    
    /* update the 'real' IP data */
    error = update_mapping (
        ip_rdb, 
        ip_cond_prbs_db,
        *ip_trace,
        MAX_NUM_EM_ITERATIONS,
        max_prb_change_for_convergence,
        CHIPSEQ_LHD_RATIO_STOP_VALUE,
        update_chipseq_trace_expectation_from_location,
        update_chipseq_mapped_read_prbs
    );

    if( 0 != error ) {
        statmap_log( LOG_FATAL, "Unrecognized error code %i from update_mapping", error );
    }
    
    gettimeofday( &stop, NULL );
    long seconds  = stop.tv_sec  - start.tv_sec;
    long useconds = stop.tv_usec - start.tv_usec;
    
    statmap_log( LOG_INFO, "Maximized LHD in %.5f seconds\n",
            ((float)seconds) + ((float)useconds)/1000000
        );
    
    /* update the NC data based upon the IP data's NC */
    normalize_traces( *ip_trace );
    
    update_mapped_reads_from_trace(
            nc_rdb, nc_cond_prbs_db,
            *ip_trace, update_chipseq_mapped_read_prbs
        );

    /* update the NC reads from the ip marginal density */
    update_traces_from_mapped_reads( 
        nc_rdb, nc_cond_prbs_db, *nc_trace, 
        update_chipseq_trace_expectation_from_location
    );

    /* we need to do this because I renormalized the read sensity, 
       and it screws up the wiggle printing. Kind of a hack... */
    update_traces_from_mapped_reads( 
        ip_rdb, ip_cond_prbs_db, *ip_trace, 
        update_chipseq_trace_expectation_from_location
    );
    
    /* 
       now we have two traces - one is the NC and one is the IP. They are 
       mapped from the same marginal read density, so we should be able
       to call peaks on them individually 
    */
    goto cleanup;
    
cleanup:
    
    return 0;    
}

void
take_chipseq_sample_wnc(
    struct mapped_reads_db* chip_mpd_rds_db, 
    struct mapped_reads_db* NC_mpd_rds_db,
    struct genome_data* genome,
    
    FILE* meta_info_fp,
    int sample_index,
    
    float max_prb_change_for_convergence,
    /* true if we should use a random start - otherwise, we use uniform */
    enum bool random_start )

{
    global_fl_dist = chip_mpd_rds_db->fl_dist;
    assert( global_fl_dist != NULL );

    /* traces and track names to store the 2 marginal densities */
    struct trace_t* ip_trace;
    struct trace_t* nc_trace;

    struct cond_prbs_db_t* chip_cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &chip_cond_prbs_db, chip_mpd_rds_db );
    
    struct cond_prbs_db_t* NC_cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &NC_cond_prbs_db, NC_mpd_rds_db );
    
    /* jointly update the mappings */
    update_chipseq_mapping_wnc(  sample_index,
                                 chip_mpd_rds_db, chip_cond_prbs_db, &ip_trace,
                                 NC_mpd_rds_db, NC_cond_prbs_db, &nc_trace,
                                 genome, 
                                 max_prb_change_for_convergence,
                                 random_start ); 
           
    /* write the joint mappings to a single wiggle */
    /* we use a trick - writing wiggles appends to the file. So we just
       call the writing code with the same filename, and let it append 
       automatically 
    */
    if( SAVE_SAMPLES )
    {
        char wig_fname[100];
        sprintf( wig_fname, "%ssample%i.ip.bin.trace", RELAXED_SAMPLES_PATH, sample_index+1 );
        write_trace_to_file( ip_trace, wig_fname );
                
        sprintf( wig_fname, "%ssample%i.nc.bin.trace", RELAXED_SAMPLES_PATH, sample_index+1 );
        write_trace_to_file( ip_trace, wig_fname );

        /* write the lhd to the meta info folder */
        double log_lhd = calc_log_lhd( 
            chip_mpd_rds_db, chip_cond_prbs_db, 
            ip_trace, update_chipseq_mapped_read_prbs );
        fprintf( meta_info_fp, "%i,%e\n", sample_index+1, log_lhd );
        fflush( meta_info_fp );
    }
            
    if( NUM_BOOTSTRAP_SAMPLES > 0 )
    {
        char buffer[200];
        sprintf( buffer, "mkdir %ssample%i/",
                 BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1 );
        int error = system( buffer );
        if (WIFSIGNALED(error) &&
            (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            statmap_log( LOG_FATAL, "System call '%s' failed", buffer );
        }

        int j;
        for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
        {                
            if( j % (NUM_BOOTSTRAP_SAMPLES/10) == 0 )
                fprintf( stderr, " %.1f%%...", (100.0*j)/NUM_BOOTSTRAP_SAMPLES );

            /* actually perform the bootstrap */
            bootstrap_traces_from_mapped_reads( 
                chip_mpd_rds_db, chip_cond_prbs_db, ip_trace, 
                update_chipseq_trace_expectation_from_location );
                    
            /* actually perform the bootstrap */
            bootstrap_traces_from_mapped_reads( 
                NC_mpd_rds_db, NC_cond_prbs_db, nc_trace,
                update_chipseq_trace_expectation_from_location );
                    
            if( SAVE_BOOTSTRAP_SAMPLES )
            {
                /* write out the IP */
                sprintf( buffer, "%ssample%i/bssample%i.ip.bin.trace", 
                         BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1, j+1 );

                write_trace_to_file( ip_trace, buffer );

                /* write out the NC */
                sprintf( buffer, "%ssample%i/bssample%i.nc.bin.trace", 
                         BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1, j+1 );

                write_trace_to_file( nc_trace, buffer );                        
            }    
        }
        fprintf( stderr, " 100%%\n" );
                
        close_traces( ip_trace );
        close_traces( nc_trace );
    }     
    
    free_cond_prbs_db( chip_cond_prbs_db );
    free_cond_prbs_db( NC_cond_prbs_db );
    
    return;
}

int
take_chipseq_sample( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    
    int sample_index,
    
    FILE* ss_mi,
    FILE* s_mi,
    
    int max_num_iterations,
    float max_prb_change_for_convergence
) {
    int trace_size = 2;

    char** track_names = NULL;
    track_names = malloc( trace_size*sizeof(char*) );
    track_names[0] = "fwd_strnd_read_density"; 
    track_names[1] = "rev_strnd_read_density";

    struct cond_prbs_db_t* cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, rdb );

    sample_random_trace(
        rdb, cond_prbs_db, genome, trace_size, track_names, sample_index,
        ss_mi, s_mi, 
        max_num_iterations, 
        max_prb_change_for_convergence,
        CHIPSEQ_LHD_RATIO_STOP_VALUE,
        update_chipseq_trace_expectation_from_location,
        update_chipseq_mapped_read_prbs
    );
    
    free( track_names );
    free_cond_prbs_db( cond_prbs_db );

    return 0;
}

int
take_cage_sample( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    
    int sample_index,
    
    FILE* ss_mi,
    FILE* s_mi,
    
    int max_num_iterations,
    double max_prb_change_for_convergence
) {
    int trace_size = 2;

    char** track_names = NULL;
    track_names = malloc( trace_size*sizeof(char*) );
    track_names[0] = "fwd_strnd_read_density"; 
    track_names[1] = "rev_strnd_read_density";

    struct cond_prbs_db_t* cond_prbs_db;
    init_cond_prbs_db_from_mpd_rdb( &cond_prbs_db, rdb );

    sample_random_trace(
        rdb, cond_prbs_db, genome, trace_size, track_names, sample_index,
        ss_mi, s_mi, 
        max_num_iterations, 
        max_prb_change_for_convergence,
        CAGE_LHD_RATIO_STOP_VALUE,
        update_CAGE_trace_expectation_from_location,
        update_CAGE_mapped_read_prbs
    );

    free( track_names );
    free_cond_prbs_db( cond_prbs_db );
    
    return 0;
}


int
generic_update_mapping(  struct mapped_reads_db* rdb, 
                         struct genome_data* genome,
                         enum assay_type_t assay_type,
                         int num_samples,
                         float max_prb_change_for_convergence)
{
    /* Set the global fl dist */
    global_fl_dist = rdb->fl_dist;
    global_genome = genome;

    /* Store meta information about the samples and starting samples */
    FILE* ss_mi = NULL;
    FILE* s_mi = NULL;
    
    /* Create the meta data csv's */
    if( SAVE_STARTING_SAMPLES )
    {
        ss_mi = fopen(STARTING_SAMPLES_META_INFO_FNAME, "w");
        fprintf( ss_mi, "sample_number,log_lhd\n" );
    }

    s_mi = fopen(RELAXED_SAMPLES_META_INFO_FNAME, "w");
    fprintf( s_mi, "sample_number,log_lhd\n" );
    
    /* build bootstrap samples. Then take the max and min. */
    int i;
    for( i = 0; i < num_samples; i++ )
    {
        if( assay_type == CAGE )
        {
            take_cage_sample(
                rdb, genome,
                i, ss_mi, s_mi,
                MAX_NUM_EM_ITERATIONS,
                max_prb_change_for_convergence
             );
        } else if ( assay_type == CHIP_SEQ ) {
             take_chipseq_sample(
                rdb, genome,
                i, ss_mi, s_mi,
                MAX_NUM_EM_ITERATIONS,
                max_prb_change_for_convergence
             );
        }
    }
    
    /* Create the meta data csv's */
    if( SAVE_STARTING_SAMPLES )
        fclose(ss_mi);

    fclose(s_mi);

    return 0;    
}

void
update_cond_prbs_from_trace_and_assay_type(  
    struct mapped_reads_db* rdb, 
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces,
    struct genome_data* genome,
    enum assay_type_t assay_type
)
{
    struct update_mapped_read_rv_t 
        (*update_reads)( struct cond_prbs_db_t* cond_prbs_db, 
                         const struct trace_t* const traces, 
                         const mapped_read_t* const r  )
        = NULL;

    switch( assay_type )
    {
    case CAGE:
        update_reads = update_CAGE_mapped_read_prbs;
        break;
    
    case CHIP_SEQ:
        update_reads = update_chipseq_mapped_read_prbs;
        break;    

    // case STRANDED_RNASEQ:

    default:
        statmap_log( LOG_FATAL, "Unrecognized assay type in iterative mapping." );
        assert(false);
        exit(-1);
    }
    
    if( NULL == update_reads )
    {
        statmap_log( LOG_ERROR, "Unrecognized assay type (%i). Did you mix statmap versions?",  assay_type );
        exit(1);

    }
    
    /* BUG!!! */
    /* Set the global fl dist */
    global_fl_dist = rdb->fl_dist;
    global_genome = genome;

    update_mapped_reads_from_trace(
        rdb, cond_prbs_db, traces, update_reads
    );
    
    goto cleanup;

cleanup:
    
    return;        
}


/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/
