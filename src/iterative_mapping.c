/* Copyright (c) 2009-2010, Nathan Boley */

#define NDEBUG

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
#include <omp.h>
#include <errno.h>
#include <float.h>

#include "statmap.h"
#include "iterative_mapping.h"
#include "mapped_read.h"
#include "trace.h"
#include "wiggle.h"
#include "mapped_location.h"
#include "snp.h"
#include "genome.h"
#include "fragment_length.h"
#include "sam.h"

#define USE_VEC_OPERATIONS
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

static inline void
update_stranded_read_start_density_from_location( 
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc )
{
    const int chr_index = get_chr_from_mapped_read_location( loc  );
    const unsigned char flag = get_flag_from_mapped_read_location( loc  );
    const unsigned int start = get_start_from_mapped_read_location( loc  );
    const float cond_prob = get_cond_prob_from_mapped_read_location( loc  );
    
    if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
    {
        /*
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_lock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #endif
        */
        
        #pragma omp atomic
        traces->traces[1][chr_index][start] += cond_prob;
        // traces->traces[1][chr_index][start] += cond_prob;

        /*
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_unlock( traces->locks[1][ chr_index ] + start/TM_GRAN );
        #endif
        */
    } else {
        /*
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_lock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #endif
        */

        #pragma omp atomic
        traces->traces[0][chr_index][start] += cond_prob;        
        // traces->traces[0][chr_index][start] += cond_prob;

        /*
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #else
        pthread_spin_unlock( traces->locks[0][ chr_index ] + start/TM_GRAN );
        #endif
        */
    }
    
    return;
}


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
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc)
)
{        
    global_fl_dist = reads_db->fl_dist;
    
    /* First, zero the trace to be summed into */
    zero_traces( traces );
    
    if( RAND_MAX < reads_db->num_mmapped_reads )
    {
        fprintf( stderr, "ERROR     : cannot bootstrap random reads - number of reads exceeds RAND_MAX. ( This needs to be fixed. PLEASE file a bug report about this ) \n" );
        return;
    }

    int num_error_reads = 0;
    
    /* Update the trace from reads chosen randomly, with replacement */
    unsigned int i;
    for( i = 0; i < reads_db->num_mmapped_reads; i++ )
    {
        unsigned int read_index = rand()%(reads_db->num_mmapped_reads);
        assert( read_index < reads_db->num_mmapped_reads );
        char* read_start = reads_db->mmapped_reads_starts[ read_index ];
                
        /* read a mapping into the struct */
        struct mapped_read_t r;
        r.read_id = *((unsigned long*) read_start);

        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);

        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (struct mapped_read_location*) read_start;
        
        r.free_locations = false;

        if( i > 0 && i%1000000 == 0 )
        {
            fprintf( stderr, "NOTICE      :  Read %i reads in bootstrap\n", i );
        }

        /* If there are no mappings, then this is pointless */
        if( r.num_mappings > 0 )
        {
            /* Choose a random location, proportional to the normalized probabilities */
            unsigned int j = 0;
            /* if there is exactly 1 mapping, then the randomness is pointless */
            if( r.num_mappings > 1 )
            {
                float random_num = (float)rand()/(float)RAND_MAX;
                double cum_dist = 0;
                for( j = 0; j < r.num_mappings; j++ ) 
                {
                    cum_dist += 
                        get_cond_prob_from_mapped_read_location( r.locations + j );
                    if( random_num <= cum_dist )
                        break;
                }

                /* 
                 * rounding error makes it possible for j == r.num_mappings. If this is 
                 * the case, then make sure the differnece isn't too big, and then chnage 
                 * it. We could save a few cycles by taking this into account earlier, but
                 * for now we're going to leave the sanity check.
                 */
                assert( j <= r.num_mappings );
                if( j == r.num_mappings )
                {
                    j = r.num_mappings - 1;
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
            /* We do this by temp changing cond_prb, and then changing it back */
            float tmp_cond_prb = 
                get_cond_prob_from_mapped_read_location( r.locations + j );
            set_cond_prob_in_mapped_read_location( 
                r.locations + j, 1.0  );
            update_trace_expectation_from_location( traces, r.locations + j );
            set_cond_prob_in_mapped_read_location( 
                r.locations + j, tmp_cond_prb  );
        }
    }
    
    if( num_error_reads > 0 )
    {
        fprintf( stderr, "ERROR      :  %i reads ( out of %lu ) had errors in which the cum dist didnt add up to 1. \n", 
                 num_error_reads, reads_db->num_mmapped_reads );
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
    struct trace_t* traces;
    void (* update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc);
};

void*
update_traces_from_mapped_reads_worker( void* params )
{
    struct trace_t* traces = ( (struct update_traces_param*) params)->traces;
    struct mapped_reads_db* rdb = 
        ( (struct update_traces_param*) params)->reads_db;
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc)
        = ( (struct update_traces_param*) 
            params)->update_trace_expectation_from_location;
    
    struct mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) )     
    {        
        /* Update the trace from this mapping */        
        unsigned int j;
        for( j = 0; j < r->num_mappings; j++ ) {
            // update_stranded_read_start_density_from_location( traces, r->locations + j );
            update_trace_expectation_from_location( traces, r->locations + j );
        }
        
        free_mapped_read( r );
    }
    
    free_mapped_read( r );
    
    pthread_exit( NULL );
    
    return 0;
}

void
update_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc)
)
{        
    zero_traces( traces );

    //struct trace_t* locs_trace;
    //copy_trace_structure( &locs_trace, traces );

    /* Move the database ( this should be a cursor ) back to the first read */
    rewind_mapped_reads_db( reads_db );
    
    /* If the number of threads is one, then just do everything in serial */
    if( num_threads == 1 )
    {
        struct mapped_read_t* r;
        int read_num = 0;
        while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) )     
        {
            read_num++;
            if( read_num%1000000 == 0 )
                fprintf( stderr, "DEBUG       :  Updated traces from mapped reads for %i reads\n", read_num );
            
            /* Update the trace from this mapping */        
            unsigned int j;
            for( j = 0; j < r->num_mappings; j++ ) {
                // update_stranded_read_start_density_from_location( traces, r->locations + j );
                update_trace_expectation_from_location( traces, r->locations + j );
            }
            
            free_mapped_read( r );
        }

        free_mapped_read( r );
        
    } 
    /* otherwise, if we are expecting more than one thread */
    else {
        /* initialize the thread parameters structure */
        struct update_traces_param params;
        
        params.reads_db = reads_db;
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
                fprintf(stderr, 
                        "ERROR; return code from pthread_create() is %d\n", 
                        rc
                    );
                exit(-1);
            }
        }

        /* wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            if (rc) {
                fprintf( stderr, 
                         "ERROR; return code from pthread_join() is %d\n", 
                         rc
                    );
                exit(-1);
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
    unsigned int* curr_read_index;
    
    struct mapped_reads_db* reads_db;
    struct trace_t* traces;
    struct update_mapped_read_rv_t 
        (* update_mapped_read_prbs)( const struct trace_t* const traces, 
                                     const struct mapped_read_t* const r  );
    
};

void*
update_mapped_reads_from_trace_worker( void* params )
{
    struct trace_t* traces = ( (struct update_mapped_reads_param*) params)->traces;
    struct mapped_reads_db* reads_db = 
        ( (struct update_mapped_reads_param*) params)->reads_db;

    struct update_mapped_read_rv_t 
        (* update_mapped_read_prbs)( const struct trace_t* const traces, 
                                     const struct mapped_read_t* const r  )
        = ( (struct update_mapped_reads_param*)
            params)->update_mapped_read_prbs;


    struct mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) ) 
    {
        /* Update the read */
        struct update_mapped_read_rv_t tmp_rv 
            = update_mapped_read_prbs( traces, r );

        /* Hand reduce the max */
        if( tmp_rv.max_change > ( (struct update_mapped_reads_param*) params)->rv.max_change )
            ( (struct update_mapped_reads_param*) params)->rv.max_change = tmp_rv.max_change;
        
        /* Update the lhd */
        ( (struct update_mapped_reads_param*) params)->rv.log_lhd += tmp_rv.log_lhd;

        free_mapped_read( r );
    }
    
    free_mapped_read( r );

    pthread_exit( NULL );
    
    return 0;
}


struct update_mapped_read_rv_t
update_mapped_reads_from_trace( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
    )
{    
    /* store the toal accumulated error */
    struct update_mapped_read_rv_t rv = { 0, 0.0 };
    
    /* reset the cursor */
    rewind_mapped_reads_db( reads_db );

    int read_cnt = 0;

    /* If the number of threads is one, then just do everything in serial */
    if( num_threads == 1 )
    {
        struct mapped_read_t* r;
        
        while( EOF != get_next_read_from_mapped_reads_db( reads_db, &r ) ) 
        {
            read_cnt++;
            if( read_cnt%1000000 == 0 )
                fprintf( stderr, "DEBUG       :  Updated reads from traces for %i reads\n", read_cnt );
            
            /* Update the read */
            struct update_mapped_read_rv_t tmp_rv 
                = update_mapped_read_prbs( traces, r );
            
            /* Hand reduce the max */
            if( tmp_rv.max_change > rv.max_change )
            {
                rv.max_change = tmp_rv.max_change;
            }
            
            /* Update the lhd */
            rv.log_lhd += tmp_rv.log_lhd;

            free_mapped_read( r );
        }

        free_mapped_read( r );
        
    } else {
        /* initialize the thread parameters structure */
        struct update_mapped_reads_param params;
        
        params.rv.max_change = 0;
        params.rv.log_lhd = 0;
                
        params.reads_db = reads_db;
        /* traces should be read only, so they are shared */
        params.traces = traces;
        params.update_mapped_read_prbs = update_mapped_read_prbs;
        
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
                                 update_mapped_reads_from_trace_worker, 
                                 (void *)(tds + t) 
                ); 
            if (rc) {
                fprintf(stderr, 
                        "ERROR; return code from pthread_create() is %d\n", 
                        rc
                    );
                exit(-1);
            }
        }

        /* wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            pthread_attr_destroy( attrs + t );
            
            if (rc) {
                fprintf( stderr, 
                         "ERROR; return code from pthread_join() is %d\n", 
                         rc
                    );
                exit(-1);
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
    struct trace_t* traces,
    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
    )
{    
    struct update_mapped_read_rv_t rv;
    
    rv = update_mapped_reads_from_trace( 
        reads_db, traces, update_mapped_read_prbs );
    
    return rv.log_lhd;
}

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
    )
{
    /* make sure the mapped reads are mmapped */
    assert( rdb->write_locked );

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
            rdb, starting_trace, 
            update_mapped_read_prbs
        );
        gettimeofday( &umr_stop, NULL );   
        
        gettimeofday( &ut_start, NULL );   
        update_traces_from_mapped_reads( 
            rdb, starting_trace, 
            update_trace_expectation_from_location
        );
        gettimeofday( &ut_stop, NULL );

        #define LHD_RATIO_STOP_VAL 1.05
        
        if( num_iterations > 0 &&
            ( 
                ( 
                  num_iterations == 1 
                  || num_iterations%25 == 0 
                  || (ut_stop.tv_sec-nt_start.tv_sec) > 30 
                )
                || rv.max_change < max_prb_change_for_convergence
                || ( pow( 10, rv.log_lhd - prev_log_lhd ) < LHD_RATIO_STOP_VAL 
                     && pow( 10, rv.log_lhd - prev_log_lhd ) > 0.95
                   )
                )
            )
        {
            fprintf( stderr, "Iter %i: \tError: %e \tLog Lhd: %e (ratio %e) \tNorm Trace:  %.5f sec\t Read UT:  %.5f sec\tTrace UT:  %.5f sec\n", 
                 num_iterations, rv.max_change, rv.log_lhd,
                     pow( 10, rv.log_lhd - prev_log_lhd ),
                     (float)(nt_stop.tv_sec - nt_start.tv_sec) + ((float)(nt_stop.tv_usec - nt_start.tv_usec))/1000000,
                     (float)(umr_stop.tv_sec - umr_start.tv_sec) + ((float)(umr_stop.tv_usec - umr_start.tv_usec))/1000000,
                     (float)(ut_stop.tv_sec - ut_start.tv_sec) + ((float)(ut_stop.tv_usec - ut_start.tv_usec))/1000000
            );
            
            if( num_iterations > 1 
                && ( rv.max_change < max_prb_change_for_convergence 
                || ( pow( 10, rv.log_lhd - prev_log_lhd ) < LHD_RATIO_STOP_VAL 
                     && pow( 10, rv.log_lhd - prev_log_lhd ) > 0.95
                   )
                )
              )
                break;
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
    
    void (* const update_trace_expectation_from_location)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc),

    struct update_mapped_read_rv_t 
        (* const update_mapped_read_prbs)( const struct trace_t* const traces, 
                                           const struct mapped_read_t* const r  )
    )
{
    global_fl_dist = rdb->fl_dist;
    global_genome = genome;
    
    /* int traces to a low consant. This should heavily favor prior reads. */
    set_trace_to_uniform( traces, EXPLORATION_PRIOR );

    rewind_mapped_reads_db( rdb );
    struct mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) ) 
    {        
        /* reset the read cond prbs under a uniform prior */
        reset_read_cond_probs( r );
        
        /* update the read conditional probabilities from the trace */
        update_mapped_read_prbs( traces, r );        
        
        if( r->num_mappings > 0 )
        {
            /* Choose a random loc, proportional to the normalized probabilities */
            unsigned int j = 0;
            /* if there is exactly 1 mapping, then the randomness is pointless */
            if( r->num_mappings > 1 )
            {
                float random_num = (float)rand()/(float)RAND_MAX;
                double cum_dist = 0;
                for( j = 0; j < r->num_mappings; j++ ) 
                {
                    cum_dist += 
                        get_cond_prob_from_mapped_read_location( r->locations + j );
                    if( random_num <= cum_dist )
                        break;
                }

                /* 
                 * rounding error makes it possible for j == r.num_mappings. If this is 
                 * the case, then make sure the differnece isn't too big, and then chnage 
                 * it. We could save a few cycles by taking this into account earlier, but
                 * for now we're going to leave the sanity check.
                 */
                assert( j <= r->num_mappings );
                if( j == r->num_mappings )
                {
                    // printf( "WARNING - rounding error! %e\n", cum_dist );
                    j = r->num_mappings - 1;
                    assert( cum_dist > 0.999 );
                }
            }
                        
            /* Add this location to the trace, with prb 1 */
            /* We do this by temp changing cond_prb, and then changing it back */
            float tmp_cond_prb = 
                get_cond_prob_from_mapped_read_location( r->locations + j );
            set_cond_prob_in_mapped_read_location( r->locations + j, 1.0 );
            update_trace_expectation_from_location( traces, r->locations + j );
            set_cond_prob_in_mapped_read_location( r->locations + j, tmp_cond_prb );
        }

        free_mapped_read( r );
    }
    
    free_mapped_read( r );
    
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
                          
)
{
    /* Store meta information about the samples and starting samples */
    FILE* ss_mi;
    FILE* s_mi;
    
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
        printf( "Starting Sample %i\n", i+1 );

        struct trace_t* sample_trace;
        init_trace( genome, &sample_trace, num_tracks, track_names );

        build_random_starting_trace( 
            sample_trace, genome, rdb, 
            update_trace_expectation_from_location,
            update_mapped_read_prbs
        );
        
        if( SAVE_STARTING_SAMPLES )
        {
            char buffer[100];
            sprintf( buffer, "%ssample%i.wig", STARTING_SAMPLES_PATH, i+1 );

            write_wiggle_from_trace(
                sample_trace, 
                buffer, max_prb_change_for_convergence );
            
            double log_lhd = calc_log_lhd( rdb, sample_trace, update_mapped_read_prbs );
            fprintf( ss_mi, "%i,%e\n", i+1, log_lhd );
            fflush( ss_mi );
        }
        
        /* maximize the likelihood */
        update_mapping( 
            rdb, sample_trace, max_num_iterations,
            max_prb_change_for_convergence,
            update_trace_expectation_from_location,
            update_mapped_read_prbs
        );
        
        if( SAVE_SAMPLES )
        {
            char buffer[100];
            sprintf( buffer, "%ssample%i.wig", RELAXED_SAMPLES_PATH, i+1 );

            write_wiggle_from_trace( 
                sample_trace, 
                buffer, max_prb_change_for_convergence );
            
            double log_lhd = calc_log_lhd( rdb, sample_trace, update_mapped_read_prbs );
            fprintf( s_mi, "%i,%e\n", i+1, log_lhd );
            fflush( s_mi );
        }
        
        if( NUM_BOOTSTRAP_SAMPLES > 0 )
        {
            fprintf( stderr, "Bootstrapping %i samples.", NUM_BOOTSTRAP_SAMPLES );
            char buffer[200];
            sprintf( buffer, "mkdir %ssample%i/",
                     BOOTSTRAP_SAMPLES_ALL_PATH, i+1 );
            int error = system( buffer );
            if (WIFSIGNALED(error) &&
                (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
            {
                fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
                perror( "System Call Failure");
                assert( false );
                exit( -1 );
            }
            
            int j;
            for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
            {                
                
                if(  j%(NUM_BOOTSTRAP_SAMPLES/10)  == 0 )
                    fprintf( stderr, " %.1f%%...", (100.0*j)/NUM_BOOTSTRAP_SAMPLES );

                /* actually perform the bootstrap */
                bootstrap_traces_from_mapped_reads( 
                    rdb, sample_trace, update_trace_expectation_from_location );
                
                if( SAVE_BOOTSTRAP_SAMPLES )
                {
                    sprintf( buffer, "%ssample%i/bssample%i.wig", 
                             BOOTSTRAP_SAMPLES_ALL_PATH, i+1, j+1 );
                    
                    write_wiggle_from_trace( 
                        sample_trace, 
                        buffer, max_prb_change_for_convergence );
                }    
            }
            fprintf( stderr, " 100%%\n");
        }
        
        if( SAVE_AGGREGATED_BOOTSTRAP_SAMPLES )
        {
            char buffer[200];
            FILE** fps = malloc(sizeof(FILE*)*NUM_BOOTSTRAP_SAMPLES);
            
            int j;
            for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
            {
                sprintf( buffer, "%ssample%i/bssample%i.wig", 
                         BOOTSTRAP_SAMPLES_ALL_PATH, i+1, j+1 );
                fps[j] = fopen( buffer, "r" );
            }

            /* make the MIN aggregates */
            sprintf( buffer, "%ssample%i.wig", 
                     BOOTSTRAP_SAMPLES_MIN_PATH, i+1 );
            FILE* ofp = fopen( buffer, "w" );
            aggregate_over_wiggles( 
                fps, NUM_BOOTSTRAP_SAMPLES, ofp, FLT_EPSILON, wig_lines_min  );
            fclose( ofp );
            for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
                fseek( fps[j], 0, SEEK_SET );
            
            /* make the MAX aggregates */
            sprintf( buffer, "%ssample%i.wig", 
                     BOOTSTRAP_SAMPLES_MAX_PATH, i+1 );
            ofp = fopen( buffer, "w" );
            aggregate_over_wiggles( 
                fps, NUM_BOOTSTRAP_SAMPLES, ofp, FLT_EPSILON, wig_lines_max  );
            fclose( ofp );
            for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
                fclose( fps[j] );
            
            free( fps );
        }
        
        close_traces( sample_trace );
    }

    if( SAVE_AGGREGATED_BOOTSTRAP_SAMPLES )
    {
        char buffer[200];
        FILE** fps = malloc(sizeof(FILE*)*num_samples);
        
        int j;

        /* make the MIN aggregates */
        for( j = 0; j < num_samples; j++ )
        {
            sprintf( buffer, "%ssample%i.wig", BOOTSTRAP_SAMPLES_MIN_PATH, j+1 );
            fps[j] = fopen( buffer, "r" );
            assert( fps[j] != NULL );
        }
        FILE* ofp = fopen( MIN_TRACE_FNAME, "w" );
        aggregate_over_wiggles( fps, num_samples, ofp, FLT_EPSILON, wig_lines_min  );
        fclose( ofp );
        for( j = 0; j < num_samples; j++ )
            fclose( fps[j] );
        
        /* make the MAX aggregates */
        for( j = 0; j < num_samples; j++ )
        {
            sprintf( buffer, "%ssample%i.wig", BOOTSTRAP_SAMPLES_MAX_PATH, j+1 );
            fps[j] = fopen( buffer, "r" );
            assert( fps[j] != NULL );
        }
        ofp = fopen( MAX_TRACE_FNAME, "w" );
        aggregate_over_wiggles( fps, num_samples, ofp, FLT_EPSILON, wig_lines_max  );
        fclose( ofp );
        for( j = 0; j < num_samples; j++ )
            fclose( fps[j] );
        
        free( fps );
    }

    /* Create the meta data csv's */
    if( SAVE_STARTING_SAMPLES )
        fclose(ss_mi);

    fclose(s_mi);

    
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
    const struct mapped_read_location* const loc )
{
    const unsigned char flag = get_flag_from_mapped_read_location( loc  );
    const int chr_index = get_chr_from_mapped_read_location( loc  );
    const unsigned int start = get_start_from_mapped_read_location( loc  );
    const unsigned int stop = get_stop_from_mapped_read_location( loc  );
    const float cond_prob = get_cond_prob_from_mapped_read_location( loc  );
    
    assert( cond_prob >= 0.0 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_chrs );
    assert( traces->chr_lengths[chr_index] >= stop );
    
    /* iteration variable */
    unsigned int k;

    /* update the trace */
    /* If the reads are paired */
    if( (flag&IS_PAIRED) != 0 )
    {
        const float scale = (1.0/(stop-start));

        int trace_index;
        if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
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
        assert( 0 == (flag&IS_PAIRED) ); 
        /* the binding site must be on the five prime side of the read */
        unsigned int window_start;
        unsigned int window_stop;
        /* determine which trace to lock */
        int trace_index;
        if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
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
            if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
            {
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
update_chipseq_mapped_read_prbs( const struct trace_t* const traces, 
                                 const struct mapped_read_t* const r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };
    
    /* allocate space to store the temporary values */
    float* new_prbs = malloc( sizeof(float)*(r->num_mappings) );
    
    /* Update the reads from the trace */
    /* set to flt_min to avoid division by 0 */
    double prb_sum = ML_PRB_MIN;
    double error_prb_sum = ML_PRB_MIN;

    unsigned int i;
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* calculate the mean density */
        double window_density = 0;

        int chr_index = get_chr_from_mapped_read_location( r->locations + i  );
        unsigned char flag = get_flag_from_mapped_read_location( r->locations + i );
        unsigned int start = get_start_from_mapped_read_location( r->locations + i );
        unsigned int stop = get_stop_from_mapped_read_location( r->locations + i );

        /* If the reads are paired */
        unsigned int j = 0;
        if( (flag&IS_PAIRED) > 0 )
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
                get_fl_from_mapped_read_location( r->locations + i ) 
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
            if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
            {
                window_start = MAX( stop-global_fl_dist->max_fl, 0 );
                window_stop = MIN( stop, global_genome->chr_lens[chr_index] );
            } else {
                window_start = start;
                window_stop = MIN( start + global_fl_dist->max_fl, 
                                   global_genome->chr_lens[chr_index] );
            }

            if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
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
            get_seq_error_from_mapped_read_location( r->locations + i )
               *window_density;
        prb_sum += new_prbs[i];        
        error_prb_sum += get_seq_error_from_mapped_read_location( r->locations + i );
    }    
    
    /* renormalize the read probabilities */
    if( prb_sum > 0 )
    {
        for( i = 0; i < r->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* absolute value error */
            double normalized_density = new_prbs[i]/prb_sum; 
            
            rv.max_change += MAX( 
                normalized_density 
                    - get_cond_prob_from_mapped_read_location( r->locations + i ),
                get_cond_prob_from_mapped_read_location( r->locations + i ) 
                    - normalized_density 
            );
            
            set_cond_prob_in_mapped_read_location( 
                r->locations + i, normalized_density );
        }
        /* BUG - Can we be smarter here? */
        /* if the window density sums to 0, then it is impossible to normalize. 
           so we assume uniformity over the reads. I dont know if this is the right
           thing to do ( in fact, im worried that it may destroy convergence. However,
           it's the only thing to think of except ignoring the read, so I do. Note 
           that this should not happen with any trace that is generated by the same 
           stream. ie, this is really only a problem for the starting location trace.
        */
    } else {
        for( i = 0; i < r->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* absolute value error */
            double normalized_density = 
                get_seq_error_from_mapped_read_location( r->locations + i )
                /error_prb_sum; 
            set_cond_prob_in_mapped_read_location( 
                r->locations + i, normalized_density );
        }
    }
    
    rv.log_lhd = log10( prb_sum );
    
    /* Free the tmp array */
    free( new_prbs );

    return rv;
}

/*****************************************************************************
 * 
 * CAGE specific functions 
 *
 *****************************************************************************/

void update_CAGE_trace_expectation_from_location(
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc )
{
    int chr_index = get_chr_from_mapped_read_location( loc  );
    unsigned char flag = get_flag_from_mapped_read_location( loc  );
    unsigned int start = get_start_from_mapped_read_location( loc  );
    // unsigned int stop = get_stop_from_mapped_read_location( loc  );
    float cond_prob = get_cond_prob_from_mapped_read_location( loc  );
    
    assert( cond_prob >= 0.0 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_chrs );
    assert( traces->chr_lengths[chr_index] >= stop );

    assert( cond_prob >= 0.0 );
    
    /* update the trace */
    /* If the reads are paired */
    if( (flag&IS_PAIRED) != 0 )
    {
        fprintf( stderr, "FATAL: paired cage reads are not supported" );
        exit( -1 );                
    } 
    /* If the read is *not* paired */
    else {
        int trace_index;

        /* If this is in the fwd strnanded transcriptome */
        if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
        {
            trace_index = 0;
        } 
        /* We are in the 3' ( negative ) transcriptome */
        else {
            trace_index = 1;
        }

        /* lock the spinlock */
        #ifdef USE_MUTEX
        pthread_mutex_lock( traces->locks[ trace_index ][ chr_index ] + (start/TM_GRAN) );
        #else
        pthread_spin_lock( traces->locks[ trace_index ][ chr_index ] + (start/TM_GRAN) );
        #endif        
        traces->traces[ trace_index ][ chr_index ][ start ] += cond_prob; 
        #ifdef USE_MUTEX
        pthread_mutex_unlock( traces->locks[ trace_index ][ chr_index ] + (start/TM_GRAN) );
        #else
        pthread_spin_unlock( traces->locks[ trace_index ][ chr_index ] + (start/TM_GRAN) );
        #endif
    }
    
    return;
}

/* Returns the probability of observing the read, conditional on the trace */
inline struct update_mapped_read_rv_t 
update_CAGE_mapped_read_prbs( 
    const struct trace_t* const traces, 
    const struct mapped_read_t* const r  )
{
    struct update_mapped_read_rv_t rv = { 0, 0 };
    
    /* allocate space to store the temporary values */
    float* new_prbs = malloc( sizeof(float)*(r->num_mappings) );
    
    /* Update the reads from the trace */
    double density_sum = 0;
    unsigned int i;
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* calculate the mean density */
        /* We set this to 2*DBL_EPSILON to prevent the division by 0 */
        double window_density = 2*DBL_EPSILON;

        int chr_index = get_chr_from_mapped_read_location( r->locations + i  );
        unsigned char flag = get_flag_from_mapped_read_location( r->locations + i );
        unsigned int start = get_start_from_mapped_read_location( r->locations + i );
        
        /* If the reads are paired */
        unsigned int j = 0;
        if( (flag&IS_PAIRED) > 0 )
        {
            fprintf( stderr, "FATAL: paired cage reads are not supported\n" );
            exit( -1 );
        } 
        /* If the read is *not* paired */
        else {
            /* store the trace that we care about */
            TRACE_TYPE** trace;
            
            /* If this is in the fwd strnanded transcriptome */
            if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
            {
                trace = traces->traces[0];
            } 
            /* We are in the 3' ( negative ) transcriptome */
            else {
                trace = traces->traces[1];
            }
            
            unsigned int stop = MIN( traces->chr_lengths[chr_index], start + WINDOW_SIZE );
            for( j = start; j < stop; j++ )
            {
                window_density += trace[chr_index][j];
            }
        }
        
        new_prbs[i] = 
            get_seq_error_from_mapped_read_location( r->locations + i )
                *window_density;
        
        density_sum += new_prbs[i];
    }
    
    /* renormalize the read probabilities */
    if( r->num_mappings > 0 )
    {
        assert( density_sum > 0 );
        
        for( i = 0; i < r->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* quadratic error */
            // abs_error += pow(new_prbs[i]/density_sum 
            //                 - r.locations[i].cond_prob, 2 ) ;
            /* absolute value error */
            double normalized_density = new_prbs[i]/density_sum; 

            rv.max_change += MAX( 
                normalized_density - new_prbs[i], new_prbs[i] - normalized_density 
            );
            
            set_cond_prob_in_mapped_read_location( 
                r->locations + i, normalized_density );
            
            assert( get_cond_prob_from_mapped_read_location( r->locations + i ) >= 0 );
        }
    }      
    
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
    struct mapped_reads_db* ip_rdb, 
    struct trace_t** ip_trace,
    
    struct mapped_reads_db* nc_rdb,
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
    reset_all_read_cond_probs( ip_rdb );
    reset_all_read_cond_probs( nc_rdb );
            
    /* iteratively map from a uniform prior */
    gettimeofday( &start, NULL );
    fprintf(stderr, "NOTICE      :  Starting iterative mapping.\n" );
    
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
            *ip_trace, genome, ip_rdb, 
            update_chipseq_trace_expectation_from_location,
            update_chipseq_mapped_read_prbs
        );
    }
    
    /* reset the read cond prbs under a uniform prior */
    reset_all_read_cond_probs( ip_rdb );
    
    /* update the 'real' IP data */
    error = update_mapping (
        ip_rdb, 
        *ip_trace,
        MAX_NUM_EM_ITERATIONS,
        max_prb_change_for_convergence,
        update_chipseq_trace_expectation_from_location,
        update_chipseq_mapped_read_prbs
    );
    
    gettimeofday( &stop, NULL );
    long seconds  = stop.tv_sec  - start.tv_sec;
    long useconds = stop.tv_usec - start.tv_usec;
    
    fprintf(stderr, "PERFORMANCE :  Maximized LHD in %.5f seconds\n", ((float)seconds) + ((float)useconds)/1000000 );
    
    /* update the NC data based upon the IP data's NC */
    normalize_traces( *ip_trace );
    
    struct update_mapped_read_rv_t rv;
    rv = update_mapped_reads_from_trace(
            nc_rdb, *ip_trace, 
            update_chipseq_mapped_read_prbs
        );
    
    /* update the NC reads from the ip marginal density */
    update_traces_from_mapped_reads( 
        nc_rdb, *nc_trace, 
        update_chipseq_trace_expectation_from_location
    );

    /* we need to do this because I renormalized the read sensity, 
       and it screws up the wiggle printing. Kind of a hack... */
    update_traces_from_mapped_reads( 
        ip_rdb, *ip_trace, 
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
    
    /* jointly update the mappings */
    update_chipseq_mapping_wnc(  chip_mpd_rds_db, &ip_trace,
                                 NC_mpd_rds_db, &nc_trace,
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
        sprintf( wig_fname, "%ssample%i.ip.wig", RELAXED_SAMPLES_PATH, sample_index+1 );
        write_wiggle_from_trace( ip_trace, wig_fname, MAX_PRB_CHANGE_FOR_CONVERGENCE );
                
        sprintf( wig_fname, "%ssample%i.nc.wig", RELAXED_SAMPLES_PATH, sample_index+1 );
        write_wiggle_from_trace( nc_trace, wig_fname, MAX_PRB_CHANGE_FOR_CONVERGENCE );

        /* write the lhd to the meta info folder */
        double log_lhd = calc_log_lhd( 
            chip_mpd_rds_db, ip_trace, update_chipseq_mapped_read_prbs );
        fprintf( meta_info_fp, "%i,%e\n", sample_index+1, log_lhd );
        fflush( meta_info_fp );
    }
            
    if( NUM_BOOTSTRAP_SAMPLES > 0 )
    {
        fprintf( stderr, "Bootstrapping %i samples.", NUM_BOOTSTRAP_SAMPLES );
        char buffer[200];
        sprintf( buffer, "mkdir %ssample%i/",
                 BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1 );
        int error = system( buffer );
        if (WIFSIGNALED(error) &&
            (WTERMSIG(error) == SIGINT || WTERMSIG(error) == SIGQUIT))
        {
            fprintf(stderr, "FATAL     : Failed to call '%s'\n", buffer );
            perror( "System Call Failure");
            assert( false );
            exit( -1 );
        }

        int j;
        for( j = 0; j < NUM_BOOTSTRAP_SAMPLES; j++ )
        {                
            if(  j%(NUM_BOOTSTRAP_SAMPLES/10)  == 0 )
                fprintf( stderr, " %.1f%%...", (100.0*j)/NUM_BOOTSTRAP_SAMPLES );
                    
            /* actually perform the bootstrap */
            bootstrap_traces_from_mapped_reads( 
                chip_mpd_rds_db, ip_trace, 
                update_chipseq_trace_expectation_from_location );
                    
            /* actually perform the bootstrap */
            bootstrap_traces_from_mapped_reads( 
                NC_mpd_rds_db, nc_trace,
                update_chipseq_trace_expectation_from_location );
                    
            if( SAVE_BOOTSTRAP_SAMPLES )
            {
                /* write out the IP */
                sprintf( buffer, "%ssample%i/bssample%i.ip.wig", 
                         BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1, j+1 );                        

                write_wiggle_from_trace( 
                    ip_trace, buffer, MAX_PRB_CHANGE_FOR_CONVERGENCE );

                /* write out the NC */
                sprintf( buffer, "%ssample%i/bssample%i.nc.wig", 
                         BOOTSTRAP_SAMPLES_ALL_PATH, sample_index+1, j+1 );                        
                        
                write_wiggle_from_trace( 
                    nc_trace, buffer, MAX_PRB_CHANGE_FOR_CONVERGENCE );
            }    
        }
        fprintf( stderr, " 100%%\n");
                
        close_traces( ip_trace );
        close_traces( nc_trace );
    }     

    return;
}


int
generic_update_mapping(  struct mapped_reads_db* rdb, 
                         struct genome_data* genome,
                         enum assay_type_t assay_type,
                         int num_samples,
                         float max_prb_change_for_convergence)
{
    /* tell whether or not we *can* iteratively map ( ie, do we know the assay? ) */
    clock_t start, stop;
    
    int error = 0;
    
    void (*update_expectation)(
        const struct trace_t* const traces, 
        const struct mapped_read_location* const loc) 
        = NULL;
    
    struct update_mapped_read_rv_t 
        (*update_reads)( const struct trace_t* const traces, 
                         const struct mapped_read_t* const r  )
        = NULL;

    char** track_names = NULL;

    int trace_size = -1;

    switch( assay_type )
    {
    case CAGE:
        update_expectation = update_CAGE_trace_expectation_from_location;
        update_reads = update_CAGE_mapped_read_prbs;
        trace_size = 2;
        track_names = malloc( trace_size*sizeof(char*) );
        track_names[0] = "fwd_strnd_read_density"; 
        track_names[1] = "rev_strnd_read_density";
        break;
    
    case CHIP_SEQ:
        update_expectation = update_chipseq_trace_expectation_from_location;
        update_reads = update_chipseq_mapped_read_prbs;
        trace_size = 2;
        track_names = malloc( trace_size*sizeof(char*) );
        track_names[0] = "fwd_strand_read_density";
        track_names[1] = "bkwd_strand_read_density";
        break;
    
    default:
        fprintf( stderr, "WARNING     :  Can not iteratively map for assay type '%u'. Returning marginal mappings.\n", assay_type);
        return 0;
    }
    
    /* BUG!!! */
    /* Set the global fl dist */
    global_fl_dist = rdb->fl_dist;
    global_genome = genome;
    
    /* reset the read cond prbs under a uniform prior */
    reset_all_read_cond_probs( rdb );
    
    struct trace_t* uniform_trace;
    
    /* iteratively map from a uniform prior */
    start = clock();
    fprintf(stderr, "NOTICE      :  Starting iterative mapping.\n" );
    
    /* initialize the trace that we will store the expectation in */
    init_trace( genome, &uniform_trace, trace_size, track_names );
    set_trace_to_uniform( uniform_trace, 1 );
    
    error = update_mapping (
        rdb, 
        uniform_trace,
        MAX_NUM_EM_ITERATIONS,
        max_prb_change_for_convergence,
        update_expectation,
        update_reads
    );
    
    write_wiggle_from_trace( 
        uniform_trace, 
        "relaxed_mapping.wig", max_prb_change_for_convergence );
    
    stop = clock();
    fprintf(stderr, "PERFORMANCE :  Maximized LHD in %.2lf seconds\n", 
            ((float)(stop-start))/CLOCKS_PER_SEC );
    
    /* If users want to do this, why not just call the utility? */
    #if 0 
        /* write the mapped reads to SAM */
        start = clock();
        fprintf(stderr, "NOTICE      :  Writing relaxed mapped reads to SAM file.\n" );
        
        FILE* sam_ofp = fopen( "mapped_reads_relaxed.sam", "w+" );
        write_mapped_reads_to_sam( 
            rawread_db, rdb, genome, false, false, sam_ofp );
        fclose( sam_ofp );    
        
        stop = clock();
        fprintf(stderr, 
                "PERFORMANCE :  Wrote relaxed mapped reads to sam in %.2lf seconds\n", 
                ((float)(stop-start))/CLOCKS_PER_SEC );
    #endif
    
    error = sample_random_traces(
        rdb, genome, 
        trace_size, track_names,
        num_samples, MAX_NUM_EM_ITERATIONS, 
        max_prb_change_for_convergence,
        update_expectation, update_reads
    );
    
    goto cleanup;

cleanup:

    free( track_names );

    close_traces( uniform_trace );
    
    return 0;    
}


/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/
