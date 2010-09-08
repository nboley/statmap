#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <fcntl.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <errno.h>
#include <float.h>
#include <pthread.h>
#include <omp.h>

#include "config.h"
#include "trace.h"
#include "genome.h"

/* get the correct mutex index to access the specified position */
/* this is slow - it should probably be done as a macro in hot areas */
int
get_mutex_index( unsigned int bp)
{
    return bp/TM_GRAN;
}

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size  )
{
    TRACE_TYPE* chr_trace = malloc(sizeof(TRACE_TYPE)*size);
    if( NULL == chr_trace )
    {
        fprintf(stderr, "Can not allocate '%zu' bytes for a trace,\n", sizeof(TRACE_TYPE)*size );
        exit( -1 );
    }
    
    return chr_trace;
}

/* Build mmapped arrays for all of the chrs ***/
void
init_traces( struct genome_data* genome,
             struct trace_t** traces,
             int num_traces )
{
    /* Allocate space for the struct */
    *traces = malloc( sizeof( struct trace_t ) );
    
    /* set the number of traces */
    (*traces)->num_traces = num_traces;
    (*traces)->traces = malloc(sizeof(TRACE_TYPE**)*num_traces);
    assert( (*traces)->traces != NULL );

    /* set the number of chrs */
    (*traces)->num_chrs = genome->num_chrs;
    (*traces)->trace_lengths = malloc(sizeof(unsigned int)*genome->num_chrs);
    assert( (*traces)->trace_lengths != NULL );

    int i;
    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->trace_lengths[i] = genome->chr_lens[i];
    }

    /* allocate space for the locks */
    #ifdef USE_MUTEX
    (*traces)->locks = malloc(sizeof(pthread_mutex_t**)*num_traces);
    #else
    (*traces)->locks = malloc(sizeof(pthread_spinlock_t**)*num_traces);
    #endif
    
    assert( (*traces)->locks != NULL );
    
    /* Allocate space for the pointers to the chr's individual traces and the locks */
    for( i = 0; i < num_traces; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->traces[i] = malloc( (*traces)->num_chrs*sizeof(TRACE_TYPE*) );
        assert( (*traces)->traces[i] != NULL );

        /* Store the trace mutex pointers */
        #ifdef USE_MUTEX
        (*traces)->locks[i] = malloc( (*traces)->num_chrs*sizeof(pthread_mutex_t*) );
        #else
        (*traces)->locks[i] = malloc( (*traces)->num_chrs*sizeof(pthread_spinlock_t*) );
        #endif        
        assert( (*traces)->locks[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {            
            /* initialize the forward trace */
            (*traces)->traces[i][j] = init_trace( (*traces)->trace_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memset( (*traces)->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(*traces)->trace_lengths[j] );
            
            /* initialize the locks */
            int locks_len = (*traces)->trace_lengths[j]/TM_GRAN;
            if( (*traces)->trace_lengths[j] % TM_GRAN > 0 )
                locks_len += 1;
            
            #ifdef USE_MUTEX
            (*traces)->locks[i][j] = malloc( locks_len*sizeof(pthread_mutex_t) );
            #else
            (*traces)->locks[i][j] = malloc( locks_len*sizeof(pthread_spinlock_t) );
            #endif        
            
            int k;
            for( k = 0; k < locks_len; k++ )
            {
                #ifdef USE_MUTEX
                int error = pthread_mutex_init( (*traces)->locks[i][j] + k, 0 );
                #else
                int error = pthread_spin_init( (*traces)->locks[i][j] + k, 0 );
                #endif        
                
                
                if( error != 0 )
                {
                    perror( "Failed to initialize lock in init_trace" );
                    exit( -1 );
                }
            }
        }
    }
    
    return;
}

void
copy_trace_data( struct trace_t* traces,
                 struct trace_t* original )
{
    assert( original->num_traces == traces->num_traces );
    int i;
    for( i = 0; i < original->num_traces; i++ )
    {
        /* Store the pointers for the chrs */
        assert( original->num_chrs == traces->num_chrs );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < traces->num_chrs; j++ )
        {            
            assert( traces->trace_lengths[j] == original->trace_lengths[j] );
            
            memcpy( traces->traces[i][j], original->traces[i][j],
                    sizeof(TRACE_TYPE)*(traces->trace_lengths[j])   );
        }
    }
    
    return;
}

/* init a trace that has the same structure as original, but is inited to 0 */
void
copy_trace( struct trace_t** traces,
            struct trace_t* original )
{
    /* Allocate space for the struct */
    *traces = malloc( sizeof( struct trace_t ) );
    
    /* set the number of traces */
    (*traces)->num_traces = original->num_traces;
    (*traces)->traces = malloc(sizeof(TRACE_TYPE**)*original->num_traces);
    assert( (*traces)->traces != NULL );

    (*traces)->num_chrs = original->num_chrs;
    (*traces)->trace_lengths = malloc(sizeof(unsigned int)*original->num_chrs);
    assert( (*traces)->trace_lengths != NULL );

    int i;
    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->trace_lengths[i] = original->trace_lengths[i];
    }

    /* allocate space for the locks */
    #ifdef USE_MUTEX
    (*traces)->locks = malloc(sizeof(pthread_mutex_t**)*original->num_traces);
    #else
    (*traces)->locks = malloc(sizeof(pthread_spinlock_t**)*original->num_traces);
    #endif

    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < original->num_traces; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->traces[i] = malloc( (*traces)->num_chrs*sizeof(TRACE_TYPE**) );
        assert( (*traces)->traces[i] != NULL );

        /* Store the trace mutex pointers */
        #ifdef USE_MUTEX
        (*traces)->locks[i] = malloc( (*traces)->num_chrs*sizeof(pthread_mutex_t*) );
        #else
        (*traces)->locks[i] = malloc( (*traces)->num_chrs*sizeof(pthread_spinlock_t*) );
        #endif        
        assert( (*traces)->locks[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {            
            (*traces)->traces[i][j] = init_trace( (*traces)->trace_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memcpy( (*traces)->traces[i][j], original->traces[i][j],
                    sizeof(TRACE_TYPE)*(*traces)->trace_lengths[j]   );

            /* initialize the locks */
            int locks_len = (*traces)->trace_lengths[j]/TM_GRAN;
            if( (*traces)->trace_lengths[j] % TM_GRAN > 0 )
                locks_len += 1;
            
            #ifdef USE_MUTEX
            (*traces)->locks[i][j] = malloc( locks_len*sizeof(pthread_mutex_t) );
            #else
            (*traces)->locks[i][j] = malloc( locks_len*sizeof(pthread_spinlock_t) );
            #endif        
            
            int k;
            for( k = 0; k < locks_len; k++ )
            {
                #ifdef USE_MUTEX
                int error = pthread_mutex_init( (*traces)->locks[i][j] + k, 0 );
                #else
                int error = pthread_spin_init( (*traces)->locks[i][j] + k, 0 );
                #endif        
                
                
                if( error != 0 )
                {
                    perror( "Failed to initialize lock in init_trace" );
                    exit( -1 );
                }
            }
        }
    }
    
    return;
}

void
copy_trace_structure( struct trace_t** traces,
                      struct trace_t* original )
{
    copy_trace( traces, original );
    zero_traces( *traces );
    return;
};

void
close_traces( struct trace_t* traces )
{
    int i, j; 
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            free( traces->traces[i][j] );                          

            /* initialize the locks */
            int locks_len = traces->trace_lengths[j]/TM_GRAN;
            if( traces->trace_lengths[j] % TM_GRAN > 0 )
                locks_len += 1;
            
            int k;
            for( k = 0; k < locks_len; k++ )
            {
                #ifdef USE_MUTEX
                int error = pthread_mutex_destroy( traces->locks[i][j] + k );
                #else
                int error = pthread_spin_destroy( traces->locks[i][j] + k );
                #endif        
                              
                if( error != 0 )
                {
                    perror( "Failed to destroy lock in close_trace" );
                    exit( -1 );
                }
            }

            free( traces->locks[i][j] );
        }
        free( traces->traces[i] );                          
        free( traces->locks[i] );
    }

    free( traces->trace_lengths );
    free( traces->traces );
    free( traces->locks );
    free( traces );

    return;
}

void
divide_trace_by_sum( struct trace_t* traces, double value )
{
    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->trace_lengths[j]; k++ )
            {
                traces->traces[i][j][k] /= value;
            }
        }
    }
    return;
}

void
normalize_traces( struct trace_t* traces )
{
    double sum = sum_traces( traces );
    assert( sum > 0 );
    divide_trace_by_sum( traces, sum );
    return;
}

double
sum_traces( struct trace_t* traces )
{
    double sum = 0;

    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->trace_lengths[j]; k++ )
            {
                sum += traces->traces[i][j][k];
            }
        }
    }
    return sum;
}

void
zero_traces( struct trace_t* traces )
{
    /* Zero out the trace for the update */
    int i, j;
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            memset( traces->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(traces->trace_lengths[j]) );
        }
    }
}

void
set_trace_to_uniform( struct trace_t* traces, double value )
{
    /* Zero out the trace for the update */
    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->trace_lengths[j]; k++ )
            {
                traces->traces[i][j][k] = value;
            }
        }
    }
}

void
apply_to_trace( struct trace_t* traces, double (*fun)(double) )
{
    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->trace_lengths[j]; k++ )
            {
                traces->traces[i][j][k] = fun( traces->traces[i][j][k] );
            }
        }
    }
    
    return;
}

/* traces must be the same dimension */
/* This function applies an aggregate to every basepair in the traces */
void
aggregate_over_traces(  struct trace_t* update_trace, 
                        const struct trace_t* const other_trace,
                        TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
                     )
{
    int trace_index, chr;
    unsigned int bp;
    assert( update_trace->num_traces = other_trace->num_traces );
    assert( update_trace->num_chrs = other_trace->num_chrs );
    for( trace_index = 0; trace_index < update_trace->num_traces; trace_index++ )
    {
        for( chr = 0; chr < update_trace->num_chrs; chr++ )
        {
            assert( update_trace->trace_lengths[chr] == other_trace->trace_lengths[chr] );
            for( bp = 0; bp < update_trace->trace_lengths[chr]; bp++ )
            {
                /* update the trace at the correct basepair */
                update_trace->traces[trace_index][chr][bp] 
                    = aggregate( 
                        update_trace->traces[trace_index][chr][bp],   
                        other_trace->traces[trace_index][chr][bp]
                );
            }
        }
    }

    return;
}

