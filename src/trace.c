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

#include "config.h"
#include "trace.h"
#include "genome.h"

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size  )
{
    void* chr_trace
        = mmap( NULL, sizeof(TRACE_TYPE)*size, 
                PROT_READ|PROT_WRITE, MAP_ANON|MAP_POPULATE|MAP_PRIVATE, 0, 0 );
    
    if( NULL == chr_trace || chr_trace == (void*) -1 )
    {
        fprintf(stderr, "Can not anonymous mmap of size '%zu'\n", size );
        exit( -1 );
    }
    
    return (TRACE_TYPE*) chr_trace;
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

    (*traces)->num_chrs = genome->num_chrs;
    (*traces)->trace_lengths = malloc(sizeof(unsigned int)*genome->num_chrs);
    assert( (*traces)->trace_lengths != NULL );

    int i;
    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->trace_lengths[i] = genome->chr_lens[i];
    }


    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < num_traces; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->traces[i] = malloc( (*traces)->num_chrs*sizeof(TRACE_TYPE**) );
        assert( (*traces)->traces[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {            
            /* initialize the forward trace */
            (*traces)->traces[i][j] = init_trace( (*traces)->trace_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memset( (*traces)->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(*traces)->trace_lengths[j] );
        }
    }
    
    return;
}

/* init a trace that has the same structure as original, but is inited to 0 */
void
copy_trace_structure( struct trace_t** traces,
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


    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < original->num_traces; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->traces[i] = malloc( (*traces)->num_chrs*sizeof(TRACE_TYPE**) );
        assert( (*traces)->traces[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {            
            /* initialize the forward trace */
            (*traces)->traces[i][j] = init_trace( (*traces)->trace_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memset( (*traces)->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(*traces)->trace_lengths[j] );
        }
    }
    
    return;
}


void
close_traces( struct trace_t* traces )
{
    int i, j, error; 
    for( i = 0; i < traces->num_traces; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            error = munmap( traces->traces[i][j], 
                            sizeof(TRACE_TYPE)*traces->trace_lengths[j] );
            if( error == -1 )
            {
                perror( "Problem closing mmaped trace" );
                assert( false );
                exit( -1 );
            }
        }
    }

    free( traces->trace_lengths );
    free( traces );

    return;
}

double
sum_traces( struct trace_t* traces )
{
    double sum = 0;

    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_traces; i++ )
        for( j = 0; j < traces->num_chrs; j++ )
            for( k = 0; k < traces->trace_lengths[i]; k++ )
                sum += traces->traces[i][j][k];
    
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
                    sizeof(TRACE_TYPE)*(traces->trace_lengths[i]) );
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
            for( k = 0; k < traces->trace_lengths[i]; k++ )
            {
                traces->traces[i][j][k] = value;
            }
        }
    }
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


extern void
write_wiggle_from_trace( struct trace_t* traces,

                         /* These are usually chr names */
                         char** scaffold_names,
                         /* The names of the various tracks */
                         /* Use null for the indexes */
                         char** track_names,
                         
                         const char* output_fname,                           
                         const double filter_threshold )
{    
    FILE* wfp = fopen( output_fname, "w" );
    if( wfp == NULL )
    {
        perror( "FATAL        : Could not open wiggle file for writing " );
        assert( 0 );
        exit( -1 );
    }

    int track_index, j;
    unsigned int k;
    for( track_index = 0; track_index < traces->num_traces; track_index++ )
    {
        /* Print out the header */
        fprintf( wfp, "track type=wiggle_0 name=%s\n", track_names[track_index] );

        for( j = 0; j < traces->num_chrs; j++ )
        {
            /* Print out the new chr start line */
            if( scaffold_names == NULL ) {
                fprintf( wfp, "variableStep chrom=%i\n", j );
            } else {
                fprintf( wfp, "variableStep chrom=%s\n", scaffold_names[j] );
            }
            
            for( k = 0; k < traces->trace_lengths[j]; k++ )
            {
                if( traces->traces[track_index][j][k] > filter_threshold )
                    fprintf( wfp, "%i\t%e\n", k+1, 
                             traces->traces[track_index][j][k] );
            }
        }
    }
    
    fclose( wfp );
}
