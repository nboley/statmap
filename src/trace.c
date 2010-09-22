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

static FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s\n'", fname );
        assert( false );
        exit( -1 );
    }
    return tmp;
}

/* get the correct mutex index to access the specified position */
/* this is slow - it should probably be done as a macro in hot areas */
int
get_mutex_index( unsigned int bp)
{
    return bp/TM_GRAN;
}

/* build a array to store the density */
TRACE_TYPE*
init_trace_array( size_t size  )
{
    TRACE_TYPE* chr_trace = malloc(sizeof(TRACE_TYPE)*size);
    if( NULL == chr_trace )
    {
        fprintf(stderr, "Can not allocate '%zu' bytes for a trace,\n", sizeof(TRACE_TYPE)*size );
        exit( -1 );
    }
    
    return chr_trace;
}

void
init_trace_locks( struct trace_t* trace )
{
    /* allocate space for the locks */
    #ifdef USE_MUTEX
    trace->locks = malloc(sizeof(pthread_mutex_t**)*trace->num_tracks);
    #else
    trace->locks = malloc(sizeof(pthread_spinlock_t**)*trace->num_tracks);
    #endif
    
    assert( trace->locks != NULL );
    
    /* Allocate space for the pointers to the chr's individual traces and the locks */
    int i;
    for( i = 0; i < trace->num_tracks; i++ )
    {
        /* Store the trace mutex pointers */
        #ifdef USE_MUTEX
        trace->locks[i] = malloc( trace->num_chrs*sizeof(pthread_mutex_t*) );
        #else
        trace->locks[i] = malloc( trace->num_chrs*sizeof(pthread_spinlock_t*) );
        #endif        
        assert( trace->locks[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < trace->num_chrs; j++ )
        {            
            /* initialize the locks */
            int locks_len = trace->chr_lengths[j]/TM_GRAN;
            if( trace->chr_lengths[j] % TM_GRAN > 0 )
                locks_len += 1;
            
            #ifdef USE_MUTEX
            trace->locks[i][j] = malloc( locks_len*sizeof(pthread_mutex_t) );
            #else
            trace->locks[i][j] = malloc( locks_len*sizeof(pthread_spinlock_t) );
            #endif        
            
            int k;
            for( k = 0; k < locks_len; k++ )
            {
                #ifdef USE_MUTEX
                int error = pthread_mutex_init( trace->locks[i][j] + k, 0 );
                #else
                int error = pthread_spin_init( trace->locks[i][j] + k, 0 );
                #endif        
                
                
                if( error != 0 )
                {
                    perror( "Failed to initialize lock in init_trace" );
                    exit( -1 );
                }
            }
        }
    }

}

/* Build mmapped arrays for all of the chrs ***/
void
init_trace( struct genome_data* genome,
            struct trace_t** traces,
            int num_tracks,
            char** track_names )
{
    /* Allocate space for the struct */
    *traces = malloc( sizeof( struct trace_t ) );
    
    /* set the number of traces */
    (*traces)->num_tracks = num_tracks;
    
    int i;
    (*traces)->track_names = malloc( sizeof(char*)*num_tracks );
    for( i = 0; i < num_tracks; i++ )
    {
        (*traces)->track_names[i] = malloc( (strlen(track_names[i])+1)*sizeof(char) );
        strcpy( (*traces)->track_names[i], track_names[i] );
    }
    
    (*traces)->traces = malloc(sizeof(TRACE_TYPE**)*num_tracks);
    assert( (*traces)->traces != NULL );

    /* set the number of chrs */
    (*traces)->num_chrs = genome->num_chrs;
    (*traces)->chr_lengths = malloc(sizeof(unsigned int)*genome->num_chrs);
    (*traces)->chr_names = malloc(sizeof(char*)*genome->num_chrs);
    assert( (*traces)->chr_lengths != NULL );

    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->chr_lengths[i] = genome->chr_lens[i];

        (*traces)->chr_names[i] = malloc( (strlen(genome->chr_names[i])+1)*sizeof(char) );
        strcpy( (*traces)->chr_names[i], genome->chr_names[i] );
        
    }
    
    /* Allocate space for the pointers to the chr's individual traces and the locks */
    for( i = 0; i < num_tracks; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->traces[i] = malloc( (*traces)->num_chrs*sizeof(TRACE_TYPE*) );
        assert( (*traces)->traces[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {            
            /* initialize the forward trace */
            (*traces)->traces[i][j] = init_trace_array( (*traces)->chr_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memset( (*traces)->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(*traces)->chr_lengths[j] );            
        }
    }
    
    init_trace_locks( *traces );
    
    return;
}

void
copy_trace_data( struct trace_t* traces,
                 struct trace_t* original )
{
    assert( original->num_tracks == traces->num_tracks );
    int i;
    for( i = 0; i < original->num_tracks; i++ )
    {
        /* Store the pointers for the chrs */
        assert( original->num_chrs == traces->num_chrs );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < traces->num_chrs; j++ )
        {            
            assert( traces->chr_lengths[j] == original->chr_lengths[j] );
            
            memcpy( traces->traces[i][j], original->traces[i][j],
                    sizeof(TRACE_TYPE)*(traces->chr_lengths[j])   );
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
    (*traces)->num_tracks = original->num_tracks;
    (*traces)->traces = malloc(sizeof(TRACE_TYPE**)*original->num_tracks);
    assert( (*traces)->traces != NULL );

    (*traces)->num_chrs = original->num_chrs;
    (*traces)->chr_lengths = malloc(sizeof(unsigned int)*original->num_chrs);
    assert( (*traces)->chr_lengths != NULL );

    int i;
    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->chr_lengths[i] = original->chr_lengths[i];
    }

    /* allocate space for the locks */
    #ifdef USE_MUTEX
    (*traces)->locks = malloc(sizeof(pthread_mutex_t**)*original->num_tracks);
    #else
    (*traces)->locks = malloc(sizeof(pthread_spinlock_t**)*original->num_tracks);
    #endif

    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < original->num_tracks; i++ )
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
            (*traces)->traces[i][j] = init_trace_array( (*traces)->chr_lengths[j] );
            
            assert( (*traces)->traces[i][j] != NULL );
            memcpy( (*traces)->traces[i][j], original->traces[i][j],
                    sizeof(TRACE_TYPE)*(*traces)->chr_lengths[j]   );

            /* initialize the locks */
            int locks_len = (*traces)->chr_lengths[j]/TM_GRAN;
            if( (*traces)->chr_lengths[j] % TM_GRAN > 0 )
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
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            free( traces->traces[i][j] );                          

            /* initialize the locks */
            int locks_len = traces->chr_lengths[j]/TM_GRAN;
            if( traces->chr_lengths[j] % TM_GRAN > 0 )
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

            // BUG ( memory leak ) should this be freed??!?!?
            // free( traces->locks[i][j] );
        }
        free( traces->traces[i] );                          
        free( traces->locks[i] );
    }

    free( traces->chr_lengths );
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
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->chr_lengths[j]; k++ )
            {
                traces->traces[i][j][k] /= value;
            }
        }
    }
    return;
}

void
multiply_trace_by_scalar( struct trace_t* traces, double value )
{
    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->chr_lengths[j]; k++ )
            {
                traces->traces[i][j][k] *= value;
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
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->chr_lengths[j]; k++ )
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
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            memset( traces->traces[i][j], 0, 
                    sizeof(TRACE_TYPE)*(traces->chr_lengths[j]) );
        }
    }
}

void
set_trace_to_uniform( struct trace_t* traces, double value )
{
    /* Zero out the trace for the update */
    int i, j;
    unsigned int k;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->chr_lengths[j]; k++ )
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
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            for( k = 0; k < traces->chr_lengths[j]; k++ )
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
    assert( update_trace->num_tracks = other_trace->num_tracks );
    assert( update_trace->num_chrs = other_trace->num_chrs );
    for( trace_index = 0; trace_index < update_trace->num_tracks; trace_index++ )
    {
        for( chr = 0; chr < update_trace->num_chrs; chr++ )
        {
            assert( update_trace->chr_lengths[chr] == other_trace->chr_lengths[chr] );
            for( bp = 0; bp < update_trace->chr_lengths[chr]; bp++ )
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

/******************************************************************
 *
 * Wiggle output code
 *
 */

extern void
write_wiggle_from_trace_to_stream( 
    struct trace_t* traces,
    
    /* These are usually chr names */
    // char** scaffold_names,
    /* The names of the various tracks */
    /* Use null for the indexes */
    // char** track_names,
    
    FILE* os, /* output stream */                           
    const double filter_threshold )
{    
    unsigned int global_counter = 0;

    int track_index, j;
    unsigned int k;
    for( track_index = 0; track_index < traces->num_tracks; track_index++ )
    {
        /* Print out the header */
        fprintf( os, "track type=wiggle_0 name=%s\n", 
                 traces->track_names[track_index] );

        for( j = 0; j < traces->num_chrs; j++ )
        {
            /* skip the pseudo chr */
            if( j == PSEUDO_LOC_CHR_INDEX )
                continue;
            
            /* Print out the new chr start line */
            if( traces->chr_names == NULL ) {
                fprintf( os, "variableStep chrom=%i\n", j );
            } else {
                fprintf( os, "variableStep chrom=%s\n", traces->chr_names[j] );
            }
            
            for( k = 0; k < traces->chr_lengths[j]; k++ )
            {
                global_counter += 1;
                
                if( traces->traces[track_index][j][k] > filter_threshold )
                    fprintf( os, "%i\t%e\n", k+1, 
                             traces->traces[track_index][j][k] );
                
                if( global_counter > 0  && global_counter%10000000 == 0 )
                    fprintf( stderr, "NOTICE        :  Written %i positions to trace.\n", global_counter );
            }
        }
    }
    
    return;
}

extern void
write_wiggle_from_trace( struct trace_t* traces,
                         const char* output_fname,                           
                         const double filter_threshold )
{    
    FILE* wfp = fopen( output_fname, "a" );
    if( wfp == NULL )
    {
        perror( "FATAL        : Could not open wiggle file for writing " );
        fprintf( stderr, "Filename: %s\n", output_fname );
        assert( 0 );
        exit( -1 );
    }
    
    write_wiggle_from_trace_to_stream( traces, wfp, filter_threshold );
    
    fclose( wfp );
}


/******************************************************************
 *
 * Trace marshalling/unmarshalling code
 *
 */

static void
write_trace_header_to_stream( struct trace_t* trace, FILE* os )
{
    /* a counter for for loops  */
    int i;
    
    /* write the magic number */
    fwrite( TRACE_MAGIC_NUMBER, 6, sizeof(char), os );
    
    /** write the number of tracks */
    fwrite( &trace->num_tracks, sizeof(int), 1, os );

    /** write the track data */
    /* write the track name lengths */
    for( i = 0; i < trace->num_tracks; i++ )
    {
        size_t track_name_len = strlen(trace->track_names[i]);
        fwrite( &track_name_len, sizeof(size_t), 1, os );    
    }

    /* write the track names */
    for( i = 0; i < trace->num_tracks; i++ )
        fwrite( trace->track_names[i], sizeof(char), strlen(trace->track_names[i]), os );    
    
    /** write the chromosomes data */
    /* write the num of chromosomes */
    fwrite( &trace->num_chrs, sizeof(int), 1, os );    

    /* write the chr name lengths */
    for( i = 0; i < trace->num_chrs; i++ )
    {
        size_t chr_name_len = strlen(trace->chr_names[i]);
        fwrite( &chr_name_len, sizeof(size_t), 1, os );    
    }

    /* write the chr names */
    for( i = 0; i < trace->num_chrs; i++ )
        fwrite( trace->chr_names[i], sizeof(char), strlen(trace->chr_names[i]), os );    
    
    /* write the chr lengths */
    fwrite( trace->chr_lengths, sizeof(unsigned int), trace->num_chrs, os );    

    return;
}

static void
load_trace_header_from_stream( struct trace_t* trace, FILE* is )
{
    /* return values of gets */
    int rv;
    
    /* a counter for for loops  */
    int i;
    
    /* write the magic number */
    char MN[6];
    rv = fread( MN, sizeof(char), 6, is );
    assert( 6 == rv );
    assert( 0 == strcmp( MN, TRACE_MAGIC_NUMBER ) );
    
    /** find the number of tracks */
    rv = fread( &(trace->num_tracks), sizeof(int), 1, is );
    assert( 1 == rv );
    
    /** read the track data */
    /* read the track name lengths */
    size_t track_name_lengths[trace->num_tracks];
    rv = fread( track_name_lengths, sizeof(size_t), trace->num_tracks, is );
    assert( trace->num_tracks == rv );

    trace->track_names = calloc( sizeof(char*), trace->num_tracks );
    for( i = 0; i < trace->num_tracks; i++ )
    {
        trace->track_names[i] = calloc( sizeof(char), track_name_lengths[i] + 1 );
        rv = fread( trace->track_names[i], sizeof(char), track_name_lengths[i], is );
        /* make sure that the terminating NULL is present */
        trace->track_names[i][track_name_lengths[i]] = '\0';
        assert( track_name_lengths[i] == (unsigned int) rv );
    }
    
    /** read the chromosomes data */
    /* read the num of chromosomes */
    rv = fread( &(trace->num_chrs), sizeof(int), 1, is );
    assert( rv == 1 );

    /* read the chr name lengths */
    size_t chr_name_lengths[trace->num_chrs];
    rv = fread( chr_name_lengths, sizeof(size_t), trace->num_chrs, is );
    assert( trace->num_chrs == rv );

    /* read the chr names */
    trace->chr_names = calloc( sizeof(char*), trace->num_chrs );
    for( i = 0; i < trace->num_chrs; i++ )
    {
        trace->chr_names[i] = calloc( sizeof(char), chr_name_lengths[i] + 1 );
        rv = fread( trace->chr_names[i], sizeof(char), chr_name_lengths[i], is );
        /* sxplicitly set the NULL terminating char */
        trace->chr_names[i][chr_name_lengths[i]] = '\0';
        assert( chr_name_lengths[i] == (unsigned int) rv );
    }
    
    /* read the chr lengths */
    trace->chr_lengths = calloc( sizeof(unsigned int), trace->num_chrs );
    rv = fread( trace->chr_lengths, sizeof(unsigned int), trace->num_chrs, is );
    assert( trace->num_chrs == rv );

    return;
}

void
write_trace_to_stream( struct trace_t* trace, FILE* os )
{
    write_trace_header_to_stream( trace, os );

    int i;
    for( i = 0; i < trace->num_tracks; i++ )
    {
        int j;
        for( j = 0; j < trace->num_chrs; j++ )
        {
            /* write the array from track i, chr j */
            fwrite( trace->traces[i][j], 
                    sizeof(TRACE_TYPE), trace->chr_lengths[j],
                    os );
        }
    }
    
    return;
}

void
write_trace_to_file( struct trace_t* trace, char* fname )
{
    FILE* fp = fopen( fname, "w" );
    write_trace_to_stream( trace, fp );
    fclose(fp);
}

void
load_trace_from_stream( struct trace_t** trace, FILE* is )
{
    *trace = malloc( sizeof( struct trace_t ) );
    
    load_trace_header_from_stream( *trace, is );

    init_trace_locks( *trace );

    int i;
    (*trace)->traces = calloc( sizeof(TRACE_TYPE*), (*trace)->num_tracks  );
    for( i = 0; i < (*trace)->num_tracks; i++ )
    {
        int j;
        (*trace)->traces[i] = calloc( sizeof(TRACE_TYPE*), (*trace)->num_chrs  );
        for( j = 0; j < (*trace)->num_chrs; j++ )
        {
            /* read the array from track i, chr j */            
            (*trace)->traces[i][j] = 
                malloc( sizeof(TRACE_TYPE)*((*trace)->chr_lengths[j]) );
            fread( (*trace)->traces[i][j], 
                    sizeof(TRACE_TYPE), (*trace)->chr_lengths[j],
                    is );
        }
    }
    
    return;
}

void
load_trace_from_file( struct trace_t** trace, char* fname )
{
    FILE* fp = open_check_error( fname, "r" );
    load_trace_from_stream( trace, fp );
    fclose(fp);
}
