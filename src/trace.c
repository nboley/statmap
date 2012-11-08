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
#include <errno.h>
#include <float.h>
#include <pthread.h>

#include "config.h"
#include "trace.h"
#include "genome.h"
#include "util.h"
#include "log.h"

void
init_trace_segment_t(
        struct trace_segment_t *ts,
        int real_track_id,
        int real_chr_id,
        int real_start,
        int length
    )
{
    /* ts is a pointer to an allocated trace_segment_t */
    assert(ts != NULL);

    ts->real_track_id = real_track_id;
    ts->real_chr_id = real_chr_id;
    ts->real_start = real_start;'
    '
    ts->length = length;

    ts->data = malloc(sizeof(TRACE_TYPE)*length);
    assert(ts->data != NULL);
    memset(ts->data, 0, sizeof(TRACE_TYPE)*length);

    int rv;
    rv = pthread_mutex_init(ts->data_lock, NULL); // NULL uses attr defaults
    assert(rv == 0);
}

void
free_trace_segment_t(
        struct trace_segment_t* ts
    )
{
    /* for now, require our design to be so clean that we never try to free an
       unallocated trace_segment. Can be relaxed (and return on NULL) if this
       utopian vision is not achieved. */
    assert(ts != NULL);

    /* free memory allocated in the trace_segment_t */
    free(ts->data);

    /* free the lock */
    int rv;
    rv = pthread_mutex_destroy(ts->data_lock);
    assert(rv == 0);
    free((void*)ts->data_lock);
}

struct trace_segment_t*
add_trace_segment_to_trace_segments(
        struct trace_segments_t* segs,
        int real_track_id,
        int real_chr_id,
        int real_start,
        int length
    )
{
    segs->num_segments++;
    segs->segments = realloc(segs->segments,
        sizeof(struct trace_segment_t)*segs->num_segments);
    assert(segs->segments != NULL);

    /* initialize the new trace segment with the given values */
    struct trace_segment_t* new_tseg = segs->segments + segs->num_segments-1;
    init_trace_segment_t(new_tseg, real_track_id, real_chr_id, real_start,
        length);

    return new_tseg;
}

void
init_trace_segments_t(
        struct trace_segments_t *segs
    )
{
    /* segs is already allocated as part of a contiguous array */
    assert(segs != NULL);

    segs->num_segments = 0;
    segs->segments = NULL;
}

void
free_trace_segments_t(
    struct trace_segments_t* tsegs
)
{
    /* free memory allocated in the trace segments */
    int i;
    for(i=0; i < tsegs->num_segments; tsegs++)
    {
        free_trace_segment_t(tsegs->segments + i);
    }
    /* free the array of trace segments */
    free(tsegs->segments);
}

struct trace_segment_t*
find_trace_segment_t(
    struct trace_t* traces,
    int track_index,
    int chr_index,
    int bp
)
{
    struct trace_segments_t* trace_segments
        = &(traces->segments[track_index][chr_index]);

    /* assume the trace_segments are sorted by start position and non-overlapping. */

    /* special case if there's only one trace segment - skip bisect */
    assert(trace_segments->num_segments > 0);
    if(trace_segments->num_segments == 1) {
        return trace_segments->segments; // + 0
    }

    /* otherwise, bisect to find the segment containing this location */
    int lo = 0
    int hi = trace_segments->num_segments;
    while(lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if trace_segments->segments[mid].real_start < bp {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    /* make sure the binary search is working */
    assert(lo <= hi);
    assert(lo <= trace_segments->num_segments);
    assert(lo >= 0);

    int bisect_result;
    if(trace_segments->segments[lo].real_start == bp) {
        bisect_result = lo;
    } else {
        assert(lo > 0);
        bisect_result = lo - 1;
    }

    return trace_segments->segments + bisect_result;
}

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
        (*traces)->track_names[i]
            = malloc( (strlen(track_names[i])+1)*sizeof(char) );
        strcpy( (*traces)->track_names[i], track_names[i] );
    }
    
    /* Allocate pointers to lists of segments for each track */
    (*traces)->segments = malloc(sizeof(struct trace_segments_t*)*num_tracks);
    assert( (*traces)->segments != NULL );

    /* set the number of chrs */
    (*traces)->num_chrs = genome->num_chrs;
    (*traces)->chr_lengths = malloc(sizeof(unsigned int)*genome->num_chrs);
    (*traces)->chr_names = malloc(sizeof(char*)*genome->num_chrs);
    assert( (*traces)->chr_lengths != NULL );
    assert( (*traces)->chr_names != NULL );

    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        (*traces)->chr_lengths[i] = genome->chr_lens[i];
        (*traces)->chr_names[i]
            = malloc( (strlen(genome->chr_names[i])+1)*sizeof(char) );
        strcpy( (*traces)->chr_names[i], genome->chr_names[i] );
    }
    
    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < num_tracks; i++ )
    {
        /* Store the pointers for the chrs */
        (*traces)->segments[i]
            = malloc( (*traces)->num_chrs*sizeof(struct trace_segments_t) );
        assert( (*traces)->segments[i] != NULL );
        
        /* initialize each chr */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {
            /* initialize the trace segments object for each chr */
            init_trace_segments_t( &((*traces)->segments[i][j]) );
            /* for now, initialize a single trace_segment_t for the entire contig */
            add_trace_segment_to_trace_segments( &((*traces)->segments[i][j]),
                i, j, 0. (*traces)->chr_lengths[j] );
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
    /* TODO copy track_names ? the original code here didn't */

    (*traces)->segments = malloc(
        sizeof(struct trace_segments_t*)*original->num_tracks);
    assert( (*traces)->segments != NULL );

    (*traces)->num_chrs = original->num_chrs;
    (*traces)->chr_lengths = malloc(sizeof(unsigned int)*original->num_chrs);
    assert( (*traces)->chr_lengths != NULL );

    int i;
    for( i = 0; i < (*traces)->num_chrs; i++ )
    {
        /* copy the trace lengths, in bps */
        (*traces)->chr_lengths[i] = original->chr_lengths[i];
        /* TODO copy chromsome names? the original code here didn't */
    }

    /* Allocate space for the pointers to the chr's individual traces */
    for( i = 0; i < original->num_tracks; i++ )
    {
        /* Allocate an array of trace segments objects for each chromosome */
        (*traces)->segments[i] = malloc( 
            (*traces)->num_chrs*sizeof(struct trace_segments_t));
        assert( (*traces)->segments[i] != NULL );
        
        /* copy the trace_segments_t */
        int j;
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {
            struct trace_segments_t *new_tsegs = &((*traces)->segments[i][j]);
            init_trace_segments_t(new_tsegs);
            assert(new_tsegs != NULL);

            /* copy each trace segment */
            struct trace_segments_t *original_tsegs = &(original->segments[i][j]);
            int k;
            for(k = 0; k < original_tsegs->num_segments; k++)
            {
                /* TODO this interface could maybe be cleaner - need to
                   1) realloc the contiguous array of trace_segments to add the
                      new one
                   2) initialize the trace_segment struct (mostly just need to
                      allocate memory for ts->data)
                   3) copy other fields, and possibly data
                */
                struct trace_segment_t *original_tseg 
                    = original_tsegs->segments + k;
                struct trace_segment_t *new_tseg
                    = add_trace_segment_to_trace_segments(new_tsegs,
                        original_tseg->real_track_id, original_tseg->real_chr_id,
                        original_tseg->real_start, original_tseg->length);

                /* copy data from original trace segment */
                memcpy(new_tseg->data, original_tseg->data, sizeof(TRACE_TYPE)*
                    original_tseg->length);
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
            free_trace_segments_t(&(traces->segments[i][j]));
            free( traces->chr_names[j] );
        }
        free( traces->segments[i] );                          
        free( traces->track_names[i] );
    }

    free( traces->track_names );
    free( traces->chr_names );
    free( traces->chr_lengths );
    free( traces->segments );
    
    free( traces );

    return;
}

void
divide_trace_by_sum( struct trace_t* traces, double value )
{
    int i, j, k, l;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            /* Loop over each trace segment in the set of trace segments for
               this track, chr combination */
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                for( l = 0; l < tseg->length; l++ )
                {
                    tseg->data[l] /= value;
                }
            }
        }
    }
    return;
}

void
multiply_trace_by_scalar( struct trace_t* traces, double value )
{
    int i, j, k, l;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                for( l = 0; l < tseg->length; l++ )
                {
                    tseg->data[l] *= value;
                }
            }
        }
    }
    return;
}

double
sum_traces( struct trace_t* traces )
{
    double sum = 0;
    int i, j, k, l;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                for( l = 0; l < tseg->length; l++ )
                {
                    sum += tseg->data[l];
                }
            }
        }
    }
    
    return sum;
}

void
normalize_traces( struct trace_t* traces )
{
    double sum = sum_traces( traces );
    assert( sum > 0 );
    divide_trace_by_sum( traces, sum );
    return;
}

void
zero_traces( struct trace_t* traces )
{
    /* Zero out the trace for the update */
    int i, j, k;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                memset(tseg->data, 0, sizeof(TRACE_TYPE)*tseg->length);
            }
        }
    }
}

void
set_trace_to_uniform( struct trace_t* traces, double value )
{
    int i, j, k, l;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                for( l = 0; l < tseg->length; l++ )
                {
                    tseg->data[l] = value;
                }
            }
        }
    }
}

void
apply_to_trace( struct trace_t* traces, double (*fun)(double) )
{
    int i, j, k, l;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            struct trace_segments_t *tsegs = &(traces->segments[i][j]);
            for( k = 0; k < tsegs->num_segments; k++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + k;
                for( l = 0; l < tseg->length; l++ )
                {
                    tseg->data[l] = fun(tseg->data[l]);
                }
            }            
        }
    }
}

/* traces must be the same dimension */
/* This function applies an aggregate to every basepair in the traces */
void
aggregate_over_trace_pairs(  struct trace_t* update_trace, 
                             const struct trace_t* const other_trace,
                             TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
                          )
{
    assert( update_trace->num_tracks = other_trace->num_tracks );
    assert( update_trace->num_chrs = other_trace->num_chrs );

    int trace_index, chr, segment_index, bp;
    for( trace_index = 0; trace_index < update_trace->num_tracks; trace_index++ )
    {
        for( chr = 0; chr < update_trace->num_chrs; chr++ )
        {
            assert( update_trace->chr_lengths[chr] == other_trace->chr_lengths[chr] );
            
            struct trace_segments_t *update_tsegs
                = &(update_trace->segments[trace_index][chr]);
            struct trace_segments_t *other_tsegs
                = &(other_trace->segments[trace_index][chr]);
            /* make sure both traces have the same sets of segments */
            assert( update_tsegs->num_segments == other_tsegs->num_segments );

            for( segment_index = 0;
                 segment_index < update_tsegs->num_segments;
                 segment_index++ )
            {
                struct trace_segment_t *update_tseg
                    = update_tsegs->segments + segment_index;
                struct trace_segment_t *other_tseg
                    = other_tsegs->segments + segment_index;

                /* make sure the segments match */
                assert( update_tseg->length == other_tseg->length );

                for( bp = 0; bp < update_tseg->length; bp++ )
                {
                    update_tseg->data[bp] = aggregate(update_tseg->data[bp],
                        other_tseg->data[bp] );
                }
            }
        }
    }

    return;
}

void
aggregate_over_traces_from_fnames(  
    char** fnames,
    int num_files,
    struct trace_t** agg_trace,
    TRACE_TYPE (*aggregate)( const TRACE_TYPE, const TRACE_TYPE )
)
{
    /* initalize the trace that we will aggregate over from the first file */    
    load_trace_from_file( agg_trace, fnames[0] );
    
    int i;
    for( i = 1; i < num_files; i++ )
    {
        char* fname = fnames[i];
        
        struct trace_t* curr_trace;
        load_trace_from_file( &curr_trace, fname );
        aggregate_over_trace_pairs( *agg_trace, curr_trace, aggregate );
        close_traces( curr_trace );
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
    FILE* os, /* output stream */                           
    const double filter_threshold )
{    
    unsigned int global_counter = 0;

    int track, chr, segment, bp;
    for( track = 0; track < traces->num_tracks; track++ )
    {
        /* Print out the header */
        fprintf( os, "track type=wiggle_0 name=%s\n", 
                 traces->track_names[track] );

        for( chr = 0; chr < traces->num_chrs; chr++ )
        {
            /* skip the pseudo chr */
            if( chr == PSEUDO_LOC_CHR_INDEX )
                continue;
            
            /* Print out the new chr start line */
            if( traces->chr_names == NULL ) {
                fprintf( os, "variableStep chrom=%i\n", chr );
            } else {
                fprintf( os, "variableStep chrom=%s\n", traces->chr_names[chr] );
            }
            
            struct trace_segments_t *tsegs = &(traces->segments[track][chr]);
            for( segment = 0; segment < tsegs->num_segments; segment++ )
            {
                struct trace_segment_t *tseg = tsegs->segments + segment;
                for( bp = 0; bp < tseg->length; bp++ )
                {
                    global_counter += 1;
                
                    if( tseg->data[bp] > filter_threshold )
                        fprintf( os, "%i\t%e\n", bp+1, tseg->data[bp] );
                    
                    if( global_counter > 0  && global_counter%10000000 == 0 )
                        statmap_log(LOG_NOTICE, "Written %i positions to trace.",
                            global_counter );
                }
            }
        }
    }
    
    return;
}

/* Wrapper to write wiggle to stdout from Python API */
extern void
write_wiggle_from_trace_to_stdout(
    struct trace_t* traces,
    const double filter_threshold
)
{
    write_wiggle_from_trace_to_stream( traces, stdout, filter_threshold );
}

extern void
write_wiggle_from_trace( struct trace_t* traces,
                         const char* output_fname,                           
                         const double filter_threshold )
{    
    FILE* wfp = open_check_error( output_fname, "a" );
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
    char MN[7];
    memset( MN, 0, sizeof(char)*7 );
    rv = fread( MN, sizeof(char), 6, is );
    if( 6 != rv || 0 != strcmp( MN, TRACE_MAGIC_NUMBER ))
    {
        statmap_log( LOG_FATAL, "Mismatched bin trace header ( rv: %i, header: '%s' )",  rv, MN  );
        assert( 0 );
    }
    
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
write_trace_segment_to_stream(
    struct trace_segment_t *tseg,
    FILE* os )
{
    /* write out the trace segment metadata */
    fwrite(tseg->real_track_id, sizeof(int), 1, os);
    fwrite(tseg->real_chr_id, sizeof(int), 1, os);
    fwrite(tseg->real_start, sizeof(int), 1, os);

    /* write out the length of the trace */
    fwrite(tseg->length, sizeof(int), 1, os);

    /* write out the trace */
    fwrite(tseg->data, sizeof(TRACE_TYPE), tseg->length, os);
}

void
write_trace_segments_to_stream(
    struct trace_segments_t *tsegs,
    FILE* os )
{
    /* write out the number of segments */
    fwrite( tsegs->num_segments, sizeof(int), 1, os );

    /* write out the trace segments */
    int i;
    for( i = 0; i < tsegs->num_segments; i++ )
    {
        write_trace_segment_to_stream(tsegs->segments + i, os);
    }
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
            write_trace_segments_to_stream(&(trace->segments[i][j]), os);
        }
    }
    
    return;
}

void
write_trace_to_file( struct trace_t* trace, char* fname )
{
    FILE* fp = open_check_error( fname, "w" );
    write_trace_to_stream( trace, fp );
    fclose(fp);
}

void
load_trace_segment_from_stream(
    /* append the loaded trace segment to tsegs */
    struct trace_segments_t* tsegs,
    FILE* is )
{
    int rv;

    /* load the metadata for segment */
    int real_track_id, real_chr_id, real_start, length;
    fread(&real_track_id, sizeof(int), 1, is);
    assert(rv == 1);

    fread(&real_chr_id, sizeof(int), 1, is);
    assert(rv == 1);

    fread(&real_start, sizeof(int), 1, is);
    assert(rv == 1);

    fread(&length, sizeof(int), 1, is);
    assert(rv == 1);

    /* add the trace segment to set of segments */
    struct trace_segment_t *new_segment
        = add_trace_segment_to_trace_segments(tsegs, real_track_id, real_chr_id, 
            real_start, length);
    /* load the data for this trace segment */
    rv = fread(&(new_segment->data), sizeof(TRACE_TYPE), new_segment->length, is);
    assert( rv == new_segment->length );
}

void
load_trace_segments_from_stream(
    struct trace_segments_t* tsegs,
    FILE* is )
{
    int rv;

    init_trace_segments_t(tsegs);

    /* load the number of segments from the stream */
    rv = fread(&(tsegs->num_segments), sizeof(int), 1, is);
    assert(rv == 1);
    assert(tsegs->num_segments > 0);

    /* load the segments */
    int i;
    for( i = 0; i < tsegs->num_segments; i++ )
    {
        load_trace_segment_from_stream(tsegs, is);
    }
}

void
load_trace_from_stream( struct trace_t** trace, FILE* is )
{
    size_t rv;
    
    *trace = malloc( sizeof( struct trace_t ) );
    
    load_trace_header_from_stream( *trace, is );

    (*trace)->segments 
        = malloc(sizeof(struct trace_segments_t*)*(*trace)->num_tracks);
    int i, j;
    for( i = 0; i < (*trace)->num_tracks; i++ )
    {
        (*trace)->segments[i]
            = malloc(sizeof(struct trace_segments_t)*(*trace)->num_chrs);
        
        for( j = 0; j < (*trace)->num_chrs; j++ )
        {
            load_trace_segments_from_stream(&((*trace)->segments[i][j]), is);
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

/********************************************************************
 * Trace aggregate functions
 *
 ********************************************************************/

TRACE_TYPE
trace_agg_sum( const TRACE_TYPE a, const TRACE_TYPE b )
{
    return a+b;
}

TRACE_TYPE
trace_agg_min( const TRACE_TYPE a, const TRACE_TYPE b )
{
    if ( a < b )
        return a;
    
    return b;
}

TRACE_TYPE
trace_agg_max( const TRACE_TYPE a, const TRACE_TYPE b )
{
    if ( a > b )
        return a;
    
    return b;
}
