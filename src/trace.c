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
#include "mapped_read.h"
#include "iterative_mapping.h"
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
    ts->real_start = real_start;

    ts->length = length;

    ts->data = malloc( sizeof(TRACE_TYPE)*length );
    assert(ts->data != NULL);
    memset( ts->data, 0, sizeof(TRACE_TYPE)*length );

    int rv;
    ts->data_lock = malloc(sizeof(pthread_mutex_t));
    assert(ts->data_lock != NULL);
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
    for(i=0; i < tsegs->num_segments; i++)
    {
        free_trace_segment_t(tsegs->segments + i);
    }
    /* free the array of trace segments */
    free(tsegs->segments);
}

int
find_start_of_trace_segments_to_update(
    struct trace_segments_t* trace_segments,
    int start_bp
)
{
    /* assume the trace_segments are sorted by start position and non-overlapping. */

    /* special case if there's only one trace segment - skip bisect */
    assert(trace_segments->num_segments > 0);
    if(trace_segments->num_segments == 1) {
        return 0;
    }

    /* otherwise, bisect to find the segment containing this location */
    int lo = 0;
    int hi = trace_segments->num_segments;
    while(lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if(trace_segments->segments[mid].real_start < start_bp) {
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
    if(trace_segments->segments[lo].real_start == start_bp) {
        bisect_result = lo;
    } else {
        assert(lo > 0);
        bisect_result = lo - 1;
    }

    return bisect_result;
}

void
init_trace(
        struct genome_data* genome,
        struct trace_t** traces,
        int num_tracks,
        char** track_names
    )
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
        }
    }
    
    return;
}

void
init_full_trace(
        struct genome_data* genome,
        struct trace_t** traces,
        int num_tracks,
        char** track_names
    )
{
    /* Builds the trace metadata with init_trace */
    init_trace( genome, traces, num_tracks, track_names );

    /* Create a single trace segment for each contig in the genome */
    int i, j;
    for( i = 0; i < (*traces)->num_tracks; i++ )
    {
        for( j = 0; j < (*traces)->num_chrs; j++ )
        {
            /* HACK: don't add a segment for the pseudo chromosome */
            if( (*traces)->chr_lengths[j] > 0 )
            {
                add_trace_segment_to_trace_segments(
                    &((*traces)->segments[i][j]),
                    i, j, 0, (*traces)->chr_lengths[j] );
                assert( (*traces)->segments[i][j].num_segments == 1 );
            }
        }
    }
}

void
close_traces( struct trace_t* traces )
{
    int i, j; 
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            free_trace_segments_t( &(traces->segments[i][j]) );
        }
        free( traces->segments[i] );                          
        free( traces->track_names[i] );
    }

    /* free the chr names */
    for( i = 0; i < traces->num_chrs; i++ )
    {
        free( traces->chr_names[i] );
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
                memset( tseg->data, 0, sizeof(TRACE_TYPE)*tseg->length );
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

    int track, chr, segment_index, bp;
    for( track = 0; track < update_trace->num_tracks; track++ )
    {
        for( chr = 0; chr < update_trace->num_chrs; chr++ )
        {
            assert( update_trace->chr_lengths[chr] == other_trace->chr_lengths[chr] );
            
            struct trace_segments_t *update_tsegs
                = &(update_trace->segments[track][chr]);
            struct trace_segments_t *other_tsegs
                = &(other_trace->segments[track][chr]);
                
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

                /* make sure the segments match - this isn't necessary, but
                   simplifies things for now */
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
    int rv;

    /* write out the trace segment metadata */
    rv = fwrite(&(tseg->real_track_id), sizeof(int), 1, os);
    assert(rv == 1);
    rv = fwrite(&(tseg->real_chr_id), sizeof(int), 1, os);
    assert(rv == 1);
    rv = fwrite(&(tseg->real_start), sizeof(int), 1, os);
    assert(rv == 1);

    /* write out the length of the trace */
    rv = fwrite(&(tseg->length), sizeof(int), 1, os);
    assert(rv == 1);

    /* write out the trace */
    rv = fwrite(tseg->data, sizeof(TRACE_TYPE), tseg->length, os);
    assert(rv == tseg->length);
}

void
write_trace_segments_to_stream(
    struct trace_segments_t *tsegs,
    FILE* os )
{
    /* write out the number of segments */
    fwrite(&(tsegs->num_segments), sizeof(int), 1, os );

    /* write out the trace segments */
    int i;
    for( i = 0; i < tsegs->num_segments; i++ )
    {
        write_trace_segment_to_stream(tsegs->segments + i, os);
    }
}

void
write_merged_trace_segments_to_stream(
    struct trace_segments_t* trace_segments,
    FILE* os
)
{
    /* Once iterative mapping is done, it is no longer necessary to mantain
       segmented traces. Therefore, we merge all trace segments into a single
       trace segment when we write the trace out to the file */
    int num_segments = trace_segments->num_segments;

    if( num_segments == 0 )
    {
        /* Write out that there are 0 segments and return (no segments to write) */
        fwrite( &num_segments, sizeof(int), 1, os );
        return;
    }

    assert( trace_segments->num_segments > 0 );
    /* or - if no segments, write out 0 and return? */

    /* For now, we'll maintain the same file format as the previous code, which
       wrote the trace out explicitly. */

    /* write out the number of segments (always 1) */
    num_segments = 1;
    fwrite( &num_segments, sizeof(int), 1, os );

    /* build a new trace segment that merges all of the trace segments */
    int merged_length = 0;
    int i;
    for( i = 0; i < trace_segments->num_segments; i++ )
    {
        merged_length += trace_segments->segments[i].length;
    }

    struct trace_segment_t *merged_segment
        = malloc( sizeof(struct trace_segment_t) );
    init_trace_segment_t( merged_segment,
        trace_segments->segments[0].real_track_id,
        trace_segments->segments[0].real_chr_id,
        0, merged_length );

    /* copy data from each trace segment into the contiguous data field */
    TRACE_TYPE* data_ptr = merged_segment->data;
    for( i = 0; i < trace_segments->num_segments; i++ )
    {
        struct trace_segment_t* current_segment
            = trace_segments->segments + i;
        memcpy( data_ptr, current_segment->data,
            current_segment->length*sizeof(TRACE_TYPE) );
        data_ptr += current_segment->length;
    }

    write_trace_segment_to_stream( merged_segment, os );
    free_trace_segment_t( merged_segment );
    free( merged_segment );
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
            write_merged_trace_segments_to_stream(&(trace->segments[i][j]), os);
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
    rv = fread(&real_track_id, sizeof(int), 1, is);
    assert(rv == 1);

    rv = fread(&real_chr_id, sizeof(int), 1, is);
    assert(rv == 1);

    rv = fread(&real_start, sizeof(int), 1, is);
    assert(rv == 1);

    rv = fread(&length, sizeof(int), 1, is);
    assert(rv == 1);

    /* add the trace segment to set of segments */
    struct trace_segment_t *new_segment
        = add_trace_segment_to_trace_segments(tsegs, real_track_id, real_chr_id, 
            real_start, length);
    /* load the data for this trace segment */
    rv = fread(new_segment->data, sizeof(TRACE_TYPE), new_segment->length, is);
    assert( rv == new_segment->length );
}

void
load_trace_segments_from_stream(
    struct trace_segments_t* tsegs,
    FILE* is )
{
    init_trace_segments_t(tsegs);

    int rv;
    int num_segments;
    rv = fread( &num_segments, sizeof(int), 1, is );
    assert(rv == 1);
    assert(num_segments >= 0);

    /* load the segments */
    int i;
    for( i = 0; i < num_segments; i++ )
    {
        load_trace_segment_from_stream(tsegs, is);
    }

    /* Check that expected number of segments actually loaded */
    assert( num_segments == tsegs->num_segments );
}

void
load_trace_from_stream( struct trace_t** trace, FILE* is )
{
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

void
update_trace_segments_from_mapped_read_array(
    struct trace_segments_t* trace_segments,
    float* update_vals,
    float scale_factor,
    int start,
    int stop
)
{
    /* start, stop are relative to the contig start */
    /* update_vals is from start->stop */

    /* Find the index of the first trace segment to update */
    int st_index = find_start_of_trace_segments_to_update(trace_segments, start);

    int bp = start; // need to maintain last bp updated over multiple segments
    //int update_vals_index = 0;
    /* loop over all the segments that the given mapped fragment overlaps */
    int i;
    for( i = st_index; i < trace_segments->num_segments; i++ )
    {
        struct trace_segment_t* trace_segment = trace_segments->segments + i;
        
        /* lock while we update this trace segment */
        pthread_mutex_lock(trace_segment->data_lock);

        for( ; bp < stop; bp++ )
        {
            /* check if we've reached the end of the current segment - this will
               happen if a mapping spans multiple segments */
            if( bp == trace_segment->real_start + trace_segment->length ) {
                break;
            }
            assert(bp >= trace_segment->real_start);

            float update_val;
            if( update_vals == NULL ) {
                /* update from uniform distribution */
                update_val = 1;
            } else {
                /* update_vals starts at the update value for this first bp in 
                   the trace - convert bp to 0-indexed */
                update_val = update_vals[bp - start];
            }

            /* update the current trace value */
            /* bp - trace_segment->real_start gives index relative to the start
               of the trace segment */
            trace_segment->data[bp - trace_segment->real_start]
                += update_val*scale_factor;
        }

        pthread_mutex_unlock(trace_segment->data_lock);

        /* check if finished updating the given range */
        if(bp-1 == stop)
            break;
    }

    return;
}

/* Wrapper function to update from a uniform kernel */
void
update_trace_segments_from_uniform_kernel(
    struct trace_segments_t* trace_segments,
    float scale_factor,
    int start,
    int stop
)
{
    update_trace_segments_from_mapped_read_array(trace_segments, NULL,
        scale_factor, start, stop);
}

double
accumulate_from_trace(
    const struct trace_t* const traces,
    int track_index,
    int chr_index,
    int start,
    int stop
)
{
    double acc = 0;

    struct trace_segments_t* trace_segments
        = &(traces->segments[track_index][chr_index]);

    int ss_i = find_start_of_trace_segments_to_update(trace_segments, start);

    int bp = start; // need to maintain last bp updated over multiple segments

    int si;
    for( si = ss_i; si < trace_segments->num_segments; si++ )
    {
        struct trace_segment_t* trace_segment = trace_segments->segments + si;

        /* lock while we update this trace segment */
        pthread_mutex_lock(trace_segment->data_lock);

        for( ; bp < stop; bp++ )
        {
            /* check if we've reached the end of the current segment - this will
               happen if a mapping spans multiple segments */
            if( bp == trace_segment->real_start + trace_segment->length ) {
                break;
            }
            assert(bp >= trace_segment->real_start);

            acc += trace_segment->data[bp - trace_segment->real_start];
        }

        pthread_mutex_unlock(trace_segment->data_lock);

        /* check if we've covered the given trace */
        if(bp-1 == stop)
            break;
    }

    return acc;
}

/* Wrapper over accumulate_from_trace - accumulates values from the trace across
   all tracks */
double
accumulate_from_traces(
    const struct trace_t* const traces,
    int chr_index,
    int start,
    int stop
)
{
    double acc = 0;

    int ti;
    for( ti = 0; ti < traces->num_tracks; ti++ )
    {
        acc += accumulate_from_trace(traces, ti, chr_index, start, stop);
    }

    return acc;
}

struct trace_t*
copy_trace_metadata(
    struct trace_t* original
)
{
    /* Return a new trace with the same metadata as the original */
    struct trace_t* new = malloc( sizeof(struct trace_t) );

    new->num_tracks = original->num_tracks;

    /* Copy the track names */
    int i;
    new->track_names = malloc( sizeof(char*)*new->num_tracks );
    for (i = 0; i < new->num_tracks; i++)
    {
        new->track_names[i]
            = malloc( strlen(original->track_names[i])*sizeof(char) );
        strcpy( new->track_names[i], original->track_names[i] );
    }

    /* Copy the chr info */
    new->num_chrs = original->num_chrs;
    new->chr_lengths = malloc( sizeof(unsigned int)*new->num_chrs );
    new->chr_names = malloc( sizeof(char*)*new->num_chrs );
    assert( new->chr_lengths != NULL );
    assert( new->chr_names != NULL );

    for( i = 0; i < new->num_chrs; i++ )
    {
        new->chr_lengths[i] = original->chr_lengths[i];
        new->chr_names[i]
            = malloc( strlen(original->chr_names[i])*sizeof(char) );
        strcpy( new->chr_names[i], original->chr_names[i] );
    }


    /* Allocate pointers to lists of segments for each track */
    new->segments = malloc( sizeof(struct trace_segments_t*)*new->num_tracks );
    assert( new->segments != NULL );

    /* Allocate space for the list of segments for each chr */
    for( i = 0; i < new->num_tracks; i++ )
    {
        new->segments[i]
            = malloc( new->num_chrs*sizeof(struct trace_segments_t) );
        assert( new->segments[i] );

        int j;
        for( j = 0; j < new->num_chrs; j++ )
        {
            init_trace_segments_t( &(new->segments[i][j]) );
        }
    }

    return new;
}

struct segments_list*
init_segments_list()
{
    struct segments_list* sl = malloc( sizeof(struct segments_list) );

    sl->length = 0;
    sl->segments = NULL;

    return sl;
}

void
free_segments_list(
        struct segments_list* sl
    )
{
    if( sl->segments != NULL ) {
        free( sl->segments );
    }
    free( sl );
}

void
add_segment_to_segments_list(
    int track_index,
    int chr_index,

    int start,
    int stop,

    struct segments_list* sl
)
{
    /* Reallocate segments array to add new segment */
    sl->length++;
    sl->segments = realloc( sl->segments, sizeof(struct segment)*sl->length );

    struct segment *new_segment = sl->segments + sl->length - 1;
    new_segment->track_index = track_index;
    new_segment->chr_index = chr_index;
    new_segment->start = start;
    new_segment->stop = stop;

    return;
}

struct segments_list*
build_trace_segments_list(
        struct trace_t* traces
    )
{
    struct segments_list* segments_list = init_segments_list();

    int i, j;
    for( i = 0; i < traces->num_tracks; i++ )
    {
        for( j = 0; j < traces->num_chrs; j++ )
        {
            /* make sure the original trace was not segmented */
            struct trace_segments_t* original_segments
                = &( traces->segments[i][j] );

            /* make sure the original trace was not pseudo, or segmented */
            if( original_segments->num_segments == 0 ) {
                continue;
            }
            assert( original_segments->num_segments == 1 );

            struct trace_segment_t* original_segment
                = original_segments->segments + 0;

            int bp = 0;
            while( bp < original_segment->length )
            {
                int segment_start = bp;
                /* update bp so the segment is at least the minimum length (is
                   bounded by the end of the trace) */
                bp = MIN( bp + MIN_TRACE_SEGMENT_SIZE, original_segment->length );

                for( ; bp < original_segment->length; bp++ )
                {
                    /* increment bp until we find a base with score zero */
                    if( original_segment->data[bp] == 0 )
                        break;
                }

                add_segment_to_segments_list( i, j, segment_start, bp, 
                    segments_list );
            }
        }
    }

    return segments_list;
}

void
log_segments_list(
        struct segments_list* sl
    )
{
    int i;
    for( i = 0; i < sl->length; i++ )
    {
        struct segment* s = sl->segments + i;
        statmap_log( LOG_DEBUG, "Track: %i\tChr:%i\t(%i, %i)",
            s->track_index, s->chr_index, s->start, s->stop );
    }

    return;
}

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