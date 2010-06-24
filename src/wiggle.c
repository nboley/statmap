/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "wiggle.h"
#include "trace.h"

#include "genome.h"
#include "mapped_read.h"

// move this 
void
naive_update_trace_expectation_from_location( 
    const struct trace_t* const traces, 
    const struct mapped_read_location* const loc )
{
    unsigned int i;
    
    const int chr_index = loc->chr;
    const unsigned char flag = loc->flag; 
    const unsigned int start = loc->start_pos;
    const unsigned int stop = loc->stop_pos;
    const ML_PRB_TYPE cond_prob = loc->cond_prob;

    if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
    {
        /* lock the locks */
        for( i = start/TM_GRAN; i <= stop/TM_GRAN; i++ )
        {
            #ifdef USE_MUTEX
            pthread_mutex_lock( traces->locks[1][ chr_index ] + i );
            #else
            pthread_spin_lock( traces->locks[1][ chr_index ] + i );
            #endif
        }
        
        for( i = start; i <= stop; i++ )
            traces->traces[1][chr_index][i] += cond_prob;

        /* unlock the locks */
        for( i = start/TM_GRAN; i <= stop/TM_GRAN; i++ )
        {
            #ifdef USE_MUTEX
            pthread_mutex_unlock( traces->locks[1][ chr_index ] + i );
            #else
            pthread_spin_unlock( traces->locks[1][ chr_index ] + i );
            #endif
        }
    } else {
        /* lock the locks */
        for( i = start/TM_GRAN; i <= stop/TM_GRAN; i++ )
        {
            #ifdef USE_MUTEX
            pthread_mutex_lock( traces->locks[0][ chr_index ] + i );
            #else
            pthread_spin_lock( traces->locks[0][ chr_index ] + i );
            #endif
        }
        
        for( i = start; i <= stop; i++ )
            traces->traces[0][chr_index][i] += cond_prob;

        /* unlock the locks */
        for( i = start/TM_GRAN; i <= stop/TM_GRAN; i++ )
        {
            #ifdef USE_MUTEX
            pthread_mutex_unlock( traces->locks[0][ chr_index ] + i );
            #else
            pthread_spin_unlock( traces->locks[0][ chr_index ] + i );
            #endif
        }
    }
    
    return;
}


float 
wig_lines_min( const struct wig_line_info* lines, const int ub, const int num_wigs  )
{
    assert( num_wigs > 0 );
    
    float min = lines[0].value;
    int i;
    for( i = 1; i <= ub; i++ )
    {
        if( lines[i].value < min )
            min = lines[i].value;
    }
    
    if( (ub < num_wigs-1) && (min > 0) )
        min = 0;

    return min;
}

float 
wig_lines_max( const struct wig_line_info* lines, const int ub, const int num_wigs  )
{
    assert( num_wigs > 0 );
    
    float max = lines[0].value;
    int i;
    for( i = 1; i <= ub; i++ )
    {
        if( lines[i].value > max )
            max = lines[i].value;
    }

    if( (ub < num_wigs-1) && (max < 0) )
        max = 0;
    
    return max;
}

float 
wig_lines_sum( const struct wig_line_info* lines, const int ub, const int num_wigs  )
{
    assert( num_wigs > 0 );
    
    float sum = 0;
    int i;
    for( i = 0; i <= ub; i++ )
        sum += lines[i].value;
    
    return sum;
}


/* 
 *  This makes tons of assumptions about the 
 *  wiggle file format, order, that wont hold in 
 *  general. IE, this is a poor general solution,
 *  but it is fine to be used on files generated 
 *  by statmap 
 */

static int 
cmp_wig_line_info( const void* a, const void* b)
{
    /* compare on the file pointer. if it is null, the file is empty, 
       and it is always smaller */
    if( ((struct wig_line_info*) a)->fp == NULL )
    {
        if( ((struct wig_line_info*) b)->fp == NULL )
        {
            return 0;
        } else {
            return -1;
        }
    }
    if( ((struct wig_line_info*) b)->fp == NULL )
    {
        return 1;
    }

    /* compare on chr index */
    if( ((struct wig_line_info*) b)->chr_index != ((struct wig_line_info*) a)->chr_index )
    {
        return ((struct wig_line_info*) a)->chr_index - ((struct wig_line_info*) b)->chr_index;
    }
    
    return ((struct wig_line_info*) a)->position - ((struct wig_line_info*) b)->position;
}

static void
parse_next_line( struct wig_line_info* lines, char** chr_names, int index )
{
    char* rv;
    char buffer[500];

    /* loop until we get a valid line */
    while( 1 )
    {
        rv = fgets( buffer, 500, lines[index].fp );
        if( rv == NULL ) 
        {
            lines[index].fp = NULL;
            return;
        }
        if( buffer[0] == 'f' )
        {
            fprintf(stderr, "FATAL    : Wiggle parser does not support fixed step lines\n");
            assert( 0 );
            exit( -1 );
        } else if ( buffer[0] == 'v') {
            /* increment the chr index */
            lines[index].chr_index += 1;
            /* get the chr name */
            char* chr_name = strstr( buffer, "chrom=" );
            chr_name += 6;
            /* set the chr name */
            if( chr_names[lines[index].chr_index] != NULL )
            {
                if( 0 != strncmp( chr_name, chr_names[lines[index].chr_index], strlen( chr_names[lines[index].chr_index] ) ) )
                {
                    fprintf( stderr, "ERROR     : Chr names ( new: %.*s and old: %s ) are out of sync", 
                             (int)strlen(chr_names[lines[index].chr_index]), chr_name, chr_names[lines[index].chr_index] );
                    assert( 0 );
                }
            } else {
                /* remoive the trailing newline */
                int chr_name_len = strlen(chr_name)-1;
                chr_names[lines[index].chr_index] = calloc( chr_name_len+1, sizeof(char)  );
                memcpy( chr_names[lines[index].chr_index], chr_name, chr_name_len );
                chr_names[lines[index].chr_index][ chr_name_len ] = '\0';
            }
        /* Assuming this is a variable step numeric line */
        } else {
            int rv = sscanf( buffer, "%i\t%e\n", 
                             &(lines[index].position), &(lines[index].value)  );
            if( rv != 2 )
            {
                assert( rv == EOF );
                lines[index].fp = NULL;
            }

            return;
        }
    }
    
    /* this should never happen */
    assert( 0 );
    return;
}

extern void
aggregate_over_wiggles(
    FILE** wig_fps,
    int num_wigs,
    FILE* ofp,
    float agg_fn( const struct wig_line_info*, const int, const int  )
)
{
    char* rv;

    int curr_chr_index = -1;

    /* loop over num wigs */
    int i;

    /* BUG!!! HORRIBLE HACK */
    char** chr_names = malloc(sizeof(char*)*100);
    
    /* we implement this as a merge sort. 
       First, we read lines from each until we have a chr and position.
    */

    struct wig_line_info* lines = malloc( sizeof(struct wig_line_info)*num_wigs );
    
    /* get the track names */
    char buffer[500];
    for( i = 0; i < num_wigs; i++ )
    {
        lines[i].fp = wig_fps[i];
        lines[i].chr_index = -1;
        rv = fgets( buffer, 500, lines[i].fp );
        if( rv == NULL )
        {
            perror( "FATAL    : Could not read line from wiggle file" );
            exit( -1 );
        }
    }
    /* make sure that the header is not too long */
    if( strlen( buffer ) > 495 )
    {
        fprintf( stderr, "ERROR     : wiggle header too long ( %zu )", strlen( buffer ) );
        assert( 0 );
        exit( -1 );
    }
    /* write the header to the output file */
    fprintf( ofp, "%s", buffer );
    
    /* initialize the data arrays */
    for( i = 0; i < num_wigs; i++ )
    {
        parse_next_line( lines, chr_names, i );
    }
    
    qsort( lines, num_wigs, sizeof( struct wig_line_info ), cmp_wig_line_info );
    while( NULL != lines[0].fp )
    {
        if( lines[0].chr_index > curr_chr_index )
        {    
            fprintf( ofp, "variableStep chrom=%s\n", chr_names[lines[0].chr_index] );
            curr_chr_index = lines[0].chr_index;
        }
        
        unsigned int position = lines[0].position;
        int chr_index = lines[0].chr_index;
        int ub = 0;
        i = 1;
        while( position == lines[i].position 
               && chr_index == lines[i].chr_index
               && NULL != lines[i].fp
               && i < num_wigs )
        {
            ub++;
            i++;
        }

        float value = agg_fn( lines, ub, num_wigs );
        if( value > FLT_EPSILON )
            fprintf( ofp, "%i\t%e\n", position, value  );
        
        for( i = 0; i <= ub; i++ )
            parse_next_line( lines, chr_names, i );
        
        qsort( lines, num_wigs, sizeof( struct wig_line_info ), cmp_wig_line_info );
    }
    
    return;
}

/* this isnt implemented yet - I have no need */
#if 0
extern void
load_wiggle_into_trace(
    /* the trace to init and store the data into */
    struct trace_t* trace,
    /* contains the chr names and lengths */
    struct genome_data* genome;

    /* fp's for each wiggle - each must contain exactly one track */
    FILE** ifs,
    /* the number of wiggle files */
    int num_wigs
)
{
    /* Not yet implemented - ( I dont need it for anything ) */
    assert( false );
}
#endif

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


void
write_marginal_mapped_reads_to_stranded_wiggles( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    FILE* fwd_wfp, FILE* bkwd_wfp )
{
    unsigned int read_num = 1;
    
    const double filter_threshold = 1e-6;

    /* build and update the chr traces */
    struct trace_t* traces;
    init_traces( genome, &traces, 2 );
    
    struct mapped_read_t* rd;

    /* reset the file pointers in the mapped reads db */
    rewind_mapped_reads_db( rdb );

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &rd  ) )
    {
        if( read_num%1000000 == 0 )
            fprintf( stderr, "NOTICE       : Added %u reads into trace\n", read_num );
        
        set_read_fl_probs( rd, rdb->fl_dist );
        reset_read_cond_probs( rd );
        
        /* Update the trace from this mapping */
        unsigned int j;
        for( j = 0; j < rd->num_mappings; j++ ) {
            naive_update_trace_expectation_from_location( 
                traces, rd->locations + j );
        }
        
        free_mapped_read( rd );
        
        read_num++;
    }
    
    int i;
    unsigned int j;
    /* fwd stranded reads */
    fprintf( fwd_wfp, "track type=wiggle_0 name=%s\n", "fwd_strnd" );
    for( i = 0; i < traces->num_chrs; i++ )
    {
        fprintf( fwd_wfp, "variableStep chrom=%s\n", genome->chr_names[i] );
        for( j = 0; j < traces->trace_lengths[i]; j++ )
            if( traces->traces[0][i][j] >= filter_threshold )
                fprintf( fwd_wfp, "%i\t%e\n", j+1, traces->traces[0][i][j] );
    }

    /* rev stranded reads ( technically, reads whereas the rev comp of the first
       read mapped ) */
    fprintf( bkwd_wfp, "track type=wiggle_0 name=%s\n", "bkwd_strnd" );
    for( i = 0; i < traces->num_chrs; i++ )
    {
        fprintf( bkwd_wfp, "variableStep chrom=%s\n", genome->chr_names[i] );
        for( j = 0; j < traces->trace_lengths[i]; j++ )
            if( traces->traces[1][i][j] >= filter_threshold )
                fprintf( bkwd_wfp, "%i\t%e\n", j+1, -traces->traces[1][i][j] );
    }
    
    close_traces( traces );
}

