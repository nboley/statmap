/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "trace.h"

#include "wiggle.h"


/* 
 *  This makes tons of assumptions about the 
 *  wiggle file format, order, that wont hold in 
 *  general. IE, this is a poor general solution,
 *  but it is fine to be used on files generated 
 *  by statmap 
 */
struct wig_line_info {
    FILE* fp;
    int chr_index;
    unsigned int position;
    float value;
};

extern void
aggregate_over_wiggles(
    FILE** wig_fps,
    int num_wigs,
    FILE* ofp
)
{
    /* loop over num wigs */
    int i;

    char** chr_names = malloc(sizeof(char*)*num_wigs);
    
    /* we implement this as a merge sort. 
       First, we read lines from each until we have a chr and position.
    */

    struct wig_line_info* lines = malloc( sizeof(struct wig_line_info)*num_wigs );
    
    /* get the track names */
    char buffer[500];
    for( i = 0; i < num_wigs; i++ )
    {
        lines[i].fp = wig_fps[i];
        fgets( buffer, 500, lines[i].fp );
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
        fgets( buffer, 500, &(lines[i].fp[i]) );
        if( buffer[0] == 'f' )
        {
            fprintf(stderr, "FATAL    : Wiggle parser does not support fixed step lines\n");
            assert( 0 );
            exit( -1 );
        } else if ( buffer[0] == 'w') {
            /* increment the chr index */
            lines[i].chr_index += 1;
            /* get the chr name */
            char* chr_name = strstr( buffer, "chrom=" );
            chr_name += 6;
            /* set the chr name */
            if( chr_names[i] != NULL )
            {
                if( 0 != strncmp( chr_name, chr_names[i], strlen( chr_names[i] ) ) )
                {
                    fprintf( stderr, "ERROR     : Chr names ( new: %.*s and old: %s ) are out of sync", 
                             (int)strlen(chr_names[i]), chr_name, chr_names[i] );
                    assert( 0 );
                }
            } else {
                chr_names[i] = calloc( 500, sizeof(char)  );
                sscanf( chr_name, "%s", chr_names[i] );
                chr_names = realloc( chr_names, sizeof(char)*strlen(chr_names[i]) );
            }
            
            /* move to the next line */
            i -= 1;
            continue;
        /* Assuming this is a variable step numeric line */
        } else {
            int rv = fscanf( lines[i].fp, "%i\t%e\n", 
                             &(lines[i].position), &(lines[i].value)  );
            if( rv != 2 )
            {
                assert( rv == EOF );
                lines[i].fp = NULL;
            }
            
            printf( "%i\t%e\n", lines[i].position, lines[i].value );
        }
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
