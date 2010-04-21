#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "fragment_length.h"
#include "config.h"

void
init_fl_dist( struct fragment_length_dist_t** fl_dist, int min_fl, int max_fl )
{
    *fl_dist = malloc(sizeof(struct fragment_length_dist_t));
    (*fl_dist)->min_fl = min_fl;
    (*fl_dist)->max_fl = max_fl;
    
    (*fl_dist)->density = calloc( max_fl - min_fl + 1, sizeof(float)  );
    if( NULL == (*fl_dist)->density )
    {
        fprintf( stderr, "FATAL      : Failed to allocate memory for the fragment length dist." );
    }
    
    return;
}

void
free_fl_dist( struct fragment_length_dist_t* fl_dist )
{
    free( fl_dist->density );
    free( fl_dist );
    
    return;
}

void
init_fl_dist_from_file( struct fragment_length_dist_t** fl_dist, FILE* fp )
{
    /* make sure we are at the beginning of the file */
    fseek( fp, 0, SEEK_SET );
    /* determine the max, and minimum fragment lengths */
    int max = 0;
    int min = INT_MAX;
    while( !feof( fp ) )
    {
        int curr_val;
        fscanf( fp, "%i\t%*f\n", &curr_val );
        max = MAX( max, curr_val );
        min = MIN( min, curr_val );
    }

    /* init the data structure, and populate it */
    init_fl_dist( fl_dist, min, max );
    fseek( fp, 0, SEEK_SET );
    while( !feof( fp ) )
    {
        int curr_val;
        float curr_density;
        fscanf( fp, "%i\t%f\n", &curr_val, &curr_density );
        (*fl_dist)->density[curr_val-min] = curr_density;
    }
  
    fseek( fp, 0, SEEK_SET );
    return;
}

void
print_fl_dist( struct fragment_length_dist_t* fl_dist )
{
    int i;
    for( i = 0; i <= fl_dist->max_fl - fl_dist->min_fl; i++ )
    {
        printf( "%i\t%f\n", i+fl_dist->min_fl, fl_dist->density[i] );
    }
    
    return;
}

}

