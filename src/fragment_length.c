#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "fragment_length.h"
#include "config.h"
#include "mapped_read.h"

#include "log.h"

void
init_fl_dist( struct fragment_length_dist_t** fl_dist, int min_fl, int max_fl )
{
    assert( max_fl >= min_fl );

    *fl_dist = malloc(sizeof(struct fragment_length_dist_t));

    (*fl_dist)->min_fl = min_fl;
    (*fl_dist)->max_fl = max_fl;

    (*fl_dist)->chipseq_bs_density = NULL;
    (*fl_dist)->rev_chipseq_bs_density = NULL;

    (*fl_dist)->chipseq_bs_density_start = NULL;
    (*fl_dist)->rev_chipseq_bs_density_start = NULL;

    (*fl_dist)->density = calloc( 1 + max_fl - min_fl + 1, sizeof(float)  );
    
    if( NULL == (*fl_dist)->density )
    {
        statmap_log( LOG_FATAL, "Failed to allocate memory for the fragment length dist." );
        assert( 0 );
        exit( -1 );
    }
    
    /* align the fl density */
    
    return;
}

void
free_fl_dist( struct fragment_length_dist_t** fl_dist )
{
    if( *fl_dist == NULL ) return;

    free( (*fl_dist)->density );
    (*fl_dist)->density = NULL;
    
    if( NULL != (*fl_dist)->chipseq_bs_density ) 
    {
        free( (*fl_dist)->chipseq_bs_density_start );
        (*fl_dist)->chipseq_bs_density = NULL;
        (*fl_dist)->chipseq_bs_density_start = NULL;
        free( (*fl_dist)->rev_chipseq_bs_density_start );
        (*fl_dist)->rev_chipseq_bs_density = NULL;
        (*fl_dist)->rev_chipseq_bs_density_start = NULL;
    }
        
    free( (*fl_dist) );
    (*fl_dist) = NULL;
    
    return;
}

void
init_fl_dist_from_file( struct fragment_length_dist_t** fl_dist, FILE* fp )
{
    int rv;
    
    /* make sure we are at the beginning of the file */
    fseek( fp, 0, SEEK_SET );
    /* determine the max, and minimum fragment lengths */
    int max = 0;
    int min = INT_MAX;
    while( !feof( fp ) )
    {
        int curr_val;
        rv = fscanf( fp, "%i\t%*f\n", &curr_val );
        assert( rv == 1 );
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
        rv = fscanf( fp, "%i\t%f\n", &curr_val, &curr_density );
        assert( rv == 2 );
        (*fl_dist)->density[curr_val-min] = curr_density;
    }
  
    fseek( fp, 0, SEEK_SET );
    return;
}

void
build_fl_dist_from_filename( struct fragment_length_dist_t** fl_dist, char* filename )
{
    FILE* fl_fp = fopen( filename, "r" );
    if( fl_fp == NULL ) {
        statmap_log( LOG_FATAL, "Failed to open fl_dist from filename %s", filename );
    }
    init_fl_dist_from_file( fl_dist, fl_fp );
    fclose( fl_fp );
}

int 
get_frag_len( mapped_read_t* rd )
{
    mapped_read_index* rd_index;
    alloc_and_init_mapped_read_index(rd_index, rd);

    /* If there are not mapped locations, ther is no frag len */
    if( rd_index->num_mappings == 0 )
        return -1;
    
    /* if the reads aren't paired, we can't infer the frag len */
    enum bool reads_are_paired 
        = mapped_read_location_is_paired( rd_index->mappings + 0 );
    if( !reads_are_paired )
        return -1;
    
    /* initialize the fraglen to the fraglen of the first mapped location */
    int frag_len = 1 +
        get_stop_from_mapped_read_location( rd_index->mappings + 0 ) -
        get_start_from_mapped_read_location( rd_index->mappings + 0 );
    
    MPD_RD_ID_T i;
    for( i = 1; i < rd_index->num_mappings; i++ )
    {
        mapped_read_location* current_loc = rd_index->mappings[i];

        /* reads should never be a mixture of paired and unpaired reads 
         XXX IS THIS ACTUALLY TRUE? WE COULD JUST CONTINUE... */
        assert( mapped_read_location_is_paired( current_loc) );
        
        int tmp_fl = 1 +
            get_stop_from_mapped_read_location( current_loc ) -
            get_start_from_mapped_read_location( current_loc );
        
        /* if all of the fraglens don't match, then skip this rd */
        if( frag_len != tmp_fl )
            return -1;
    }

    free_mapped_read_index( rd_index );

    return frag_len;
}

void
estimate_fl_dist_from_mapped_reads(  struct mapped_reads_db* rdb )
{
    /* TODO 
     * this next line has no reason to be here - I just made a 
     * mistake and was too lazy to cleanup. 
     */
    struct fragment_length_dist_t** fl_dist = &(rdb->fl_dist);

    /* initialize a temporary fl dist */
    int total_num_reads = 0;
    int max_fl = TEMP_FL_ARRAY_GF;
    int* temp_fls = calloc(sizeof(int), max_fl+1);
    if( temp_fls == NULL )
    {
        statmap_log( LOG_FATAL, "Failed to allocate in estimate fl dist." );
        exit( -1 );
    }
    
    /* loop through each read */
    mapped_read_t* rd;
    rewind_mapped_reads_db( rdb );

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &rd ) )
    {
        int fl = get_frag_len( rd );
        if( fl == -1 )
            continue;

        /* add this length to the fl dist */
        if( fl > max_fl )
        {
            /* allocate space for the new data */
            temp_fls = realloc( temp_fls, sizeof(int)*(fl+1) );
            /* zero out the new entries */
            memset( temp_fls + max_fl + 1, 0, sizeof(int)*(fl-max_fl) );
            max_fl = fl;
        }

        total_num_reads += 1;
        temp_fls[ fl ] += 1;
    }
    
    /* check that we have looked at every read */
    assert( rdb->current_read == rdb->num_mapped_reads );

    if( total_num_reads < MIN_NUM_READS )
    {
        statmap_log( LOG_ERROR, "Too few unique reads to estimate the fl dist" );
        return;
    }

    
    /* do a pass over the data to determine the maximum and minimum fls */
    /* 
     *  NOTE THAT WE
     *  trim off the top and bottom 1% of cnts - this will bias the estimate 
     *  a bit, but I dont know what else to do. Really long fragments good
     *  really bias the results 
     */

    int new_max_fl = 0;
    int new_min_fl = max_fl;
    int curr_cnt = 0;
    int truncated_cnt = 0;
    double mean = 0;
    int fl;
    for( fl = 0; fl <= max_fl; fl++ )
    {
        if( temp_fls[fl] == 0 )
            continue;

        if( curr_cnt > (int) (UPPER_FL_CUTOFF_QUANTILE*total_num_reads+0.5) )
            continue;
        /* We update current here to always allow for the heavy boundaries */
        curr_cnt += temp_fls[fl];
        if( curr_cnt < (int) (LOWER_FL_CUTOFF_QUANTILE*total_num_reads) )
            continue;
        
        mean += fl*temp_fls[fl];
        truncated_cnt += temp_fls[fl];
        
        new_min_fl = MIN( fl, new_min_fl );
        new_max_fl = MAX( fl, new_max_fl );
    }
    /* make the sum the mean. */
    mean = mean/truncated_cnt;
    
    /* 
     * If there are too few fragments, we assume a truncated normal 
     * density with sample mean and variance. It's not very good, but 
     * it's probably better than use what is, potentially, very biased
     * empirical estiamte. 
     */
    /* insist on some level of average fragment */
    if( (new_max_fl - new_min_fl + 1)*MIN_EMPIRICAL_COVERAGE
        > total_num_reads*(UPPER_FL_CUTOFF_QUANTILE-LOWER_FL_CUTOFF_QUANTILE) )
    {
        double emp_variance = 0;
        /* estimate the variance */
        for( fl = new_min_fl; fl <= new_max_fl; fl++ )
        {
            emp_variance += temp_fls[fl]*(fl-mean)*(fl-mean);
        }
        emp_variance /= truncated_cnt;
        double std_dev = sqrt( emp_variance );
        
        /* make the new dist a normal treuncated at 2 stds */
        new_min_fl = MAX( 0, mean - 2*std_dev );
        new_max_fl = mean + 2*std_dev;
        
        double trace_sum = 0;
        init_fl_dist( fl_dist, new_min_fl, new_max_fl );
        for( fl = new_min_fl; fl <= new_max_fl; fl++ )
        {
            /* add something proportional the normal pdf */
            (*fl_dist)->density[fl-new_min_fl] 
                = (1/std_dev)*exp( 
                    -(fl-mean)*(fl-mean)/(2*emp_variance)
                );
            trace_sum += (*fl_dist)->density[fl-new_min_fl];
        }
        /* renormalize to 1 */
        for( fl = new_min_fl; fl <= new_max_fl; fl++ )
        {
            (*fl_dist)->density[fl-new_min_fl] /= trace_sum;
        }
        
        goto cleanup;
    }
    /* If we have enough for the empirical estimate */
    else {
        /* TODO - maybe some kernel smoothing? */
        init_fl_dist( fl_dist, new_min_fl, new_max_fl );
        for( fl = new_min_fl; fl <= new_max_fl; fl++ )
        {
            (*fl_dist)->density[ fl - new_min_fl  ] 
                = (float)temp_fls[fl]/total_num_reads;
        }
        
        goto cleanup;
    }

cleanup:
    free( temp_fls );
}

void
build_chipseq_bs_density( struct fragment_length_dist_t* fl_dist )
{
    /* below, we align to the 16 bit boundary so that we can make use of
       the sse vector operations in iterative mapping.
    */

    fl_dist->chipseq_bs_density = calloc( fl_dist->max_fl + 4, sizeof(float) );
    fl_dist->chipseq_bs_density_start = fl_dist->chipseq_bs_density;
    /* align to the 16 byte boundary */
    int offset = ((size_t) fl_dist->chipseq_bs_density )%16;
    if( 0 != offset )
        fl_dist->chipseq_bs_density = ( float* ) (
            ((size_t) fl_dist->chipseq_bs_density ) + (16-offset) );

    fl_dist->rev_chipseq_bs_density = calloc( fl_dist->max_fl + 4, sizeof(float) );
    fl_dist->rev_chipseq_bs_density_start = fl_dist->rev_chipseq_bs_density;
    /* align to the 16 byte boundary */
    offset = ((size_t) fl_dist->rev_chipseq_bs_density )%16;
    if( 0 != offset )
        fl_dist->rev_chipseq_bs_density = ( float* ) (
            ((size_t) fl_dist->rev_chipseq_bs_density ) + (16-offset) );
    
    int i = 0, frag_len = 0;
    for( frag_len = fl_dist->min_fl; frag_len <= fl_dist->max_fl; frag_len++ )
    {
        const double amt = 1.0/(frag_len);
        const double frag_len_prb = fl_dist->density[frag_len - fl_dist->min_fl];
        for( i = 0; i < frag_len; i++ )
        {
            fl_dist->chipseq_bs_density[i] += frag_len_prb*amt;
            fl_dist->rev_chipseq_bs_density[frag_len - i - 1] 
                = fl_dist->chipseq_bs_density[i];
        }
    }

    double sum = 0;
    for( i = 0; i < fl_dist->max_fl; i++ )
    {
        sum += fl_dist->chipseq_bs_density[i];
    }
    
    assert( sum > 0.99 && sum < 1.01 );
    
    return;
}

float
get_fl_prb( struct fragment_length_dist_t* fl_dist, int fl )
{
    if( fl_dist == NULL )
        return 1;

    if( fl < fl_dist->min_fl || fl > fl_dist->max_fl )
        return 0;
    
    return fl_dist->density[ fl - fl_dist->min_fl ];
}

float
get_fl_log_prb( struct fragment_length_dist_t* fl_dist, int fl )
{
    float fl_prb = get_fl_prb( fl_dist, fl );

    if( fl_prb <= 0 )
    {
        /* avoid log errors */
        return -FLT_MAX;
    }

    return log10( fl_prb );
}

void
fprint_fl_dist( FILE* fp, struct fragment_length_dist_t* fl_dist )
{
    int i;
    for( i = 0; i <= fl_dist->max_fl - fl_dist->min_fl; i++ )
    {
        fprintf( fp, "%i\t%f\n", i+fl_dist->min_fl, fl_dist->density[i] );
    }
    
    return;
}

/* Returns the log probability of a randomly sampled fragment length */
float
sample_fl_dist( struct fragment_length_dist_t* fl_dist )
{
    int fl_dist_len = fl_dist->max_fl - fl_dist->min_fl;
    assert( fl_dist_len > 0 );

    int random_fl = (random() % fl_dist_len) + fl_dist->min_fl;
    return get_fl_log_prb(fl_dist, random_fl);
}
