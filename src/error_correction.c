#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "error_correction.h"

void
init_error_data( 
    struct error_data_t** data,
    pthread_mutex_t* mutex
)
{
    *data = calloc(  sizeof(struct error_data_t), 1 );
    
    (*data)->num_unique_reads = 0;
    (*data)->max_read_length = 0;
    (*data)->position_mismatch_cnts = NULL;
    
    int j;
    for( j = 0; j < max_num_qual_scores; j++ )
    {
        /* add a little bias to the scores */
        (*data)->qual_score_cnts[j] = 0;
        (*data)->qual_score_mismatch_cnts[j] = 0;
    }

    /*
     * set pointer to a mutex to control concurrent access to this struct
     * set to NULL if no need for concurrent access
     */
    (*data)->mutex = mutex;
    
    return;
}

void 
free_error_data( struct error_data_t* data )
{
    if( NULL != data->position_mismatch_cnts )
        free( data->position_mismatch_cnts );
    free( data );
    
    return;
}

void
update_max_read_length(
    struct error_data_t* data,
    int new_max_read_length
)
{
    /* realloc memory */ data->position_mismatch_cnts = realloc(
        data->position_mismatch_cnts,
        sizeof(double) * new_max_read_length
    );
    /* initialize new positions to 0 */
    int i;
    for( i = data->max_read_length; i < new_max_read_length; i++ )
    {
        data->position_mismatch_cnts[ i ] = 0;
    }
    /* updata data's new_max_read_length */
    data->max_read_length = new_max_read_length;
}

/*
 * Merge src into dest with weight on src
 */
void
merge_error_data_structs_with_weight(
    struct error_data_t* dest,
    struct error_data_t* src,
    int weight
)
{
    int i;

    /* check max_read_length on dest and src, update dest if necessary */
    if( dest->max_read_length < src->max_read_length )
        update_max_read_length( dest, src->max_read_length );

    /* merge position_mismatch_cnts */
    for( i=0; i < dest->max_read_length; i++ )
    {
        dest->position_mismatch_cnts[ i ] += src->position_mismatch_cnts[ i ] * weight;
    }

    /* merge quality scores */
    for( i=0; i < max_num_qual_scores; i++ )
    {
        dest->qual_score_cnts[ i ] += src->qual_score_cnts[ i ] * weight;
        dest->qual_score_mismatch_cnts[ i ] += src->qual_score_mismatch_cnts[ i ] * weight;
    }

    /* update num_unique_reads in dest */
    dest->num_unique_reads += src->num_unique_reads;
}

void
add_error_data_structs(
    struct error_data_t* dest,
    struct error_data_t* src
)
{
    merge_error_data_structs_with_weight( dest, src, 1 );
}

/*
 * Weighted average of error data from dest and src, stored in dest
 */
void average_error_data_structs(
    struct error_data_t* dest,
    struct error_data_t* src,
    double weight
)
{
    assert( weight >= 0 && weight <= 1 );

    int i;
    /* check max_read_length on dest and src, update dest if necessary */
    if( dest->max_read_length < src->max_read_length )
        update_max_read_length( dest, src->max_read_length );

    /* average position_mismatch_cnts */
    for( i = 0; i < dest->max_read_length; i++ )
    {
        dest->position_mismatch_cnts[ i ] =
            dest->position_mismatch_cnts[ i ] * weight
                +
            src->position_mismatch_cnts[ i ] * ( 1 - weight );
    }

    /* average quality scores */
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        dest->qual_score_cnts[ i ] =
            dest->qual_score_cnts[ i ] * weight
                +
            src->qual_score_cnts[ i ] * ( 1 - weight );

        dest->qual_score_mismatch_cnts[ i ] =
            dest->qual_score_mismatch_cnts[ i ] * weight
                +
            src->qual_score_mismatch_cnts[ i ] * ( 1 - weight );
    }

    /* update num_unique_reads in dest */
    dest->num_unique_reads += src->num_unique_reads;
}

void
sync_global_with_local_error_data(
    struct error_data_t* global,
    struct error_data_t* local
)
{
    assert( global->mutex != NULL );

    /* threadsafe update of global error data */
    pthread_mutex_lock( global->mutex );

    average_error_data_structs( global, local, ERROR_WEIGHT );

    /* Write error data out to file */
    /*
     * Do it here so we can take advantage of the global error data's mutex
     * to avoid contention in writing to the log file
     */
    FILE* error_stats_log = fopen( ERROR_STATS_LOG, "a" );
    fprint_error_data( error_stats_log, global );
    fclose( error_stats_log );

    pthread_mutex_unlock( global->mutex );

    /* reset local error data */
    clear_error_data( local );
}

void
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int length
)
{
    int i;
    data->num_unique_reads += 1;

    /*
     * check length against max_read_length; if it is longer,
     * reallocate memory, initialize new memory to 0, and update max_read_length
     */
    if( data->max_read_length < length )
        update_max_read_length( data, length );
    
    for( i = 0; i < length; i++ )
    {
        if( toupper(read[i]) != toupper(genome_seq[i])  )
        {
            data->qual_score_mismatch_cnts[ (unsigned char) error_str[i] ] += 1;
            data->position_mismatch_cnts[ i ] += 1;
        }
        
        data->qual_score_cnts[ (unsigned char) error_str[i] ] += 1;
    }
    
    return;
}

void
clear_error_data( struct error_data_t* data )
{
    int j;
    /* reset quality score counts */
    for( j = 0; j < max_num_qual_scores; j++ )
    {
        data->qual_score_cnts[j] = 0;
        data->qual_score_mismatch_cnts[j] = 0;
    }
    /* reset position mismatch counts - free allocated memory and start over */
    free( data->position_mismatch_cnts );
    data->position_mismatch_cnts = NULL;
    
    data->num_unique_reads = 0;
    data->max_read_length = 0;
}

void
fprintf_error_data( FILE* stream, struct error_data_t* data )
{
    fprintf( stream, "Num Unique Reads:\t%i\n", data->num_unique_reads );
    fprintf( stream, "Max Read Length:\t%i\n", data->max_read_length );

    fprintf( stream, "Loc Error Rates:\n" );

    int i;
    for( i = 0; i < data->max_read_length; i++ )
    {
        fprintf( stream, "%i\t%e\n", i+1, position_mismatch_cnts[ i ] );
    }
    
    fprintf( stream, "Qual Score Error Rates:\n" );
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        float qs_error_rate;
        if( data->qual_score_cnts[ i ] == 0 ) // avoid NaNs
            qs_error_rate = 0;
        else
            qs_error_rate = data->qual_score_mismatch_cnts[ i ]/data->qual_score_cnts[ i ];

        fprintf( stream, "%i\t%e\n", i+1, qs_error_rate );
    }
    
    fprintf( stream, "\n" );
    
    return;
}
