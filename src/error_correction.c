#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "config.h"
#include "error_correction.h"

void
init_error_model( 
    struct error_model_t** error_model,
    enum error_model_type_t error_model_type
)
{
    *error_model = malloc( sizeof(struct error_model_t) );
    (*error_model)->error_model_type = error_model_type;
    (*error_model)->data = NULL;
    return;
}

void
update_error_model_from_error_data( 
    struct error_model_t* error_model,
    struct error_data_t* data
)
{
    if( error_model->error_model_type == MISMATCH )
    {
        assert( error_model->data == NULL );
        return;
    } else if ( error_model->error_model_type == FASTQ_MODEL ) {
        // We havn't implemented this yet...
        assert( false );
        assert( error_model->data == NULL );
        return;
    } else if ( error_model->error_model_type == ESTIMATED ) {
        error_model->data = data;
    } else {
        fprintf( stderr, "FATAL:        Unrecognized error type '%i'", 
                 error_model->error_model_type );
        exit(1);
    }
    
    return;
};

void free_error_model( struct error_model_t* error_model )
{
    if ( error_model->error_model_type == ESTIMATED ) {
        free_error_data( (struct error_data_t*) (error_model->data) );
        error_model->data = NULL;
    } else {
        /* if we didn't free the error model, it had better be NULL */
        assert( error_model->data == NULL );
    }
    
    free( error_model );
    return;
}


void
init_error_data( struct error_data_t** data )
{
    *data = calloc( sizeof(struct error_data_t), 1 );
    
    (*data)->num_unique_reads = 0;
    (*data)->max_read_length = 0;
    (*data)->position_mismatch_cnts = NULL;
    
    int j;
    for( j = 0; j < max_num_qual_scores; j++ )
    {
        (*data)->qual_score_cnts[j] = 0;
        (*data)->qual_score_mismatch_cnts[j] = 0;
    }

    /*
     * set pointer to a mutex to control concurrent access to this struct
     * set to NULL if no need for concurrent access
     */
    int rc;
    pthread_mutexattr_t mta;
    rc = pthread_mutexattr_init(&mta);
    (*data)->mutex = malloc( sizeof(pthread_mutex_t) );
    rc = pthread_mutex_init( (*data)->mutex, &mta );
    
    return;
}

void 
free_error_data( struct error_data_t* data )
{
    /* don't free unallocated memory */
    if( data == NULL )
        return;

    if( NULL != data->position_mismatch_cnts )
        free( data->position_mismatch_cnts );
    
    if( data->mutex != NULL )
        free( data->mutex );
    
    free( data );
}

void
update_max_read_length(
    struct error_data_t* data,
    int new_max_read_length
)
{
    /* realloc memory */
    data->position_mismatch_cnts = realloc(
        data->position_mismatch_cnts,
        sizeof(double) * new_max_read_length
    );
    assert( data->position_mismatch_cnts != NULL );

    /* initialize new positions to 0 */
    int i;
    for( i = data->max_read_length; i < new_max_read_length; i++ )
    {
        data->position_mismatch_cnts[ i ] = 0;
    }
    /* updata data's new_max_read_length */
    data->max_read_length = new_max_read_length;
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
add_error_data(
    struct error_data_t* dest,
    struct error_data_t* src
)
{
    int i;

    /* Don't bother trying to add error_data structs that store information about 0 reads */
    if( src->num_unique_reads == 0 )
        return;

    /* check max_read_length on dest and src, update dest if necessary */
    if( dest->max_read_length < src->max_read_length )
        update_max_read_length( dest, src->max_read_length );

    // add position_mismatch_cnts
    for( i = 0; i < dest->max_read_length; i++ )
        dest->position_mismatch_cnts[i] += src->position_mismatch_cnts[i];

    // add qual_score_cnts and qual_score_mismatch_cnts
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        dest->qual_score_cnts[i] += src->qual_score_cnts[i];
        dest->qual_score_mismatch_cnts[i] += src->qual_score_mismatch_cnts[i];
    }

    // add num_unique_reads
    dest->num_unique_reads += src->num_unique_reads;
}

void average_error_data(
    struct error_data_t* data
)
{
    /* Avoid divide by zero if data is initialized, but no reads have been stored */
    if( data->max_read_length == 0 )
        return;

    int i;

    // average position_mismatch_cnts
    for( i = 0; i < data->max_read_length; i++ )
        data->position_mismatch_cnts[i] = data->position_mismatch_cnts[i] / data->max_read_length;

    // average qual_score_cnts and qual_score_mismatch_cnts
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        data->qual_score_cnts[i] = data->qual_score_cnts[i] / data->max_read_length;
        data->qual_score_mismatch_cnts[i] = data->qual_score_mismatch_cnts[i] / data->max_read_length;
    }
}

/*
 * Weighted average of error data from dest and src, stored in dest
 */
void weighted_average_error_data(
    struct error_data_t* dest,
    struct error_data_t* src,
    double weight
)
{
    assert( weight >= 0 && weight <= 1 );

    /* If there are no unique mappers, we can't say anything new about the error estimates */
    if( src->max_read_length == 0 )
        return; // return dest unchanged

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
    dest->num_unique_reads =
        dest->num_unique_reads * weight
            +
        src->num_unique_reads * ( 1 - weight );
}

void log_error_data( struct error_data_t* ed )
{
    FILE* error_stats_log = fopen( ERROR_STATS_LOG, "a" );
    assert( error_stats_log != NULL );
    fprintf_error_data( error_stats_log, ed );
    fclose( error_stats_log );
}

void
update_global_error_data(
    struct error_data_t* global,
    struct error_data_t* local
)
{
    weighted_average_error_data( global, local, ERROR_WEIGHT );

    /* reset local error data */
    clear_error_data( local );
}

void
fprintf_error_data( FILE* stream, struct error_data_t* data )
{
    fprintf( stream, "num_unique_reads:\t%i\n", data->num_unique_reads );
    fprintf( stream, "max_read_length:\t%i\n", data->max_read_length );

    int i;
    for( i = 0; i < data->max_read_length; i++ )
    {
        fprintf( stream, "%i\t%e\n", i, data->position_mismatch_cnts[ i ] );
    }
    
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        fprintf( stream, "%i\t%e\t%e\n", i, 
                 data->qual_score_mismatch_cnts[ i ],
                 data->qual_score_cnts[i] );
    }
}

void
load_next_error_data_t_from_log_fp( struct error_data_t** ed,
                                    FILE* fp )
{
    /* clear old error data */
    if( *ed != NULL )
        free_error_data( *ed );

    /* set ed to NULL if we're at EOF */
    if( feof(fp) ) {
        *ed = NULL;
        return;
    }

    /* otherwise, the fp should be on the first line of an error_data_t entry */
    int rv, num_unique_reads, max_read_length;
    rv = fscanf( fp, "num_unique_reads:\t%i\n", &num_unique_reads );
    assert( rv == 1 ); // make sure the fp was on the first line of a new entry
    fscanf( fp, "max_read_length:\t%i\n", &max_read_length );
    assert( rv == 1 );

    /* set up error_data_t */
    init_error_data( ed );
    (*ed)->num_unique_reads = num_unique_reads;

    // load position_mismatch_cnts
    int i;
    for( i = 0; i < max_read_length; i++ )
    {
        double pmc;
        rv = fscanf( fp, "%*i %le", &pmc );
        assert( rv == 1 );
        (*ed)->position_mismatch_cnts[i] = pmc;
    }

    // load qual_score_cnts and qual_score_mismatch_cnts
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        double qsc, qsmc;
        rv = fscanf( fp, "%*i %le %le", &qsmc, &qsc );
        assert( rv == 2 );
        (*ed)->qual_score_cnts[i] = qsc;
        (*ed)->qual_score_mismatch_cnts[i] = qsmc;
    }
}
