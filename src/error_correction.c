#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include "read.h"
#include "quality.h"
#include "genome.h"
#include "error_correction.h"
#include "find_candidate_mappings.h"
#include "statmap.h"

/*
 *  Code to estiamte error frequencies. Calls R.
 *
 */
SEXP*
init_double_vector( int data_len, double* input_data )
{
    assert( data_len >= 0 );
    SEXP* output_data = malloc(sizeof(SEXP));

    PROTECT(*output_data = allocVector(REALSXP, data_len));
    double* r_data_obj = REAL( *output_data );
    
    int i;
    for( i = 0; i < data_len; i++ )
    {
        r_data_obj[i] = input_data[i];
    }
    
    return output_data;
}

SEXP*
init_int_vector( int data_len, int* input_data )
{
    assert( data_len >= 0 );
    
    SEXP* output_data = malloc(sizeof(SEXP));

    PROTECT(*output_data = allocVector(INTSXP, data_len));
    int* r_data_obj = INTEGER( *output_data );
    
    int i;
    for( i = 0; i < data_len; i++ )
    {
        r_data_obj[i] = input_data[i];
    }
    
    return output_data;
}

void
build_vectors_from_error_data( 
    struct error_data_t* data,
    SEXP** r_poss, SEXP** r_qual_scores,
    SEXP** r_mm_cnts, SEXP** r_cnts
)
{
    double* poss;
    double* qual_scores;
    int* mm_cnts;
    int* cnts;
    
    int num_qual_scores = data->max_qual_score + 1;
    int num_read_pos = data->max_read_length + 1;
    
    int flat_vector_len = num_qual_scores*num_read_pos;
    
    poss = calloc( flat_vector_len, sizeof(double) );
    qual_scores = calloc( flat_vector_len, sizeof(double) );
    mm_cnts = calloc( flat_vector_len, sizeof(int) );
    cnts = calloc( flat_vector_len, sizeof(int) );

    int record_index;
    int pos, qual_score, flat_vec_pos;
    for( record_index = 0; record_index < data->num_records; record_index++ )
    {
        flat_vec_pos = 0;
        struct error_data_record_t* record = data->records[record_index];
        for( pos = 0; pos < num_read_pos; pos++ )
        {
            for( qual_score = 0;
                 qual_score < num_qual_scores; 
                 qual_score++ )
            {
                poss[ flat_vec_pos ] = pos;
                qual_scores[ flat_vec_pos ] = qual_score;
                mm_cnts[ flat_vec_pos ] += 
                    record->base_type_mismatch_cnts[qual_score][pos];
                cnts[ flat_vec_pos ] +=
                    record->base_type_cnts[qual_score][pos];
            
                flat_vec_pos += 1;
            }
        }
    }

    *r_poss = init_double_vector( flat_vector_len, poss );
    *r_qual_scores = init_double_vector( flat_vector_len, qual_scores );
    *r_mm_cnts = init_int_vector( flat_vector_len, mm_cnts );
    *r_cnts = init_int_vector( flat_vector_len, cnts );
    
    free(poss);
    free(qual_scores);
    free(mm_cnts);
    free(cnts);
    
    return;
}


void
init_freqs_array_from_error_data( struct freqs_array** est_freqs, 
                                  struct error_data_t* error_data )
{
    // TODO - move this into a function
    *est_freqs = malloc( sizeof(struct freqs_array) );
    (*est_freqs)->max_qual_score = error_data->max_qual_score;
    (*est_freqs)->max_position = error_data->max_read_length;
    (*est_freqs)->freqs = calloc( error_data->max_qual_score+1, sizeof(double*) );
    
    int i;
    for( i = 0; i <= error_data->max_qual_score; i++ )
    {
        (*est_freqs)->freqs[i] = calloc( error_data->max_read_length+1, sizeof(double) );
    }
    
    return;
}

void
free_freqs_array( struct freqs_array* est_freqs )
{
    int i;
    for( i = 0; i <= est_freqs->max_qual_score; i++ )
    {
        free( est_freqs->freqs[i] );
    }
    
    free( est_freqs->freqs );
    free( est_freqs );
    
    return;
}

void
update_freqs_array( struct freqs_array* est_freqs, 
                    struct error_data_t* error_data, 
                    double* flat_freqs  )
{
    /* update the array with the estiamted freqs */
    int pos, qual_score, flat_vec_pos;
    flat_vec_pos = 0;
    for( pos = 0; pos <= error_data->max_read_length; pos++ )
    {
        for( qual_score = 0;
             qual_score <= error_data->max_qual_score; 
             qual_score++ )
        {
            est_freqs->freqs[qual_score][pos] = flat_freqs[flat_vec_pos];
            flat_vec_pos += 1;
        }
    }
    
    return;
}

void
predict_freqs( struct error_data_t* data, int record_index, 
               struct freqs_array* predicted_freqs )
{
    fprintf( stderr, 
             "DEBUG       :  Predicting the error models for readkeys %i - %i.\n",
             data->records[record_index]->min_readkey, 
             data->records[record_index]->max_readkey );
    
    SEXP *mm_cnts = NULL;
    SEXP *cnts = NULL; 
    SEXP *poss = NULL;
    SEXP *qual_scores = NULL;

    build_vectors_from_error_data( 
        data, &poss, &qual_scores, &mm_cnts, &cnts
    );
    
    SEXP call;
    char plot_str[500];
    sprintf( plot_str, "error_dist_BS%i_%i_%i", 
             1,
             data->records[record_index]->min_readkey, 
             data->records[record_index]->max_readkey );
    
    call = lang6( install("predict_freqs"), 
                  *mm_cnts, *cnts, *poss, *qual_scores, 
                  mkString(plot_str) );
    
    SEXP res = eval( call, R_GlobalEnv );
    double *flat_freqs = REAL( res );
    
    update_freqs_array( predicted_freqs, data, flat_freqs );
    
    goto cleanup;
    
cleanup:

    free( mm_cnts );
    free( cnts );
    free( poss );
    free( qual_scores );
    
    return;
}


/*
 *
 *  END code to estimate error frequencies 
 *
 */


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
update_estimated_error_model_from_error_data( 
    struct error_model_t* error_model,
    struct error_data_t* data
)
{    
    struct freqs_array* pred_freqs = error_model->data;

    if( pred_freqs == NULL )
    {
        init_freqs_array_from_error_data( &pred_freqs, data );    
    }
    
    predict_freqs( data, data->num_records-1, pred_freqs );
    error_model->data = pred_freqs;

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
        update_estimated_error_model_from_error_data( error_model, data );
    } else {
        fprintf( stderr, "FATAL:        Unrecognized error type '%i'", 
                 error_model->error_model_type );
        exit(1);
    }
    
    return;
};

void free_error_model( struct error_model_t* error_model )
{
    if( error_model->error_model_type == ESTIMATED )
    {
        free_freqs_array( error_model->data );
    }

    free( error_model );
    return;
}

double
calc_min_match_penalty( struct penalty_t* penalties, int penalties_len, 
                        float exp_miss_frac )
{
    /* Unused - assert added to prevent compiler warning (for now) */
    assert( exp_miss_frac > 0 );

    double mean = 0;
    double var = 0;
    
    int i;
    for( i = 0; i < penalties_len; i++ )
    {
        /* we take [0][1] because this is a guaranteed mismatch,
           but this is wrong when we move to using actual by base
           mutation rates. */
        double penalty = penalties[i].penalties[0][1];
        double mm_prb = pow( 10, penalty );
        mean += mm_prb*penalty;
        mean += (1-mm_prb)*log10((1-mm_prb));
        
        var += (1-mm_prb)*mm_prb*penalty*penalty;
    }
    
    return mean - 4*sqrt( var );
}

int
calc_effective_sequence_length( struct penalty_t* penalties, int penalties_len )
{
    double mean = 0;
    int i;
    for( i = 0; i < penalties_len; i++ )
    {
        /* we take [0][1] because this is a guaranteed mismatch,
           but this is wrong when we move to using actual by base
           mutation rates. */
        double mm_prb = pow( 10, penalties[i].penalties[0][1] );
        mean += mm_prb;
    }
    
    return penalties_len - (int)mean - 1;
}

enum bool
filter_penalty_array( 
    struct penalty_t* penalties, int penalties_len,
    double effective_genome_len)
{
    int effective_seq_len = calc_effective_sequence_length(
        penalties, penalties_len );
                
    if (pow(4, effective_seq_len)/2 <= effective_genome_len ) {
        /*
        fprintf( stderr, "Filtering Read: %i %i - %e %e\n", 
                 effective_seq_len, penalties_len, 
                 pow(4, effective_seq_len)/2, effective_genome_len );
        */
        return true;
    }
    
    return false;
}

enum bool
filter_read(
        struct read* r,
        struct mapping_params* mapping_params,
        struct genome_data* genome
    )
{
    if( MISMATCH == mapping_params->metaparams->error_model_type )
    {
        return false;
    }
    
    int effective_genome_len = calc_genome_len( genome );
    
    // loop over the subtemplates, counting hq basepairs
    int i;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        if( filter_penalty_array( 
                mapping_params->fwd_penalty_arrays[i]->array, 
                mapping_params->fwd_penalty_arrays[i]->length,
                effective_genome_len) ) {
            return true;
        }
    }
    
    return false;
}

enum bool
filter_indexable_subtemplates(
    struct indexable_subtemplates* ists,
    struct mapping_params* mapping_params,
    struct genome_data* genome
)
{
    if( MISMATCH == mapping_params->metaparams->error_model_type )
    {
        return false;
    }
    
    int effective_genome_len = calc_genome_len( genome );
    
    int i;
    for( i = 0; i < ists->length; i++ )
    {
        struct indexable_subtemplate* ist = ists->container + i;

        if( filter_penalty_array( 
                ist->fwd_penalties, ist->subseq_length, effective_genome_len ) )
        {
            return true;
        }
    }

    return false;
}

void
init_mapping_params_for_read(
        struct mapping_params** p,
        struct read* r,        
        struct mapping_metaparams* metaparams,
        struct error_model_t* error_model
    )
{
    *p = malloc( sizeof( struct mapping_params ));

    (*p)->metaparams = metaparams;
    
    /* build the penalty arrays */
    (*p)->num_penalty_arrays = r->num_subtemplates;
    
    (*p)->fwd_penalty_arrays = calloc( 
        sizeof(struct penalty_array_t*), (*p)->num_penalty_arrays );
    (*p)->rev_penalty_arrays = calloc( 
        sizeof(struct penalty_array_t*), (*p)->num_penalty_arrays );
    
    int i;
    for( i = 0; i < (*p)->num_penalty_arrays; i++ )
    {
        (*p)->fwd_penalty_arrays[i] = calloc( 
            sizeof(struct penalty_array_t), r->subtemplates[i].length );
        build_penalty_array( (*p)->fwd_penalty_arrays[i],
                             r->subtemplates[i].length, 
                             error_model, 
                             r->subtemplates[i].error_str );
        
        (*p)->rev_penalty_arrays[i] = calloc( 
            sizeof(struct penalty_array_t), r->subtemplates[i].length );        
        build_reverse_penalty_array( 
            (*p)->rev_penalty_arrays[i], 
            (*p)->fwd_penalty_arrays[i]
        );
    }

    /* calculate the total read length. This is just the sum of the read 
       subtemplate lengths */
    int total_read_len = 0;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        total_read_len += r->subtemplates[i].length;
    }
    
    /* now, calcualte the model parameters */
    if( metaparams->error_model_type == MISMATCH ) {
        int max_num_mm = -(int)(
            metaparams->error_model_params[0]*total_read_len)-1;
        int max_mm_spread = (int)(
            metaparams->error_model_params[1]*total_read_len)+1;
        
        (*p)->recheck_min_match_penalty = max_num_mm;
        (*p)->recheck_max_penalty_spread = max_mm_spread;
    } 
    /* if the error model is estiamted, then just pass the meta params
       through ( for now ) */
    else {
        assert( metaparams->error_model_type == ESTIMATED );
        (*p)->recheck_min_match_penalty = 0;
        int j;
        for( j = 0; j < (*p)->num_penalty_arrays; j++ ) 
        {
            (*p)->recheck_min_match_penalty += 
                calc_min_match_penalty( (*p)->fwd_penalty_arrays[j]->array,
                                        (*p)->fwd_penalty_arrays[j]->length,
                                        1 - metaparams->error_model_params[0] );
        }
        
         //metaparams->error_model_params[0];
        (*p)->recheck_max_penalty_spread = log10(
                1 - metaparams->error_model_params[0] );
    }

    if( _assay_type == CAGE )
    {
        /* FIXME - for now, increase the penalty so we can at least map perfect
         * reads with up to MAX_NUM_UNTEMPLATED_GS untemplated G's */
        (*p)->recheck_min_match_penalty +=
            (MAX_NUM_UNTEMPLATED_GS*UNTEMPLATED_G_MARGINAL_LOG_PRB);
    }
    
    return;
}

void
init_index_search_params(
        struct index_search_params** isp,
        struct indexable_subtemplates* ists,
        struct mapping_params* mapping_params )
{
    /* Allocate memory for an array of index_search_params, one for each index
     * probe */
    *isp = malloc( sizeof( struct index_search_params ) * ists->length );
    
    /* TODO for now, set the index search params equal to the recheck params */
    int i;
    for( i = 0; i < ists->length; i++ )
    {
        float min_match_penalty = 1;
        float max_penalty_spread = -1;
        int length = ists->container[i].subseq_length;        
        if( mapping_params->metaparams->error_model_type == MISMATCH ) {

            min_match_penalty = -(int)(
                mapping_params->metaparams->error_model_params[0]*length)-1;

            /* Let the mismatch spread be 1/2 the allowed mismatch rate (for
             * now */
            max_penalty_spread = (int)(
                mapping_params->metaparams->error_model_params[0]*0.5*length)+1;

        } else {
            assert( mapping_params->metaparams->error_model_type == ESTIMATED );
            float expected_map_rate = 
                mapping_params->metaparams->error_model_params[0];
            
            min_match_penalty 
                = calc_min_match_penalty( 
                    ists->container[i].fwd_penalties,
                    ists->container[i].subseq_length,
                    1 - expected_map_rate );
            
            max_penalty_spread = mapping_params->recheck_max_penalty_spread;
        }        
        
        (*isp)[i].min_match_penalty = min_match_penalty;
        (*isp)[i].max_penalty_spread = max_penalty_spread;
    }
    
    return;
}

void
free_mapping_params( struct mapping_params* p )
{
    if( p == NULL ) return;
    
    if( NULL != p->fwd_penalty_arrays )
    {
        int i = 0;
        for( i = 0; i < p->num_penalty_arrays; i++ )
        {
            free_penalty_array(p->fwd_penalty_arrays[i]);
            free_penalty_array(p->rev_penalty_arrays[i]);
            free(p->fwd_penalty_arrays[i]);
            free(p->rev_penalty_arrays[i]);
        }
    
        free(p->fwd_penalty_arrays );
        free(p->rev_penalty_arrays );
    }
    
    free( p );
}

/*******************************************************************************
 *
 *
 * Error data
 *
 *
 ******************************************************************************/

void
init_error_data( struct error_data_t** data )
{
    *data = calloc( sizeof(struct error_data_t), 1 );

    (*data)->num_records = 0;
    (*data)->records = NULL;
    
    (*data)->max_read_length = MAX_READ_LEN;
    (*data)->max_qual_score = MAX_QUAL_SCORE;
    
    /*
     * mutex to control concurrent access to the 'data' struct
     */
    int rc;
    pthread_mutexattr_t mta;
    rc = pthread_mutexattr_init(&mta);
    (*data)->mutex = malloc( sizeof(pthread_mutex_t) );
    rc = pthread_mutex_init( (*data)->mutex, &mta );

    assert( rc == 0 );
    
    return;
};

void 
free_error_data( struct error_data_t* data )
{
    /* don't free unallocated memory */
    if( data == NULL )
        return;
    
    int i;
    for( i = 0; i < data->num_records; i++ )
    {
        free_error_data_record( data->records[i] );
    }
    
    free( data->records );
    
    /* release the mutex */
    assert( data->mutex != NULL );

    pthread_mutex_destroy( data->mutex );
    free( data->mutex );
    
    free( data );
};

void
add_new_error_data_record( 
    struct error_data_t* data, int min_readkey, int max_readkey )
{
    assert( data->num_records > 0 || data->records == NULL );
    
    struct error_data_record_t* record;
    init_error_data_record( 
        &record, data->max_read_length, data->max_qual_score );
    record->min_readkey = min_readkey;
    record->max_readkey = max_readkey;
    
    pthread_mutex_lock( data->mutex );
    
    data->num_records += 1;
    data->records = realloc( 
        data->records, sizeof(struct error_data_record_t*)*data->num_records );
    assert( NULL != data->records );
    data->records[ data->num_records-1 ] = record;
    
    pthread_mutex_unlock( data->mutex );
    
    return;
}

void
merge_in_error_data_record( struct error_data_t* data, int record_index,
                            struct error_data_record_t* record )
{
    /* if the record index is < 0, use the last index */
    if( record_index < 0 ) {
        record_index = data->num_records - 1;
    }
    
    assert( record_index < data->num_records && record_index >= 0 );
    
    
    pthread_mutex_lock( data->mutex );
    sum_error_data_records( data->records[record_index], record );
    pthread_mutex_unlock( data->mutex );
    
    return;
}

void
find_length_and_qual_score_limits( struct error_data_t* data,
                                   int* min_qual_score, int* max_qual_score,
                                   int* max_read_length )
{
    *min_qual_score = INT_MAX;
    *max_qual_score = 0;
    *max_read_length = 0;
    
    int i, j, k;
    for( i=0; i < data->num_records; i++ )
    {
        struct error_data_record_t* record = data->records[i];
        for( j = 0; j <= record->max_qual_score; j++ )
        {
            for( k = 0; k <= record->max_read_length; k++ )
            {
                if( record->base_type_cnts[j][k] > 0 )
                {
                    *min_qual_score = MIN( *min_qual_score, j );
                    *max_qual_score = MAX( *min_qual_score, j );
                    *max_read_length = MAX( *max_read_length, k );
                }
            }
        }
    }
    
    return;
}

void log_error_data( FILE* ofp, struct error_data_t* data )
{
    /* first find the bounds of the actually observed read lengths and quality
       scores so that we don't print out a bunch of zeros */
    int min_qual_score, max_qual_score, max_read_len;
    find_length_and_qual_score_limits( 
        data, &min_qual_score, &max_qual_score, &max_read_len );

    /* print the header */
    fprintf( ofp, "nreads\tmin_readkey\tmax_readkey\tmax_readlen\tmin_qualscore\tmax_qualscore" );
    
    int i, j;
    for( i = MAX(0, min_qual_score); 
         i <= MIN( max_qual_score, data->max_qual_score); 
         i++ )
    {
        for( j = 1; j <= MIN( max_read_len, data->max_read_length ); j++ )
        {
            fprintf( ofp, "\tE%i-%i", i, j );
        }
    }

    for( i = MAX(0, min_qual_score); 
         i <= MIN( max_qual_score, data->max_qual_score); 
         i++ )
    {
        for( j = 1; j <= MIN( max_read_len, data->max_read_length ); j++ )
        {
            fprintf( ofp, "\tC%i-%i", i, j );
        }
    }
    fprintf( ofp, "\n" );
    
    for( i=0; i < data->num_records; i++ )
    {
        /* sometimes there is trailing whitespace. Make sure these
           records aren't printed out */
        if( i == data->num_records - 1 
            && 0 == data->records[i]->num_unique_reads )
            break;
        
        fprintf_error_data_record( 
            ofp, data->records[i], 
            min_qual_score, max_qual_score, max_read_len );
    }
}

/*******************************************************************************
 *
 *
 * Error record data
 *
 *
 ******************************************************************************/

/*
 *
 * Max read length and max qual score are both allowed, so the arrays
 * are allocated to size max_qual_score+1 and max_read_len+1 ( although,
 * the 0 size read length is never used. 
 *
 */
void
init_error_data_record( struct error_data_record_t** data, 
                        int max_read_len, int max_qual_score )
{
    *data = malloc( sizeof(struct error_data_record_t) );
    
    (*data)->num_unique_reads = 0;
    (*data)->max_read_length = max_read_len;
    (*data)->max_qual_score = max_qual_score;

    (*data)->min_readkey = -1;
    (*data)->max_readkey = -1;
    
    (*data)->base_type_cnts = malloc( (1+max_qual_score)*sizeof(int*) );
    (*data)->base_type_mismatch_cnts = malloc( (1+max_qual_score)*sizeof(int*));

    int i;
    for( i = 0; i <= max_qual_score; i++) 
    {
        (*data)->base_type_cnts[i] 
            = malloc( (1 + max_read_len)*sizeof(int) );
        memset( (*data)->base_type_cnts[i], 0, 
                (1 + max_read_len)*sizeof(int)          );
        
        (*data)->base_type_mismatch_cnts[i] 
            = malloc((1 + max_read_len)*sizeof(int));
        memset( (*data)->base_type_mismatch_cnts[i], 0, 
                (1 + max_read_len)*sizeof(int)          );
    }
    
    return;
}

void
free_error_data_record( struct error_data_record_t* data )
{
    /* don't free unallocated memory */
    if( data == NULL )
        return;
    
    int i;
    for( i = 0; i < data->max_qual_score+1; i++ )
    {
        free( data->base_type_cnts[i] );
        free( data->base_type_mismatch_cnts[i] );
    }
    
    free( data->base_type_cnts );
    free( data->base_type_mismatch_cnts );
    
    free( data );
}

void
update_error_data_record(
    struct error_data_record_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int read_length
)
{
    data->num_unique_reads += 1;
    
    int i;    
    for( i = 0; i < read_length; i++ )
    {
        if( toupper(read[i]) != toupper(genome_seq[i])  )
        {
            // Add 1 because read positions are 1 based
            unsigned char error_char = (unsigned char) error_str[i];
            data->base_type_mismatch_cnts[error_char][i+1] += 1;
        }
        
        data->base_type_cnts[(unsigned char) error_str[i]][i+1] += 1;
    }
    
    return;
}

void
sum_error_data_records(
    struct error_data_record_t* dest,
    struct error_data_record_t* src
)
{
    int i, j;
    
    assert( dest->max_read_length == src->max_read_length );
    assert( dest->max_qual_score == src->max_qual_score );
    
    // add position_mismatch_cnts
    for( i = 0; i <= dest->max_qual_score; i++ ) 
    {
        for( j = 1; j <= dest->max_read_length; j++ ) 
        {
            dest->base_type_cnts[i][j] += src->base_type_cnts[i][j];
            dest->base_type_mismatch_cnts[i][j]
                    += src->base_type_mismatch_cnts[i][j];
        }
    }
    
    // add num_unique_reads
    dest->num_unique_reads += src->num_unique_reads;
    
    return;
}

void
fprintf_error_data_record( 
    FILE* stream, struct error_data_record_t* data, 
    int min_qual_score, int max_qual_score, int max_read_length )
{
    fprintf( stream, "%i", data->num_unique_reads );

    fprintf( stream, "\t%i", data->min_readkey );
    fprintf( stream, "\t\t%i", data->max_readkey );
    
    fprintf( stream, "\t\t%i", MIN( max_read_length, data->max_read_length) );
    fprintf( stream, "\t\t%i", MAX( 0, min_qual_score) );
    fprintf( stream, "\t\t%i\t", 
             MIN( max_qual_score, data->max_qual_score) );
    
    int i, j;
    for( i = MAX(0, min_qual_score); 
         i <= MIN( max_qual_score, data->max_qual_score); 
         i++ )
    {
        for( j = 1; j <= MIN( max_read_length, data->max_read_length ); j++ )
        {
            fprintf( stream, "\t%i", data->base_type_mismatch_cnts[i][j] );
        }
    }

    for( i = MAX(0, min_qual_score); 
         i <= MIN( max_qual_score, data->max_qual_score); 
         i++ )
    {
        for( j = 1; j <= MIN( max_read_length, data->max_read_length ); j++ )
        {
            fprintf( stream, "\t%i", data->base_type_cnts[i][j] );
        }
    }
    
    
    fprintf( stream, "\n" );
}
