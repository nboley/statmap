#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include "log.h"
#include "rawread.h"
#include "read.h"
#include "quality.h"
#include "genome.h"
#include "error_correction.h"
#include "find_candidate_mappings.h"
#include "statmap.h"
#include "util.h"
#include "fragment_length.h"

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
    statmap_log(LOG_INFO, "Predicting the error models for readkeys %i - %i.",
            data->records[record_index]->min_readkey,
            data->records[record_index]->max_readkey
        );
    
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
        statmap_log( LOG_FATAL, "Unrecognized error type '%i'", 
                error_model->error_model_type  );
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
get_mismatch_penalty( struct penalty_t* penalty_array )
{
    double mismatch_penalty = 0;

    /* We can safely assume the penalty for a match (in log space) is > any
    penalty for a mismatch. */
    int i;
    for( i = 0; i < 4; i++ )
    {
        mismatch_penalty = MIN( penalty_array->penalties[i], mismatch_penalty );
    }

    return mismatch_penalty;
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
        /* Use the penalty of a mismatch, which we get by taking the min over
         * the four penalty values (since penalties are negative, this will be
         * the "highest" penalty score. This is wrong when we move to using
         * actual by base mutation rates. */
        double penalty = get_mismatch_penalty( penalties + i );
        double mm_prb = pow( 10, penalty );
        mean += mm_prb*penalty;
        mean += (1-mm_prb)*log10((1-mm_prb));
        
        var += (1-mm_prb)*mm_prb*penalty*penalty;
    }
    
    return mean - 4*sqrt( var );
}

/* Separate penalty mean calculation from below so we can log independently of
 * the underlying program logic */
void
log_penalties_mean( struct penalty_t* penalties, int penalties_len )
{
    double mean = 0;
    int i;
    for( i = 0; i < penalties_len; i++ ) {
        double mm_prb = pow( 10, get_mismatch_penalty( penalties + i ) );
        mean += mm_prb;
    }

    statmap_log(LOG_DEBUG, "penalties_mean %f", mean);
}

int
calc_effective_sequence_length( struct penalty_t* penalties, int penalties_len )
{
    double mean = 0;
    int i;
    for( i = 0; i < penalties_len; i++ ) {
        double mm_prb = pow( 10, get_mismatch_penalty( penalties + i ) );
        mean += mm_prb;
    }

    return penalties_len - (int)mean - 1;
}

enum bool
filter_penalty_array( 
        struct penalty_t* penalties, int penalties_len,
        double effective_genome_len
    )
{
    int effective_seq_len = calc_effective_sequence_length(
        penalties, penalties_len );

    if (pow(4, effective_seq_len)/2 <= effective_genome_len ) {
        /*
        statmap_log( LOG_DEBUG, "Filtering Read: %i %i - %e %e",
                effective_seq_len, penalties_len,
                pow(4, effective_seq_len)/2, effective_genome_len
            );
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
        /* WARNING this logging strategy will only work for unpaired reads */
        #if PROFILE_CANDIDATE_MAPPING
        log_penalties_mean( mapping_params->fwd_penalty_arrays[i]->array,
                mapping_params->fwd_penalty_arrays[i]->length );
        #endif

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

float
expected_value_of_rst_subsequence(
        struct penalty_array_t* rst_pens,
        int subseq_start,
        int subseq_length
    )
{
    float seq_ev = 0;
    int i, j;
    for(i = subseq_start; i < subseq_start + subseq_length; i++)
    {
        /* Expected value for this basepair */
        float bp_ev = 0;
        for(j = 0; j < 4; j++)
        {
            float bp_logprb = rst_pens->array[i].penalties[j];
            float bp_prb = pow(10, bp_logprb);
            bp_ev += bp_logprb * bp_prb;
        }

        seq_ev += bp_ev;
    }
    return seq_ev;
}

float
expected_value_of_rst(
        struct read_subtemplate* rst,
        struct penalty_array_t* rst_pens
    )
{
    return expected_value_of_rst_subsequence(rst_pens, 0, rst->length);
}

enum bool
check_sum_of_probabilities(
        struct penalty_t* bp_penalties
    )
{
    float sum = 0;
    int i;
    for( i = 0; i < 4; i++ )
    {
        sum += pow(10, bp_penalties->penalties[i] );
    }

    /* should be equal within a small range around 1.0 to account for rounding
     * errors */
    if( sum < (1.0 - ROUND_ERROR) || sum > (1.0 + ROUND_ERROR) ) {
        statmap_log( LOG_FATAL,
                "Sum of bp mutation probabilities (%f) is not within range [%f, %f]",
                sum, 1.0 - ROUND_ERROR, 1.0 + ROUND_ERROR );
    }

    return true;
}

float
sample_bp_penalties(
        struct penalty_t* bp_penalties
    )
{
    /* Check that the probabilities sum to 1 */
    assert( check_sum_of_probabilities( bp_penalties ) );

    /* Sample randomly from the bp mutation probabilities */
    /* Break if the sum is greater than a random number in the range [0,1] */
    double rand_val = ((double)rand() / (double)RAND_MAX);
    float sum = 0;
    int i;
    for( i = 0; i < 4; i++ )
    {
        sum += pow(10, bp_penalties->penalties[i]);
        if( sum > rand_val )
            break;
    }

    /* return the randomly sampled probability */
    return bp_penalties->penalties[i];
}

void
sort_samples_array(
        float** samples_array, // reference parameter (for realloc)
        int num_samples
    )
{
    /* resize the samples_array to match the number of samples it contains. */
    *samples_array = realloc( *samples_array, sizeof(float)*num_samples );
    assert( *samples_array != NULL );

    /* sort */
    qsort( *samples_array,
           num_samples,
           sizeof(float),
           (int(*)(const void*, const void*))cmp_floats );
}

int
compute_sampled_penalties_for_reads(
        struct rawread_db_t* rdb,
        struct error_model_t* error_model,
        struct fragment_length_dist_t* fl_dist,
        float quantile,
        struct sampled_penalties_t* sampled_penalties
    )
{
    int rv = 0; // communicate success/failure to caller

    struct rawread_db_state saved_state
        = save_rawread_db_state( rdb );

    /* Sample from the next set of READS_STAT_UPDATE_STEP_SIZE reads */
    readkey_t max_readkey = rdb->readkey + READS_STAT_UPDATE_STEP_SIZE;

    int num_samples
        = READS_STAT_UPDATE_STEP_SIZE * NUM_SAMPLES_FOR_MIN_PENALTY_COMP;

    float* read_samples = calloc( num_samples, sizeof(float) );
    float* read_subtemplate_samples = calloc( num_samples, sizeof(float) );

    struct read* r;
    int sample_index = 0;
    while( EOF != get_next_read_from_rawread_db(
                rdb, &r, max_readkey ) )
    {
        struct penalty_array_t* penalty_arrays = malloc(
                sizeof(struct penalty_array_t) * r->num_subtemplates );
        int st_i;
        for( st_i = 0; st_i < r->num_subtemplates; st_i++ )
        {
            build_penalty_array( penalty_arrays + st_i,
                                 r->subtemplates + st_i,
                                 error_model );
        }

        /* Randomly sample penalties */
        int i;
        for( i = 0; i < NUM_SAMPLES_FOR_MIN_PENALTY_COMP; i++ )
        {
            /* the full read sampled penalty is the sum of the sampled
               subtemplate penalties */
            float read_prb = 0;

            for( st_i = 0; st_i < r->num_subtemplates; st_i++ )
            {
                int bp_i;
                for( bp_i = 0; bp_i < penalty_arrays[st_i].length; bp_i++ )
                {
                    read_prb += sample_bp_penalties(
                            penalty_arrays[st_i].array + bp_i );
                }
            }
            assert( read_prb < 0 );

            read_subtemplate_samples[sample_index] = read_prb;

            /* if the read is a full fragment, also randomly sample from the
               fragment length distribution */
            if( r->prior.frag_type == FULL_GENOME_FRAGMENT ||
                r->prior.frag_type == FULL_TRANSCRIPTOME_FRAGMENT )
            {
                assert( fl_dist != NULL );
                read_prb += sample_fl_dist( fl_dist );
            }

            read_samples[sample_index++] = read_prb;
        }

        /* cleanup memory */
        for( st_i = 0; st_i < r->num_subtemplates; st_i++ ) {
            free_penalty_array( penalty_arrays + st_i );
        }
        free( penalty_arrays );
        free_read( r );
    }

    restore_rawread_db_state( rdb, saved_state );

    if( sample_index == 0 ) {
        /* we didn't process any reads - cleanup and return EOF */
        rv = EOF;
        goto cleanup;
    }

    /* resize the sample_prbs array so it's the exact size of the number of
     * sample probabilities (we could have processed fewer than
     * READS_STAT_UPDATE_STEP_SIZE reads, and don't want to mess up the qsort) */
    num_samples = sample_index;
    sort_samples_array( &read_samples, num_samples );
    sort_samples_array( &read_subtemplate_samples, num_samples );

    assert( quantile > 0.0 && quantile < 1.0 );

    /* Since these are log prbs, and qsort sorts in ascending order, the
     * indexes of the quantiles are reversed */
    int quantile_index = (int)(num_samples*(1-quantile));
    /* Add fudge factor to account for rounding error */
    sampled_penalties->read_penalty = read_samples[quantile_index] - ROUND_ERROR;
    sampled_penalties->read_subtemplate_penalty
        = read_subtemplate_samples[quantile_index] - ROUND_ERROR;

cleanup:
    /* cleanup memory */
    free( read_samples );
    free( read_subtemplate_samples );

    return rv;
}

struct mapping_params*
init_mapping_params_for_read(
        struct read* r,        
        struct mapping_metaparams* metaparams,
        struct error_model_t* error_model,
        struct sampled_penalties_t* sampled_penalties
    )
{
    struct mapping_params* p = malloc(sizeof(struct mapping_params));

    p->metaparams = metaparams;
    
    /* build the penalty arrays */
    p->num_penalty_arrays = r->num_subtemplates;
    
    p->fwd_penalty_arrays = calloc( 
        sizeof(struct penalty_array_t*), p->num_penalty_arrays );
    p->rev_penalty_arrays = calloc( 
        sizeof(struct penalty_array_t*), p->num_penalty_arrays );
    
    int i;
    for( i = 0; i < p->num_penalty_arrays; i++ )
    {
        p->fwd_penalty_arrays[i] = malloc( sizeof(struct penalty_array_t) );
        build_penalty_array( p->fwd_penalty_arrays[i], r->subtemplates + i, 
            error_model );
        
        p->rev_penalty_arrays[i] = malloc( sizeof(struct penalty_array_t) );
        build_reverse_penalty_array( p->rev_penalty_arrays[i], 
            r->subtemplates + i, error_model );
    }

    /* calculate the total read length. This is just the sum of the read 
       subtemplate lengths */
    int total_read_len = 0;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        total_read_len += r->subtemplates[i].length;
    }
    p->total_read_length = total_read_len;

    /* compute the expected value over the read. We use this to determine the
       scaling factor to get the min match penalty for an index probe from the
       min match penalty for the entire read. */
    p->read_expected_value = 0;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        p->read_expected_value += expected_value_of_rst(r->subtemplates + i,
            p->fwd_penalty_arrays[i]);
    }

    p->sampled_penalties = sampled_penalties;
    
    /* now, calcualte the model parameters */
    if( metaparams->error_model_type == MISMATCH ) {
        int max_num_mm = -(int)(
            metaparams->error_model_params[0]*total_read_len)-1;
        int max_mm_spread = (int)(
            metaparams->error_model_params[1]*total_read_len)+1;
        
        p->recheck_min_match_penalty = max_num_mm;
        p->recheck_max_penalty_spread = max_mm_spread;
    } 
    /* if the error model is estiamted, then just pass the meta params
       through ( for now ) */
    else {
        assert( metaparams->error_model_type == ESTIMATED );

        p->recheck_min_match_penalty = sampled_penalties->read_penalty;
        p->recheck_max_penalty_spread
            = -log10(1 - metaparams->error_model_params[0]);
        assert( p->recheck_max_penalty_spread >= 0 );
        //p->recheck_max_penalty_spread = 1.3;

        if( r->prior.assay == CAGE )
        {
            /* FIXME - for now, increase the penalty so we can at least map perfect
             * reads with up to MAX_NUM_UNTEMPLATED_GS untemplated G's */
            p->recheck_min_match_penalty = MAX(
                p->recheck_min_match_penalty,
                MAX_NUM_UNTEMPLATED_GS*UNTEMPLATED_G_MARGINAL_LOG_PRB);
        }
    }

    return p;
}

struct index_search_params*
init_index_search_params(
        struct indexable_subtemplates* ists,
        struct mapping_params* mapping_params )
{
    /* Allocate an array of index_search_params structs, one for each index probe */
    struct index_search_params *isp
        = malloc(sizeof(struct index_search_params)*ists->length);

    int i;
    for( i = 0; i < ists->length; i++ )
    {
        float min_match_penalty = 1;
        float max_penalty_spread = -1;

        struct indexable_subtemplate *ist = ists->container + i;

        if( mapping_params->metaparams->error_model_type == MISMATCH )
        {
            /* The first metaparam is the expected rate of mapping (assuming 
            the read came from the genome) */
            float expected_map_rate
                = mapping_params->metaparams->error_model_params[0];

            min_match_penalty = -(int)(expected_map_rate*ist->subseq_length)-1;
            /* Let mismatch spread be 1/2 the allowed mismatch rate (for now) */
            max_penalty_spread = (int)(expected_map_rate*0.5*ist->subseq_length)+1;
        }
        else {
            assert( mapping_params->metaparams->error_model_type == ESTIMATED );

            float scaling_factor 
                = ist->expected_value / mapping_params->read_expected_value;
            /* TODO: read_subtemplate_penalty is a potentially confusing name - it is
               the sampled penalty of all the read subtempaltes, without the fl dist penalty */
            min_match_penalty 
                = scaling_factor * mapping_params->sampled_penalties->read_subtemplate_penalty;
            //min_match_penalty = mapping_params->recheck_min_match_penalty;

            #if PROFILE_CANDIDATE_MAPPING
            statmap_log(LOG_DEBUG, "ist_min_match_penalty %f", min_match_penalty);
            #endif
                
            max_penalty_spread = mapping_params->recheck_max_penalty_spread;
        }    

        /* Make sure the index search parameters have the correct signs */
        assert( min_match_penalty <= 0 );
        assert( max_penalty_spread >= 0 );
        
        isp[i].min_match_penalty = min_match_penalty;
        isp[i].max_penalty_spread = max_penalty_spread;
    }
    
    return isp;
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
