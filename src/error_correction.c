#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "error_correction.h"

#include "config.h"
#include "log.h"
#include "rawread.h"
#include "read.h"
#include "quality.h"
#include "genome.h"
#include "find_candidate_mappings.h"
#include "statmap.h"
#include "util.h"

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
    struct error_data_t* data, int record_index,
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
    
    int pos, qual_score, flat_vec_pos;
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
        (*est_freqs)->freqs[i] = calloc( 
            error_data->max_read_length+1, sizeof(double) );
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
        data, record_index, &poss, &qual_scores, &mm_cnts, &cnts
    );
    
    SEXP call;
    char plot_str[500];
    sprintf( plot_str, "error_dist Rd Keys:%i-%i Strand:%i Rd Pair:%i", 
             data->records[record_index]->min_readkey, 
             data->records[record_index]->max_readkey,
             data->records[record_index]->strand, 
             data->records[record_index]->read_subtemplate_index );
    
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
    struct freqs_array** pred_freqs = error_model->data;

    if( pred_freqs == NULL )
    {
        pred_freqs = calloc(sizeof(struct freqs_array*), data->num_records);
    }
    
    int i;
    for(i = 0; i < data->num_records; i++) {
        if( pred_freqs[i] == NULL ) {
            init_freqs_array_from_error_data( &(pred_freqs[i]), data );
        }
        predict_freqs( data, i, pred_freqs[i] );
    }
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
        int i;
        for( i = 0; i < MAX_NUM_RD_SUBTEMPLATES*2; i++ ) 
        {
            free_freqs_array( ((struct freqs_array**)error_model->data)[i] );
        }
        free( (struct freqs_array**)error_model->data );        
    }

    free( error_model );
    return;
}

double
get_mismatch_probability( struct penalty_t* penalty_array )
{
    double match_penalty = -DBL_MAX;

    /* We can safely assume the penalty for a match (in log space) is > any
    penalty for a mismatch. */
    int i;
    for( i = 0; i < 4; i++ )
    {
        match_penalty = MAX( penalty_array->penalties[i], match_penalty );
    }

    return 1 - pow(10, match_penalty);
}


double*
calc_penalty_dist_moments( 
    struct penalty_array_t* penalty_array, int num_moments)
{
    /* Unused - assert added to prevent compiler warning (for now) */
    assert( num_moments == 2 );
    
    double *moments = calloc(sizeof(double), num_moments);
    double mean = 0;
    double var = 0;
    
    int i;
    for( i = 0; i < penalty_array->length; i++ )
    {
        /* Use the penalty of a mismatch, which we get by taking the min over
         * the four penalty values (since penalties are negative, this will be
         * the "highest" penalty score. This is wrong when we move to using
         * actual by base mutation rates. */
        double mm_prb = get_mismatch_probability( penalty_array->array  + i );
        double mm_penalty = log10(mm_prb) - LOG10_3;
        double match_penalty = log10(1-mm_prb);
        mean += mm_prb*mm_penalty;
        mean += (1-mm_prb)*match_penalty;
        
        var += mm_prb*(1-mm_prb)*(mm_penalty - match_penalty)*(mm_penalty - match_penalty);
    }

    moments[0] = mean;
    moments[1] = var;
    
    return moments;
}

/* Separate penalty mean calculation from below so we can log independently of
 * the underlying program logic */
void
log_penalties_mean( struct penalty_t* penalties, int penalties_len )
{
    double mean = 0;
    int i;
    for( i = 0; i < penalties_len; i++ ) {
        double mm_prb = get_mismatch_probability( penalties + i );
        mean += mm_prb;
    }

    #ifdef PROFILE_CANDIDATE_MAPPING
    statmap_log(LOG_DEBUG, "penalties_mean %f", mean);
    #endif
}

int
calc_effective_sequence_length( struct penalty_array_t* penalties )
{
    double mean = 0;
    int i;
    for( i = 0; i < penalties->length; i++ ) {
        double mm_prb = get_mismatch_probability( penalties->array + i );
        mean += mm_prb;
    }

    return penalties->length - (int)mean - 1;
}

enum bool
filter_penalty_array( struct penalty_array_t* penalty_array, 
                      double effective_genome_len )
{
    int effective_seq_len = calc_effective_sequence_length( penalty_array );
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
        log_penalties_mean( mapping_params->fwd_penalty_arrays[i]->array,
                mapping_params->fwd_penalty_arrays[i]->length );

        if( filter_penalty_array( mapping_params->fwd_penalty_arrays[i], effective_genome_len) )
        {
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

        if( filter_penalty_array(&(ist->fwd_penalty_array), effective_genome_len) )
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

float
compute_min_match_penalty_for_reads(
        struct rawread_db_t* rdb,
        struct error_model_t* error_model,
        int num_reads,
        float quantile
    )
{
    struct rawread_db_state saved_state
        = save_rawread_db_state( rdb );

    /* Sample from the next set of READS_STAT_UPDATE_STEP_SIZE reads */
    readkey_t max_readkey = rdb->readkey + num_reads;
    
    int num_sample_prbs = 
        NUM_READ_SAMPLES_FOR_MIN_PENALTY_COMP
            *NUM_BASE_SAMPLES_FOR_MIN_PENALTY_COMP;
    float* sample_prbs = calloc( num_sample_prbs, sizeof(float) );
    int sample_prbs_index = 0;

    struct read* r;
    while( sample_prbs_index < num_sample_prbs
           && EOF != get_next_read_from_rawread_db(
               rdb, &r, max_readkey )
        )
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
        
        int sample_i;
        for( sample_i = 0;
             sample_i < NUM_BASE_SAMPLES_FOR_MIN_PENALTY_COMP;
             sample_i++ )
        {
            float sample_prb = 0;

            int bp_i;
            for( st_i = 0; st_i < r->num_subtemplates; st_i++ )
            {
                for( bp_i = 0; bp_i < penalty_arrays[st_i].length; bp_i++ )
                {
                    sample_prb += sample_bp_penalties(
                            penalty_arrays[st_i].array + bp_i );
                }
            }

            assert( sample_prb < 0 );
            sample_prbs[sample_prbs_index++] = sample_prb;
            if( sample_prbs_index == num_sample_prbs )
                break;
        }

        /* cleanup memory */
        for( st_i = 0; st_i < r->num_subtemplates; st_i++ ) {
            free_penalty_array( penalty_arrays + st_i );
        }
        free( penalty_arrays );
        free_read( r );
    }

    restore_rawread_db_state( rdb, saved_state );

    float min_match_penalty = 0;

    if( sample_prbs_index == 0 ) {
        /* we didn't process any reads - cleanup and return, leaving
         * min_match_penalty set to 0 (invalid value) */
        goto cleanup;
    }

    /* resize the sample_prbs array so it's the exact size of the number of
     * sample probabilities (if there are fewer than
     * READS_STAT_UPDATE_STEP_SIZE reads in a block, then they won't match and
     * that will mess up the qsort. This happens naturally with small data sets
     * that have less than the step size reads, or in the last set from a large
     * data set.) */
    num_sample_prbs = sample_prbs_index;
    sample_prbs = realloc( sample_prbs, sizeof(float)*num_sample_prbs );

    /* sort the mismatch_prbs and take the specified quantile */
    qsort( sample_prbs, 
           num_sample_prbs,
           sizeof(float),
           (int(*)(const void*, const void*))cmp_floats );

    assert( quantile > 0.0 && quantile < 1.0 );
    /* Since these are log prbs, and qsort sorts in ascending order, the
     * indexes of the quantiles are reversed */
    min_match_penalty = sample_prbs[(int)(num_sample_prbs*(1-quantile))];

    /* Add fudge factor to account for rounding error */
    min_match_penalty -= ROUND_ERROR;

cleanup:
    /* cleanup memory */
    free( sample_prbs );

    return min_match_penalty;
}

struct mapping_params*
init_mapping_params_for_read(
        struct read* r,        
        struct mapping_metaparams* metaparams,
        struct error_model_t* error_model
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
    
    /* now, calculate the model parameters */
    if( metaparams->error_model_type == MISMATCH ) {
        int max_num_mm = -(int)(
            metaparams->error_model_params[0]*total_read_len)-1;
        int max_mm_spread = (int)(
            metaparams->error_model_params[1]*total_read_len)+1;
        
        p->recheck_min_match_penalty = max_num_mm;
        p->recheck_max_penalty_spread = max_mm_spread;
    } 
    else {
        assert( metaparams->error_model_type == ESTIMATED );
        
        double mean = 0; 
        double variance = 0; 
        
        for( i = 0; i < r->num_subtemplates; i++ )
        {
            double* fwd_moments = calc_penalty_dist_moments(
                p->fwd_penalty_arrays[i], 2);
            double* rev_moments = calc_penalty_dist_moments(
                p->rev_penalty_arrays[i], 2);
            
            mean += MIN(fwd_moments[0], rev_moments[0]);
            variance += MAX(fwd_moments[1], rev_moments[1]);
            free(fwd_moments);
            free(rev_moments);
        }

        double shape = mean*mean/variance;        
        double scale = -variance/mean;
        
        // Allow for one additional high quality mismatch, for the continuity error
        double min_match_prb = -qgamma(DEFAULT_ESTIMATED_ERROR_METAPARAMETER, 
                                       shape, scale, 1, 0) - 2.1;
        p->recheck_min_match_penalty = min_match_prb;
        p->recheck_max_penalty_spread = -log10(
            1 - DEFAULT_ESTIMATED_ERROR_METAPARAMETER);
        //printf("%e:%e:%e\n", min_match_prb, mean, variance);
        assert( p->recheck_max_penalty_spread >= 0 );
        
        if( _assay_type == CAGE )
        {
            /* FIXME - for now, increase the penalty so we can at least map perfect
             * reads with up to MAX_NUM_UNTEMPLATED_GS untemplated G's */
            p->recheck_min_match_penalty += 
                MAX_NUM_UNTEMPLATED_GS*UNTEMPLATED_G_MARGINAL_LOG_PRB;
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

        if( mapping_params->metaparams->error_model_type == MISMATCH ) {
            /* The first metaparam is the expected rate of mapping (assuming 
            the read came from the genome) */
            float max_mm_rate
                = mapping_params->metaparams->error_model_params[0];
            float max_mm_spread
                = mapping_params->metaparams->error_model_params[1];

            /* TODO
               make the multiple index probe correction. This should actually be
               qbinom( 0.5, 20, min_match_rate**n ), but this is actually nearly
               linear over reasonable probe lengths */
            double adj_p = exp(log(0.01)/ists->length);
            min_match_penalty = -ceil(
                qbinom(adj_p, ist->subseq_length, max_mm_rate, 0, 0));
            max_penalty_spread = ceil(
                qbinom(adj_p, ist->subseq_length, max_mm_spread, 0, 0));            
        } else {
            assert( mapping_params->metaparams->error_model_type == ESTIMATED );
            
            //float scaling_factor = ((double)ist->subseq_length )
            //    /mapping_params->total_read_length;
            double* fwd_moments = calc_penalty_dist_moments(
                &(ist->fwd_penalty_array), 2);
            double* rev_moments = calc_penalty_dist_moments(
                &(ist->rev_penalty_array), 2);
            double mean = MIN(fwd_moments[0], rev_moments[0]);
            double variance = MAX(fwd_moments[1], rev_moments[1]);
            free(fwd_moments);
            free(rev_moments);

            double shape = mean*mean/variance;        
            double scale = -variance/mean;

            // Adjust the p-value so that the mismatch rate is based 
            // upon every search missing
            double adj_p = exp(log(1-DEFAULT_ESTIMATED_ERROR_METAPARAMETER)/ists->length);
            min_match_penalty = -qgamma(adj_p, shape, scale, 0, 0) - 1e-3;
            
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

    (*data)->num_records = MAX_NUM_RD_SUBTEMPLATES*2;
    (*data)->records = calloc(
        sizeof(struct error_data_t*), (*data)->num_records);
    
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
    
    /* add a error data record for each strand, and for each possible 
       read template */
    int subtemplate_index, strand;
    for( subtemplate_index= 0; 
         subtemplate_index < MAX_NUM_RD_SUBTEMPLATES; 
         subtemplate_index++ )
    {
        for( strand = 1; strand < 3; strand++ ) 
        {
            init_error_data_record( 
                (*data)->records + 2*subtemplate_index + (strand-1),
                subtemplate_index, (enum STRAND)strand,
                (*data)->max_read_length, 
                (*data)->max_qual_score );
        }
    }
    
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

/*
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
*/

void
merge_in_error_data( struct error_data_t* data_to_update, 
                     struct error_data_t* data )
{
    int i;
    assert( data_to_update->num_records == data->num_records );
    pthread_mutex_lock( data_to_update->mutex );
    for( i = 0; i < data_to_update->num_records; i++ )
    {
        sum_error_data_records( 
            data_to_update->records[i], 
            data->records[i] );
    }
    pthread_mutex_unlock( data_to_update->mutex );
    
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
                        int read_subtemplate_index, enum STRAND strand,
                        int max_read_len, int max_qual_score )
{
    *data = malloc( sizeof(struct error_data_record_t) );
    
    (*data)->num_unique_reads = 0;
    
    (*data)->max_read_length = max_read_len;
    (*data)->max_qual_score = max_qual_score;

    (*data)->read_subtemplate_index = read_subtemplate_index;
    (*data)->strand = strand;

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
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int read_length,
    int subtemplate_index,
    enum STRAND strand,
    int location_offset
)
{
    /* find the correct error data record */
    struct error_data_record_t* record = data->records[
        2*subtemplate_index + ((int)strand - 1)];
    record->num_unique_reads += 1;
    
    int i;    
    for( i = 0; i < read_length; i++ )
    {
        if( toupper(read[i]) != toupper(genome_seq[i])  )
        {
            // Add 1 because read positions are 1 based
            unsigned char error_char = (unsigned char) error_str[i];
            record->base_type_mismatch_cnts[
                error_char][i+1+location_offset] += 1;
        }
        
        record->base_type_cnts[
            (unsigned char) error_str[i]][i+1+location_offset] += 1;
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

/*******************************************************************************
 *
 *
 * Error data updating code
 *
 *
 ******************************************************************************/

void*
update_error_data_from_index_search_results(
    struct read_subtemplate* rst,
    mapped_locations** search_results, 
    struct genome_data* genome, 
    struct error_data_t* error_data)
{
    /* Update the error data record */

    /* find the best candidate mapping */
    double lowest_penalty = -1e9;
    mapped_location* best_mapped_location = NULL;
    struct indexable_subtemplate* best_ist = NULL;
    int num_lowest_penalties = 0;
    
    int i;
    for( i = 0; NULL != search_results[i]; i++ ) 
    {
        /* if there aren't any results, there is nothing to do */
        if (search_results[i]->length == 0) continue;
        
        
        struct indexable_subtemplate* ist = search_results[i]->probe;
        
        int j;
        for ( j = 0; j < search_results[i]->length; j++ ) 
        {
            /* Find the location for the full read corresponding to 
               this mapped location, assuming that it is ungapped */
            mapped_location* curr_loc = &(search_results[i]->locations[j]);
            
            /* skip pseudo locations */
            if(curr_loc->chr == PSEUDO_LOC_CHR_INDEX) continue;
            
            int read_location
                = modify_mapped_read_location_for_index_probe_offset(
                    curr_loc->loc, curr_loc->chr, curr_loc->strnd,
                    ist->subseq_offset, ist->subseq_length,
                    rst->length, genome );
            if( 0 > read_location ) continue;
            
            char* genome_seq = find_seq_ptr( 
                genome, curr_loc->chr, read_location, rst->length );
            
            struct penalty_t* pa;
            if(curr_loc->strnd == FWD) {
                pa = rst->fwd_penalty_array->array;
            } else {
                pa = rst->rev_penalty_array->array;
            }
            float penalty = recheck_penalty(
                    genome_seq,
                    pa,
                    rst->length
                );
            
            /* if this is the best, then set it as such */
            if( penalty >= lowest_penalty - 1e-6 )
            {
                /* If we don't have a valid match yet or the new penalty
                   is strictly greater than then previous, then update */
                if( best_ist == NULL || 
                    penalty - 1e-6 >= lowest_penalty )
                {
                    best_mapped_location = search_results[i]->locations + j;
                    best_ist = ist;
                    lowest_penalty = penalty;
                    num_lowest_penalties = 1;
                } 
                /* if the penalties are the same, choose one randomly using
                   a resevoir sampling scheme (e.g. knuth)*/
                else {
                    num_lowest_penalties += 1;
                    if(rand()%num_lowest_penalties == num_lowest_penalties - 1)
                    {
                        best_mapped_location = search_results[i]->locations + j;
                        best_ist = ist;
                    }
                }
                
            }
        }
    }
    
    if( NULL == best_mapped_location ) {
        return NULL;
    }
    
    char* error_str = rst->error_str + best_ist->subseq_offset;
            
    char* fwd_genome_seq = find_seq_ptr( 
        genome, 
        best_mapped_location->chr, 
        best_mapped_location->loc,
        best_ist->subseq_length
        );
        
    char* genome_seq;
    if( best_mapped_location->strnd == BKWD )
    {
        genome_seq = calloc( best_ist->subseq_length+1, sizeof(char) );
        rev_complement_read(fwd_genome_seq, genome_seq, best_ist->subseq_length);
    } else {
        genome_seq = fwd_genome_seq;
    }
        
    update_error_data( 
        error_data, 
        genome_seq, 
        best_ist->char_seq, 
        error_str, 
        best_ist->subseq_length, 
        rst->pos_in_template.pos,
        best_mapped_location->strnd,
        best_ist->subseq_offset );
        
    if( best_mapped_location->strnd == BKWD )
    {
        free(genome_seq);
    }

    return NULL;
}

#ifdef INCREMENTLY_UPDATE_ERROR_MODEL

/* 
   Returns true if these candidate mappigns can be used to update the error 
   data. Basically, we just test for uniqueness. 
*/
static inline enum bool
can_be_used_to_update_error_data(
    candidate_mappings* mappings,
    struct genome_data* genome
)
{
    /* skip empty mappings */
    if( NULL == mappings )
        return false;
    
    /*** we only want unique mappers for the error estiamte updates */        
    // We allow lengths of 2 because we may have diploid locations
    if( mappings->length < 1 || mappings->length > 2 ) {
        return false;
    }
    
    /* we know that the length is at least 1 from directly above */
    assert( mappings->length >= 1 );
    
    candidate_mapping* loc = mappings->mappings + 0; 
    int mapped_length = loc->mapped_length;
    
    char* genome_seq = find_seq_ptr( 
        genome, 
        loc->chr, 
        loc->start_bp, 
        mapped_length
    );
        
    /* 
       if there are two mappings, make sure that they have the same 
       genome sequence. They actually may not be mapping to the same, 
       corresponding diploid locations, but as long as there isn't a second
       competing sequence, we really don't care because the mutation rates 
       should still be fine.
    */
    if( mappings->length > 1 )
    {
        /* this is guaranteed at the start of the function */
        assert( mappings->length == 2 );

        char* genome_seq_2 = find_seq_ptr( 
            genome, 
            mappings->mappings[1].chr, 
            mappings->mappings[1].start_bp, 
            mapped_length
        );
        
        /* if the sequences aren't identical, then return */
        if( 0 != strncmp( genome_seq, genome_seq_2, mapped_length ) )
            return false;
    }
    
    return true;
}

static inline void
update_error_data_record_from_candidate_mappings(
    struct genome_data* genome,
    candidate_mappings* mappings,
    struct read_subtemplate* rst,
    struct error_data_record_t* error_data_record
)
{
    if( !can_be_used_to_update_error_data( mappings, genome ) )
        return;
    
    /* we need at least one valid mapping */
    if( mappings->length == 0 )
        return;
    
    // use the index with the lowest penalty to update the error structure
    candidate_mapping* cm = mappings->mappings + 0;
    int i;
    for( i = 1; i < mappings->length; i++ ) {
        if( mappings->mappings[i].penalty > cm->penalty )
            cm = mappings->mappings + i;
    }
    
    // Ignore reads that are reverse complemented - we need to fix this
    if(cm->rd_strnd == BKWD) return; // XXX TODO

    int mapped_length = cm->mapped_length;
    
    char* genome_seq = find_seq_ptr( 
        genome, 
        cm->chr, 
        cm->start_bp,
        mapped_length
    );            
        
    /* get the read sequence - rev complement if on reverse strand */
    char* read_seq;
    if( cm->rd_strnd == BKWD )
    {
        read_seq = calloc( mapped_length + 1, sizeof(char) );
        rev_complement_read( rst->char_seq + cm->trimmed_length,
                read_seq, mapped_length );
    } else {
        read_seq = rst->char_seq + cm->trimmed_length;
    }
    
    /* Is this correct for rev comp? */
    char* error_str = rst->error_str + cm->trimmed_length;

    /* the offset is always 0 because this is a candidate mapping */
    update_error_data_record( 
        error_data_record, genome_seq, read_seq, error_str, mapped_length, 0 );
    
    /* free memory if we allocated it */
    if( cm->rd_strnd == BKWD )
        free( read_seq );
    
    return;
}

#endif
