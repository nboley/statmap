/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include <limits.h>
#include "math.h"

#include "statmap.h"
#include "find_candidate_mappings.h"
#include "quality.h"
#include "index_genome.h"
#include "mapped_location.h"
#include "error_correction.h"
#include "diploid_map_data.h"
#include "genome.h" // TODO: or incorporate diploid stuff into index_genome?
#include "rawread.h"

const float untemplated_g_marginal_log_prb = -1.30103;

#define MAX_NUM_UNTEMPLATED_GS 3

int
find_optimal_subseq_offset( 
    struct rawread* r,
    /* store the desired subsequences length */
    int subseq_len
) {
    if( subseq_len > r->length ) {
        fprintf( stderr, "============ %i \t\t %i \n", subseq_len, r->length );
    }
    assert( subseq_len <= r->length );
    
    
    /* XXX for now, we just use the first subseq_len characters,
       so the offset is always 0 */
    if( subseq_len == r->length )
        return 0;

    if( r->assay == CAGE )
    {
        if( r->length - MAX_NUM_UNTEMPLATED_GS < subseq_len )
        {
            fprintf(stderr, "FATAL        : CAGE experiments need indexes that have probe lengths at least 3 basepairs short, to account for templated G's." );
            return 0;
        }
        return MAX_NUM_UNTEMPLATED_GS;
    }
    
    return 0;
};

void
search_index( struct index_t* index, 
              
              float min_match_penalty,
              float max_penalty_spread,
              mapped_locations** results,

              struct rawread* r,
              float* bp_mut_rates,

              float* lookuptable_position,
              float* inverse_lookuptable_position,
              float* reverse_lookuptable_position,
              float* reverse_inverse_lookuptable_position,

              enum INDEX_SEARCH_MODE mode 
    )
{
    /**** Prepare the read for the index search */
    /* 
       first, we need to find the subseq of the read that we want to 
       probe the index for. This is controlled by the index->seq_length.
    */
    int subseq_offset = find_optimal_subseq_offset( r, index->seq_length );
    int subseq_length = index->seq_length;
    
    /* Store a copy of the read */
    /* This read has N's replaced with A's, and might be RC'd */
    char* sub_read = calloc(subseq_length + 1, sizeof(char));
    assert( sub_read != NULL );
    /* note that the NULL ending is pre-set from the calloc */
    memcpy( sub_read, r->char_seq + subseq_offset, sizeof(char)*(subseq_length) );

    /* prepare the results container */
    init_mapped_locations( results );
    (*results)->subseq_len = subseq_length;
    (*results)->subseq_offset = subseq_offset;
    
    /** Deal with the read on the fwd strand */
    /* Store the translated sequences here */
    LETTER_TYPE *fwd_seq;
    fwd_seq = translate_seq( sub_read, subseq_length, &fwd_seq );
    /* If we couldnt translate it */
    if( fwd_seq == NULL )
    {
        // fprintf(stderr, "Could Not Translate: %s\n", r->char_seq);
        return;
    }
    assert( fwd_seq != NULL );
    
    /** Deal with the read on the opposite strand */
    LETTER_TYPE *bkwd_seq;
    char* tmp_read = calloc(subseq_length + 1, sizeof(char));
    rev_complement_read( sub_read, tmp_read, subseq_length );
    bkwd_seq = translate_seq( tmp_read, subseq_length, &bkwd_seq );
    // BUG why did I have this here?
    //replace_ns_inplace( tmp_read, r->length );
    assert( bkwd_seq != NULL );
    
    /* map the full read */
    find_matches_from_root( index, 
                            
                            min_match_penalty,
                            max_penalty_spread,
                            *results,

                            /* length of the reads */
                            subseq_length,
                            
                            /* the fwd stranded sequence */
                            fwd_seq, 
                            lookuptable_position,
                            inverse_lookuptable_position,
                            
                            /* the bkwd stranded sequence */
                            bkwd_seq, 
                            reverse_lookuptable_position,
                            reverse_inverse_lookuptable_position,
                            
                            bp_mut_rates,

                            mode
        );
    
    /* Free the allocated memory */
    free( fwd_seq );
    free( bkwd_seq );

    free( sub_read );
    free( tmp_read );

    return;
};

static inline void
make_assay_specific_corrections( struct rawread* r, 
                                 mapped_locations* results )
{
    /* for now, we only deal with cage here */
    if( r->assay != CAGE )
        return;
    
    unsigned int i;
    unsigned int original_results_length = results->length;
    for( i = 0; i < original_results_length; i++ )
    {
        mapped_location loc = results->locations[i];
        
        /* deal with untemplated G's */
        int j;
        for( j = 0; 
             j < MAX_NUM_UNTEMPLATED_GS
                 && ( r->char_seq[j] == 'G' || r->char_seq[j] == 'g' );
             j++
            )
        {
            GENOME_LOC_TYPE tmp_loc = loc.location;
            
            /* A location always needs to be greater than 0 */
            /* note that we use j+1 because j == 0 corresponds
               to having an untemplated g at index 1 */
            if( tmp_loc.loc < (j+1) ) {
                continue;
            }
            /* ELSE */
            /* shift the genomic location */
            tmp_loc.loc += (j+1);
            
            /* add this new location */
            add_mapped_location( 
                results, tmp_loc, loc.strnd, j+1, 
                loc.penalty + (j+1)*untemplated_g_marginal_log_prb
            );
        }        
    }
    
    return;
}

static inline void 
recheck_location( struct genome_data* genome, 
                  struct rawread* r, 
                  candidate_mapping* loc,
                  const float* const lookuptable_position,
                  const float* const inverse_lookuptable_position,
                  const float* const reverse_lookuptable_position,
                  const float* const reverse_inverse_lookuptable_position,
                  const float* const bp_mut_rates
)
{
    if( PSEUDO_LOC_CHR_INDEX == loc->chr ) {
        return;
    }
    
    float marginal_log_prb = 0;
    
    const float* correct_lookuptable_position = NULL;
    const float* correct_inverse_lookuptable_position = NULL;
    
    /* find a pointer to the sequence at this genomic location */
    
    int mapped_length = r->length - loc->trimmed_len;
    
    char* genome_seq = find_seq_ptr( 
        genome, 
        loc->chr, 
        loc->start_bp,
        mapped_length
    );    
    
    /* if the genome_seq pointer is null, the sequence isn't valid
       for some reason ( ie, it falls off the end of the chromosome). 
       In such cases, mark the location as invalid and continue */
    if( NULL == genome_seq )
    {
        loc->recheck = INVALID;
        return;
    }
                    
    char* mut_genome_seq = NULL;
    mut_genome_seq = malloc(sizeof(char)*(r->length+1));
    assert( mut_genome_seq != NULL ); 
                    
    if( BKWD == loc->rd_strnd )
    {
        rev_complement_read( genome_seq, mut_genome_seq, r->length );
        correct_lookuptable_position = reverse_lookuptable_position;
        correct_inverse_lookuptable_position = reverse_inverse_lookuptable_position;
    } else {
        memcpy( mut_genome_seq, genome_seq, sizeof(char)*r->length );
        mut_genome_seq[r->length] = '\0';
        correct_lookuptable_position = lookuptable_position + loc->trimmed_len;
        correct_inverse_lookuptable_position = inverse_lookuptable_position + loc->trimmed_len;
    }
    
    float rechecked_penalty = recheck_penalty( 
        mut_genome_seq, 
        // char* observed,
        r->char_seq + loc->trimmed_len,
        // const int seq_length,
        mapped_length,
        correct_lookuptable_position,
        correct_inverse_lookuptable_position,
        bp_mut_rates
    );

    /*
     // * DEBUG
    printf( "%.2f\t%.2f\t%i\n", loc->penalty, rechecked_penalty, (loc->rd_strnd == BKWD) );
    printf( "%.*s\t%.*s\t%.*s\n", 
            r->length, r->char_seq, 
            r->length, genome_seq, 
            r->length, mut_genome_seq );
    */
    
    loc->penalty = rechecked_penalty + marginal_log_prb;
    
    free( mut_genome_seq );
}

static inline void
recheck_locations(
    struct genome_data* genome, 

    struct rawread* r, 
    candidate_mappings* mappings,
    
    float min_match_penalty,
    float max_penalty_spread,
    
    const float* const lookuptable_position,
    const float* const inverse_lookuptable_position,
    const float* const reverse_lookuptable_position,
    const float* const reverse_inverse_lookuptable_position,
    const float* const bp_mut_rates
)
{
    /* 
     * Currently, everything should be set *except* the gene strand. 
     * This is because we dont know what it is. Therefore, we will 
     * add these to the db with that bit unset, and then during the 
     * merging stage add to the penalty and say that it was equally 
     * likely to have come from either gene strand. This corresponds
     * with us being equally certain that the read is from either gene.
     *
     * Also, we need to check that we dont have any low quality reads.
     * ( This is possible if the path search went awry, and we found 
     *   a low quality read before a high quality read )
     */

    /* 
     * We keep track of the max observed penalty so that we can filter
     * out penalties that are too low. The index will never return 
     * results that are strictly below the minimum penalty, but it may 
     * return  results below the relative penalty. ( Read the indexing
     *  header for details )
     */
    float max_penalty = min_match_penalty;
            
    /* first, if necessary, recheck the locations */
    int k;
    for( k = 0; k < mappings->length; k++ )
    {
        recheck_location( genome, r, mappings->mappings + k,
                          lookuptable_position,
                          inverse_lookuptable_position,
                          reverse_lookuptable_position,
                          reverse_inverse_lookuptable_position,
                          bp_mut_rates
            );
                    
        /* we may need to update the max penalty */
        if( (mappings->mappings + k)->penalty > max_penalty ) {
            max_penalty = (mappings->mappings + k)->penalty;
        }
    }

    // int k declared earlier 
    for( k = 0; k < mappings->length; k++ )
    {
        /* this should be optimized out */
        candidate_mapping* loc = mappings->mappings + k;

        /* 
         * We always need to do this because of the way the search queue
         * handles poor branches. If our brnach prediction fails we could
         * add a low probability read, and then go back and find a very 
         * good read which would invalidate the previous. Then, the read
         * may not belong, but it isnt removed in the index searching method.
         */
                                
        /* we check max_penalty_spread > -0.00001 to make sure that the
           likelihood ratio threshold is actually active */

        /* I set the safe bit to 1e-6, which is correct for most floats */
        if(  max_penalty_spread > -0.00001  ) 
        {
            assert( max_penalty_spread >= 0.0 );
            if(   loc->penalty < ( max_penalty - max_penalty_spread ) )
            {
                loc->recheck = INVALID;
            } 
        }
                
        /* make sure that the penalty isn't too low */
        if( loc->penalty < min_match_penalty  )
        {
            loc->recheck = INVALID;
        } 

        /* if it's passed all of the tests, it must be valid */
        if( loc->recheck != INVALID )
            loc->recheck = VALID;                
    }

    return;
}

static void
add_candidate_mapping_from_haploid (
        struct rawread*     r, 
        mapped_location*    result,
        mapped_locations*   results,
        candidate_mapping   cm,
        candidate_mappings** mappings,
        struct genome_data* genome
    )
{
    /* set the chr */
    cm.chr = (result->location).chr;

    /* set the location. */
    /* We need to play with this a bit to account
       for index probes that are shorter than the read. */
    /* Skip the pseudo chr, this wil be modified later, ( if ever actually ) */
    int read_location = (result->location).loc;
    if( (result->location).chr != PSEUDO_LOC_CHR_INDEX )
    {
        read_location = modify_mapped_read_location_for_index_probe_offset(
            (result->location).loc, (result->location).chr, result->strnd, 
            results->subseq_offset, results->subseq_len, r->length,
            genome
        );
    }
    if( read_location < 0 ) // the read location was invalid; skip this mapped_location
        return;
    cm.start_bp = read_location;

    /* add the candidate mapping */
    add_candidate_mapping( *mappings, &cm );
}

static void
add_candidate_mapping_from_diploid (
        struct rawread*     r, 
        mapped_location*    result,
        mapped_locations*   results,
        candidate_mapping   cm,
        candidate_mappings** mappings,
        struct genome_data* genome
    )
{
    /* add the paternal candidate mapping. */
    add_candidate_mapping_from_haploid(
        r, result, results, cm, mappings, genome
    );

    /*
     * If this read mapped to a pseudo location, don't add another candidate mapping
     * Expand it when we expand all of the pseudo locs, later
     */
    if( result->location.chr == PSEUDO_LOC_CHR_INDEX )
        return;

    /* named variables for clarity */
    int paternal_chr_index = (result->location).chr;
    int paternal_loc = (result->location).loc;

    /* look up maternal chr_index */
    char* prefix = get_chr_prefix( genome->chr_names[paternal_chr_index] );
    int maternal_chr_index = find_diploid_chr_index(
            genome, prefix, MATERNAL
        );
    assert( maternal_chr_index > 0 ); // pseudo chr is not allowed
    free( prefix );

    /* look up associated diploid map data structure */
    int map_data_index = get_map_data_index_from_chr_index(
            genome, paternal_chr_index
        );
    assert( map_data_index >= 0 );

    /* get maternal start pos from diploid index */
    /* locations offset because diploid index is 1-indexed, but statmap is 0-indexed */
    int maternal_start = find_diploid_locations(
            &(genome->index->diploid_maps->maps[map_data_index]),
            paternal_loc + 1
        ) - 1;
    assert( maternal_start >= 0 ); // if not, this isn't true shared sequence

    int read_location;
    if( (result->location).chr != PSEUDO_LOC_CHR_INDEX )
    {
        read_location = modify_mapped_read_location_for_index_probe_offset(
            maternal_start, maternal_chr_index, result->strnd, 
            results->subseq_offset, results->subseq_len, r->length,
            genome
        );
    }
    if( read_location < 0 ) // the read location was invalid; skip this mapped_location
        return;

    /* modify cm to be maternal complement of original paternal cm */
    cm.start_bp = read_location;
    cm.chr = maternal_chr_index;

    /* add maternal candidate mapping */
    add_candidate_mapping( *mappings, &cm );
}


/* build candidate mappings from mapped locations ( 
   the data structure that index lookups return  )    */
static inline void 
build_candidate_mappings_from_mapped_locations(
    struct genome_data* genome,
    struct rawread* r, 
    mapped_locations* results, 
    candidate_mappings** mappings,
    float max_penalty_spread
)
{
    int subseq_len = results->subseq_len;
    
    /****** Prepare the template candidate_mapping objects ***********/
    candidate_mapping template_candidate_mapping 
        = init_candidate_mapping_from_template( 
            r, max_penalty_spread 
        );
    
    assert( template_candidate_mapping.rd_type != 0 );

    /* Make sure the "subseq" is acutally shorter than the read */
    assert( subseq_len <= r->length );
            
    /***** COPY information from the index lookup into the result set
     * build and populate an array of candidate_mapping's. 
     */        
    init_candidate_mappings( mappings );
                        
    /* Loop over mapped locations returned by index lookup */
    int i;
    for( i = 0; i < results->length; i++ )
    {        
        /* 
         * I'm scribbling on the base relation, but it doesnt 
         * matter because the add will copy it and then I will
         * overwrite what I scribbled on anyways. 
         */
        /* hopefully this will be optimized out */
        mapped_location* result;
        result = results->locations + i;

        /*** Set read-dependent info (same for diploid and haploid) ***/

        /* set the strand */
        if( result->strnd == FWD )
        {
            template_candidate_mapping.rd_strnd = FWD;
        } else {
            assert( result->strnd == BKWD );
            template_candidate_mapping.rd_strnd = BKWD;
        }

        /* set metadata */
        template_candidate_mapping.penalty = result->penalty;
        template_candidate_mapping.subseq_offset = results->subseq_offset;
        template_candidate_mapping.trimmed_len = result->trim_offset;

        /* Build, verify, and add candidate mappings depending on type of loc */
        if( result->location.is_paternal && result->location.is_maternal )
            add_candidate_mapping_from_diploid(
                r, result, results, template_candidate_mapping, mappings, genome
            );
        else
            add_candidate_mapping_from_haploid(
                r, result, results, template_candidate_mapping, mappings, genome
            );
    }
    
    return;
}

static inline void
update_error_data_from_candidate_mappings(
    struct genome_data* genome,
    struct error_data_t* error_data,
    candidate_mappings* mappings,
    struct rawread* r
)
{
    /* skip empty mappings */
    if( NULL == mappings )
        return;
        
    /*** we only want unique mappers for the error estiamte updates */        
    if( mappings->length != 1 ) {
        return;
    }
        
    /* we know that the length is one from right above */
    assert( mappings->length == 1 );
    if( mappings->mappings[0].recheck != VALID ) {
        return;
    }
        
    // emphasize the array aspect with the + 0
    // but, since the length is exactly 1, we know that
    // we only need to deal with this read
    candidate_mapping* loc = mappings->mappings + 0; 
    assert( mappings->length == 1 && loc->recheck == VALID );
        
    int mapped_length = r->length - loc->trimmed_len;
        
    char* genome_seq = find_seq_ptr( 
        genome, 
        loc->chr, 
        loc->start_bp,
        mapped_length
    );            
        
    char* error_str = r->error_str + loc->trimmed_len;

    /* get the read sequence - rev complement if on reverse strand */
    char* read_seq;
    if( loc->rd_strnd == BKWD )
    {
        read_seq = calloc( r->length + 1, sizeof(char) );
        rev_complement_read( r->char_seq + loc->trimmed_len, read_seq, r->length );
    } else {
        read_seq = r->char_seq + loc->trimmed_len;
    }

    update_error_data( error_data, genome_seq, read_seq, error_str, mapped_length );

    /* free memory if we allocated it */
    if( loc->rd_strnd == BKWD )
        free( read_seq );

    return;
}

/* TODO - revisit the read length vs seq length distinction */
/* 
 * I use a struct for the parameters so that I can initialzie threads
 * with this function. 
 *
 */
void*
find_candidate_mappings( void* params )    
{
    /* 
     * recreate the struct parameters for readability
     * this should be optimized out 
     */

    struct single_map_thread_data* td = params;

    int thread_id = td->thread_id;

    struct genome_data* genome = td->genome;
    
    unsigned int* mapped_cnt = td->mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex = td->mapped_cnt_mutex;

    struct rawread_db_t* rdb = td->rdb;
    
    pthread_mutex_t* mappings_db_mutex = td->mappings_db_mutex;
    candidate_mappings_db* mappings_db = td->mappings_db;
    /* The minimum absolute penalty
       that a valid read can have */
    float min_match_penalty = td->min_match_penalty;
    /* The minimum difference between the lowest penalty read and a
       valid read, set <= -1 to disable */
    float max_penalty_spread = td->max_penalty_spread;

    /* Use global_error_data to compute error prbs (unless we're bootstrapping) */
    struct error_data_t* global_error_data = td->global_error_data;
    /* Store observed error data from the current thread's execution in scratch */
    struct error_data_t* scratch_error_data = td->scratch_error_data;
    /* Flag - are we bootstrapping error data, or mapping? */
    enum INDEX_SEARCH_MODE mode = td->mode;

    /* END parameter 'recreation' */

    assert( genome->index != NULL );

    clock_t start;
    start = clock();
    
    /* build the basepair mutation rate lookup table */
    float* bp_mut_rates;
    if( mode == BOOTSTRAP )
        determine_bp_mut_rates_for_bootstrap( &bp_mut_rates );
    else
        determine_bp_mut_rates( &bp_mut_rates );

    /* how often we print out the mapping status */
    #define MAPPING_STATUS_GRANULARITY 100000

    /******** cache the candidate mapping results ********/
    /* cache the candidate mappings so that we can add them ( or not ) together at the
       end of this mapping */
    /* We do this so that we can update the error estimates */
    int curr_read_index = 0;
    /* we need 2* the step size to account for paired end reads */
    candidate_mappings* candidate_mappings_cache[2*READS_STAT_UPDATE_STEP_SIZE];
    memset( candidate_mappings_cache, 0,
            sizeof( candidate_mappings* )*2*READS_STAT_UPDATE_STEP_SIZE );

    struct rawread* rawreads_cache[2*READS_STAT_UPDATE_STEP_SIZE];
    memset( rawreads_cache, 0,
            sizeof( struct rawread* )*2*READS_STAT_UPDATE_STEP_SIZE );

    readkey_t readkeys[2*READS_STAT_UPDATE_STEP_SIZE];
    memset( readkeys, 0,
            sizeof( readkey_t )*2*READS_STAT_UPDATE_STEP_SIZE );

    int max_read_length = 0;
    
    /* The current read of interest */
    readkey_t readkey;
    struct rawread *r1, *r2;
    /* 
     * While there are still mappable reads in the read DB. All locking is done
     * in the get next read functions. 
     */
    while( EOF != get_next_read_from_rawread_db(
               rdb, &readkey, &r1, &r2, td->max_readkey )  
         ) 
    {
        /* Check that this read is mappable */
        if( mode == SEARCH )
        {
            if( filter_rawread( r1, global_error_data ) ||
                filter_rawread( r2, global_error_data ) )
            {
                /* skip the unmappable read */
                continue;
            }
        }

        /* We dont memory lock mapped_cnt because it's read only and we dont 
           really care if it's wrong 
         */
        if( readkey > 0 && 0 == readkey%MAPPING_STATUS_GRANULARITY )
        {
            fprintf(stderr, "DEBUG       :  Mapped %u reads, %i successfully\n", 
                    readkey, *mapped_cnt);
        }
        
        /* consider both read pairs */
        int j = 0;
        struct rawread* reads[2] = { r1, r2 };
        for( j = 0; j < 2 && reads[j] != NULL; j++ )
        {
            struct rawread* r = reads[j];
                        
            /* Build the quality lookup tables */
            float* lookuptable_position = malloc(sizeof(float)*r->length);
            float* inverse_lookuptable_position = malloc(sizeof(float)*r->length);
            float* reverse_lookuptable_position = malloc(sizeof(float)*r->length);
            float* reverse_inverse_lookuptable_position = malloc(sizeof(float)*r->length);
            build_lookup_table_from_rawread(
                r,
                // pass error data unless we're in bootstrap mode
                (mode == BOOTSTRAP) ? NULL : global_error_data,
                lookuptable_position, 
                inverse_lookuptable_position,
                reverse_lookuptable_position, 
                reverse_inverse_lookuptable_position
            );

            /**** go to the index for mapping locations */
            mapped_locations *results = NULL;

            search_index( genome->index, 

                          min_match_penalty,
                          max_penalty_spread,

                          &results,

                          r,
                          bp_mut_rates,

                          lookuptable_position,
                          inverse_lookuptable_position,
                          reverse_lookuptable_position,
                          reverse_inverse_lookuptable_position,

                          mode
                );

            /* if bootstrapping, we only want to work with unique mappers.
             * find_matches would have terminated early for this read */
            if( mode == BOOTSTRAP && results->length == 0 )
                /* bootstrap mode terminated early - cleanup and continue */
                goto cleanup;

            /* make an assay specific changes to the results. For instance, in CAGE,
               we need to add extra reads for the untemplated g's */
            make_assay_specific_corrections( r, results );
            
            /* make a reference to the current set of mappings. This should be
               optimized out by the compiler */
            candidate_mappings* mappings;

            /* cache current readkey and reads */
            readkeys[2*curr_read_index + j] = readkey;
            rawreads_cache[ 2*curr_read_index + j ] = r;

            build_candidate_mappings_from_mapped_locations(
                genome, r, results, 
                &mappings,
                min_match_penalty
            );

            /****** Do the recheck ******/
            
            // call recheck
            recheck_locations( 
                genome, 
                               
                r, mappings,
                               
                min_match_penalty,
                max_penalty_spread,
                               
                lookuptable_position,
                inverse_lookuptable_position,
                reverse_lookuptable_position,
                reverse_inverse_lookuptable_position,
                bp_mut_rates 
            );

            /* cache candidate mappings */
            assert( 2*curr_read_index + j < 2*READS_STAT_UPDATE_STEP_SIZE );
            candidate_mappings_cache[ 2*curr_read_index + j ] = mappings;

            /* update the maximum read length */
            max_read_length = MAX( max_read_length, r->length );

cleanup:
            /* free the quality score lookup tables */
            free( lookuptable_position );
            free( inverse_lookuptable_position );
            free( reverse_lookuptable_position );
            free( reverse_inverse_lookuptable_position );

            /* free the mapped_locations */
            free_mapped_locations( results );
        }

        curr_read_index += 1;
    }

    /******* update the error estimates *******/
    struct error_data_t* error_data;
    init_error_data( &error_data, 0, NULL );
    int i;
    for( i = 0; i < 2*READS_STAT_UPDATE_STEP_SIZE; i++ ) {
        update_error_data_from_candidate_mappings(
            genome,
            error_data,
            candidate_mappings_cache[ i ],
            rawreads_cache[ i ]
        );
    }
    // add error_data to scratch_error_data
    pthread_mutex_lock( scratch_error_data->mutex );
    add_error_data( scratch_error_data, error_data );
    pthread_mutex_unlock( scratch_error_data->mutex );
    // free local copy of error data
    free_error_data( error_data );

    /* cleanup the bp mutation rates */
    free( bp_mut_rates );
    
    /******* add the results to the database *******/
    for( i = 0; i < 2*READS_STAT_UPDATE_STEP_SIZE; i++ )
    {
        /* skip empty mapping lists */
        if( NULL == candidate_mappings_cache[i] )
            continue;

        candidate_mappings* mappings = candidate_mappings_cache[i];

        /* if we're running in search mode, add reads to the candidate mappings db and update mapped_cnts */
        if( mode == SEARCH )
        {
            /* increment the number of reads that mapped, if any pass the rechecks */
            int k;
            for( k = 0; k < mappings->length; k++ ) {
                if( (mappings->mappings + k)->recheck == VALID )
                {
                    pthread_mutex_lock( mapped_cnt_mutex );
                    *mapped_cnt += 1;
                    pthread_mutex_unlock( mapped_cnt_mutex );
                    break;
                }
            }

            /* add cms to the db */
            pthread_mutex_lock( mappings_db_mutex );
            assert( thread_id < num_threads );
            /* note that we add to the DB even if there are 0 that map,
               we do this so it is easier to join with the rawreads later */
            add_candidate_mappings_to_db(
                mappings_db, mappings, readkeys[i], thread_id );
            pthread_mutex_unlock( mappings_db_mutex );
        }

        /* if mode is BOOTSTRAP, then we still loop over the cache so we can free
           the allocated memory - but we already have what we came for (error_data) */

        /* free the cached reads and mappings */
        free_rawread( rawreads_cache[i] );
        free_candidate_mappings( mappings );
    }
    
    return NULL;
}

void
spawn_threads( struct single_map_thread_data* td_template )
{
    long t;
    int rc;
    void* status;
    pthread_t thread[num_threads];

    pthread_attr_t attrs[num_threads];
    
    struct single_map_thread_data tds[num_threads];
    
    for( t = 0; t < num_threads; t++ )
    {  
        memcpy( tds+t,  td_template, sizeof(struct single_map_thread_data) );
        tds[t].thread_id = t;
        
        pthread_attr_init(attrs + t);
        pthread_attr_setdetachstate(attrs + t, PTHREAD_CREATE_JOINABLE);
        
        rc = pthread_create( thread + t, 
                             attrs + t, 
                             find_candidate_mappings, 
                             (void *)(tds + t) 
            ); 
        
        if (rc) {
            fprintf(stderr, 
                    "ERROR; return code from pthread_create() is %d\n", 
                    rc
            );
            exit(-1);
        }
    }
    
    /* Free attribute and wait for the other threads */    
    size_t num_reads2 = 0;
    for(t=0; t < num_threads; t++) {
        rc = pthread_join(thread[t], &status);
        pthread_attr_destroy(attrs+t);
        if (rc) {
            fprintf( stderr, 
                     "ERROR; return code from pthread_join() is %d\n", 
                     rc
            );
            exit(-1);
        }
        num_reads2 += (size_t) status;
    }
}


/*
 * Find all locations that a read could map to, given thresholds,
 * reads, a genome and an index.
 *
 */

void
find_all_candidate_mappings( struct genome_data* genome,
                             FILE* log_fp,
                             struct rawread_db_t* rdb,

                             candidate_mappings_db* mappings_db,
                             float min_match_penalty,
                             float max_penalty_spread,
                             float max_seq_length

    )
{
    clock_t start = clock();

    /* 
     * init mutexes to guard access to the log_fp, the fastq 'db', and
     * the mappings db. They will be stored in the same structure as the 
     * passed data
     */

    /* initialize the necessary mutex's */
    pthread_mutex_t log_fp_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mapped_cnt_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mappings_db_mutex = PTHREAD_MUTEX_INITIALIZER;

    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    td_template.genome = genome;
    
    td_template.log_fp = log_fp;
    td_template.log_fp_mutex = &log_fp_mutex;

    unsigned int mapped_cnt = 0;
    td_template.mapped_cnt = &mapped_cnt;
    td_template.mapped_cnt_mutex = &mapped_cnt_mutex;
    td_template.max_readkey = 0;

    td_template.rdb = rdb;
    
    td_template.mappings_db = mappings_db;
    td_template.mappings_db_mutex = &mappings_db_mutex;
    
    td_template.min_match_penalty = min_match_penalty;
    td_template.max_penalty_spread = max_penalty_spread;
    td_template.max_subseq_len = max_seq_length;

    /*
     * initialize global error_data and scratch error_data
     * scratch error_data has a mutex so it can be safely updated by each thread
     */
    struct error_data_t* global_error_data;
    init_error_data( &global_error_data, 0, NULL );
    td_template.global_error_data = global_error_data;

    struct error_data_t* scratch_error_data;
    pthread_mutex_t scratch_err_mutex = PTHREAD_MUTEX_INITIALIZER;
    init_error_data( &scratch_error_data, 0, &scratch_err_mutex );
    td_template.scratch_error_data = scratch_error_data;

    /*
       bootstrap global error data
     */
    /* bootstrap error data from first READS_STAT_UPDATE_STEP_SIZE reads */
    td_template.max_readkey = READS_STAT_UPDATE_STEP_SIZE;
    td_template.mode = BOOTSTRAP;
    spawn_threads( &td_template );

    /* set global_error_data to bootstrap's averaged scratch_error_data, and reset scratch */
    //average_error_data( scratch_error_data );
    add_error_data( global_error_data, scratch_error_data );
    //global_error_data->num_unique_reads = 0;
    clear_error_data( td_template.scratch_error_data );

    /*
     * check that there were some unique mappers - warn if no
     */
    if( !(global_error_data->num_unique_reads > 0) )
    {
        fprintf( stderr, "FATAL       :  Statmap could not map any unique reads in the bootstrap.\n" );
        exit(1);
    }

    /* 
       map reads with the bootstrapped error data
       reset td_template parameters to initial states
     */
    td_template.max_readkey = 0;
    td_template.mode = SEARCH;
    /* rewind rawread db to beginning and now we'll actually map */
    rewind_rawread_db( rdb );

    /* initialize the threads */
    while( false == rawread_db_is_empty( rdb ) )
    {
        // update the maximum allowable readkey
        td_template.max_readkey += READS_STAT_UPDATE_STEP_SIZE;

        /* log the error data we're using for this round of mapping */
        log_error_data( global_error_data );

        spawn_threads( &td_template );

        /* after threads are done, average scratch data over the reads processed
           by all threads, then weighted average into global_error_data */
        //average_error_data( scratch_error_data );
        update_global_error_data( global_error_data, scratch_error_data );
    }

    // free error data structs
    free_error_data( scratch_error_data );
    free_error_data( global_error_data );

    /* Find all of the candidate mappings */    
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Mapped (%i/%u) Partial Reads in %.2lf seconds ( %e/thread-hour )\n",
            mapped_cnt, rdb->readkey, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)mapped_cnt)*CLOCKS_PER_SEC*3600)/(stop-start)
        );
}
