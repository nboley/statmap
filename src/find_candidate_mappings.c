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
#include "genome.h"
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
            fprintf(stderr, "FATAL        : CAGE experiments need indexes that "
                            "have probe lengths at least 3 basepairs short, to "
                            "account for untemplated G's." );
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

              struct penalty_array_t* fwd_pa,
              struct penalty_array_t* rev_pa,

              enum bool only_find_unique_mappers
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
                            /* the bkwd stranded sequence */
                            bkwd_seq, 

                            fwd_pa,
                            rev_pa,

                            only_find_unique_mappers
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
recheck_location(
        struct genome_data* genome, 
        struct rawread* r, 
        candidate_mapping* loc,

        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa
    )
{
    if( PSEUDO_LOC_CHR_INDEX == loc->chr ) {
        return;
    }
    
    float marginal_log_prb = 0;
    
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

    struct penalty_array_t* penalty_array;
                    
    if( BKWD == loc->rd_strnd )
    {
        rev_complement_read( genome_seq, mut_genome_seq, r->length );
        penalty_array = rev_pa;
    } else {
        memcpy( mut_genome_seq, genome_seq, sizeof(char)*r->length );
        mut_genome_seq[r->length] = '\0';
        penalty_array = fwd_pa;
    }
    
    float rechecked_penalty = recheck_penalty( 
            mut_genome_seq, 
            r->char_seq + loc->trimmed_len,
            mapped_length,
            penalty_array 
        );

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
    
    struct penalty_array_t* fwd_pa,
    struct penalty_array_t* rev_pa 

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
        recheck_location( genome, r, mappings->mappings + k, fwd_pa, rev_pa );
                    
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

/* 
   Returns true if these candidate mappigns can be used to update the error 
   data. Basically, we just test for uniqueness. 
*/
static inline enum bool
can_be_used_to_update_error_data(
    struct genome_data* genome,
    candidate_mappings* mappings,
    struct rawread* r
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
    int mapped_length = r->length - loc->trimmed_len;
    
    if( loc->recheck != VALID ) {
        return false;
    }
    
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
        if( mappings->mappings[1].recheck != VALID ) {
            return false;
        }
        
        if( mappings->mappings[1].trimmed_len != loc->trimmed_len )
            return false;
        
        /* this is guaranteed at the start of the function */
        assert( mappings->length == 2 );
        char* genome_seq_2 = find_seq_ptr( 
            genome, 
            mappings->mappings[1].chr, 
            mappings->mappings[1].start_bp, 
            r->length - mappings->mappings[1].trimmed_len
        );
        
        /* if the sequences aren't identical, then return */
        if( 0 != strncmp( genome_seq, genome_seq_2, mapped_length ) )
            return false;
    }
    
    return true;
}

static inline void
update_error_data_from_candidate_mappings(
    struct genome_data* genome,
    candidate_mappings* mappings,
    struct rawread* r,
    struct error_data_t* error_data
)
{
    if( !can_be_used_to_update_error_data( genome, mappings, r ) )
        return;
    
    /* we need at least one valid mapping ( although this case should be handled
       above by can_be_used_to_update_error_data */
    if( mappings->length == 0 )
        return;
    
    // emphasize the array aspect with the + 0
    // but, since the length is exactly 1, we know that
    // we only need to deal with this read
    candidate_mapping* loc = mappings->mappings + 0;         
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
        rev_complement_read( r->char_seq+loc->trimmed_len, read_seq, r->length);
    } else {
        read_seq = r->char_seq + loc->trimmed_len;
    }

    update_error_data( 
        error_data, genome_seq, read_seq, error_str, mapped_length );

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

    /* Store observed error data from the current thread's execution in scratch */
    struct error_model_t* error_model = td->error_model;
    
    /* store statistics about mapping quality here  */
    struct error_data_t* error_data = td->error_data;
    
    /* if we only want error data, then there is not reason to find antyhing 
       except unqiue reads. */
    enum bool only_find_unique_mappers = td->only_collect_error_data;
    
    /* END parameter 'recreation' */

    assert( genome->index != NULL );

    clock_t start;
    start = clock();
    
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
        /* We dont memory lock mapped_cnt because it's read only and we dont 
           really care if it's wrong 
         */
        if( readkey > 0 && 0 == readkey%MAPPING_STATUS_GRANULARITY )
        {
            fprintf(stderr, "DEBUG       :  Mapped %u reads, %i successfully\n", 
                    readkey, *mapped_cnt);
        }
        
        if( filter_rawread( r1, error_model ) ||
            filter_rawread( r2, error_model ) )
        {
            /* skip the unmappable read */
            continue;
        }
        
        /* consider both read pairs */
        int j = 0;
        struct rawread* reads[2] = { r1, r2 };
        for( j = 0; j < 2 && reads[j] != NULL; j++ )
        {
            struct rawread* r = reads[j];
                        
            /*
               initialize read penalty array
            */
            struct penalty_array_t fwd_pa, rev_pa;
            init_penalty_array( r->length, &fwd_pa );
            init_penalty_array( r->length, &rev_pa );

            /* fn ptr to penalty array builder - build different penalty
               array depending on type of search */
            build_penalty_array( r, error_model, &fwd_pa );
            build_reverse_penalty_array( &fwd_pa, &rev_pa );
            
            /**** go to the index for mapping locations */
            mapped_locations *results = NULL;

            search_index( genome->index, 

                          min_match_penalty,
                          max_penalty_spread,

                          &results,

                          r,

                          &fwd_pa,
                          &rev_pa,

                          only_find_unique_mappers
                );

            /* cache current readkey and reads */
            readkeys[2*curr_read_index + j] = readkey;
            rawreads_cache[ 2*curr_read_index + j ] = r;
            
            /* if bootstrapping, we only want to work with unique mappers.
             * find_matches would have terminated early for this read */
            if( results->skip )
                /* bootstrap mode terminated early - cleanup and continue */
                goto cleanup_penalty_arrays;

            /* make an assay specific changes to the results. For instance, in CAGE,
               we need to add extra reads for the untemplated g's */
            make_assay_specific_corrections( r, results );
            
            /* make a reference to the current set of mappings. This should be
               optimized out by the compiler */
            candidate_mappings* mappings;
            
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
                               
                &fwd_pa,
                &rev_pa
            );

            /* cache candidate mappings */
            assert( 2*curr_read_index + j < 2*READS_STAT_UPDATE_STEP_SIZE );
            candidate_mappings_cache[ 2*curr_read_index + j ] = mappings;

            /* update the maximum read length */
            max_read_length = MAX( max_read_length, r->length );

cleanup_penalty_arrays:
            /* free the penalty arrays */
            free_penalty_array( &fwd_pa );
            free_penalty_array( &rev_pa );

            /* free the mapped_locations */
            free_mapped_locations( results );
        }

        curr_read_index += 1;
    }

    /******* update the error estimates *******/
    /* create a thread local copy of the error data to avoid excess locking */
    struct error_data_t* scratch_error_data;
    init_error_data( &scratch_error_data );

    int i;
    for( i = 0; i < 2*READS_STAT_UPDATE_STEP_SIZE; i++ ) {
        update_error_data_from_candidate_mappings(
            genome,
            candidate_mappings_cache[ i ],
            rawreads_cache[ i ],
            scratch_error_data
            );
    }
    // add error_data to scratch_error_data
    pthread_mutex_lock( error_data->mutex );
    add_error_data( error_data, scratch_error_data );
    pthread_mutex_unlock( error_data->mutex );
    
    // free local copy of error data
    free_error_data( scratch_error_data );
    
    if( td->only_collect_error_data )
        goto cleanup;
    
    /******* add the results to the database *******/
    for( i = 0; i < 2*READS_STAT_UPDATE_STEP_SIZE; i++ )
    {
        /* skip empty mapping lists */
        if( NULL == candidate_mappings_cache[i] )
            continue;

        candidate_mappings* mappings = candidate_mappings_cache[i];

        /* increment the number of reads that mapped, if 
           any pass the rechecks */
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
    
cleanup:

    for( i = 0; i < 2*READS_STAT_UPDATE_STEP_SIZE; i++ )
    {
        if( NULL == candidate_mappings_cache[i])
            continue;
        
        /* free the cached reads and mappings */
        free_rawread( rawreads_cache[i] );
        free_candidate_mappings( candidate_mappings_cache[i] );
    }
    
    return NULL;
}

void
spawn_threads( struct single_map_thread_data* td_template )
{
    if( num_threads == 1 )
    {
        find_candidate_mappings( td_template );
        return;
    }
    
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
 */

void init_td_template( struct single_map_thread_data* td_template,
                       struct genome_data* genome, 
                       struct rawread_db_t* rdb,
                       candidate_mappings_db* mappings_db,
                       struct error_model_t* error_model,
                       float min_match_penalty, float max_penalty_spread )
{
    /* initialize the necessary mutex's */
    /*
    pthread_mutex_t* mapped_cnt_mutex;
    mapped_cnt_mutex = malloc( sizeof(pthread_mutex_t) );
    (*mapped_cnt_mutex) = PTHREAD_MUTEX_INITIALIZER;

    pthread_mutex_t* mappings_db_mutex;
    mappings_db_mutex= malloc( sizeof(pthread_mutex_t) );
    *mappings_db_mutex = PTHREAD_MUTEX_INITIALIZER;
    */
    int rc;
    pthread_mutexattr_t mta;
    rc = pthread_mutexattr_init(&mta);
    
    td_template->genome = genome;
    
    td_template->mapped_cnt = malloc( sizeof(unsigned int) );
    *(td_template->mapped_cnt) = 0;
    
    td_template->mapped_cnt_mutex = malloc( sizeof(pthread_mutex_t) );
    rc = pthread_mutex_init( td_template->mapped_cnt_mutex, &mta );
    td_template->max_readkey = 0;

    td_template->rdb = rdb;
    
    td_template->mappings_db = mappings_db;
    td_template->mappings_db_mutex = malloc( sizeof(pthread_mutex_t) );
    rc = pthread_mutex_init( td_template->mappings_db_mutex, &mta );
    
    td_template->min_match_penalty = min_match_penalty;
    td_template->max_penalty_spread = max_penalty_spread;

    init_error_data( &(td_template->error_data) );

    td_template->error_model = error_model;

    td_template->max_readkey = READS_STAT_UPDATE_STEP_SIZE;
    
    td_template->thread_id = 0;
    
    /* Make sure we are finding all mappers for a normal index search */
    td_template->only_collect_error_data = false;
    
    return;
}

void 
free_td_template( struct single_map_thread_data* td_template )
{
    free( td_template->mapped_cnt );
    free( td_template->mapped_cnt_mutex );
    free( td_template->mappings_db_mutex );
    return;
}

/* bootstrap an already initialized error model */
void
bootstrap_estimated_error_model( 
    struct genome_data* genome,
    
    struct rawread_db_t* rdb,
    candidate_mappings_db* mappings_db,
    
    struct error_model_t* error_model
) 
{
    assert( error_model != NULL );
    
    struct error_model_t* bootstrap_error_model;
    init_error_model( &bootstrap_error_model, MISMATCH );
    
    /* put the search arguments into a structure */
    #define MAX_NUM_MM 5
    #define MAX_MM_SPREAD 3
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mappings_db, 
                      bootstrap_error_model, MAX_NUM_MM, MAX_MM_SPREAD );
    
    /* 
       only use unique mappers for the initial bootstrap. This is just a small
       performance optimization, it prevents us from going too deeply into the 
       index as soon as we know that a mapping isn't unique.
    */
    td_template.only_collect_error_data = true;
    spawn_threads( &td_template );
    
    update_error_model_from_error_data( error_model, td_template.error_data );
    
    free_error_model( bootstrap_error_model );
    
    free_td_template( &td_template );
    
    return;
}

void
find_all_candidate_mappings( struct genome_data* genome,

                             struct rawread_db_t* rdb,
                             candidate_mappings_db* mappings_db,
                             
                             struct error_model_t* error_model,
                             
                             float min_match_penalty,
                             float max_penalty_spread
    )
{
    clock_t start = clock();

    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mappings_db, error_model,
                      min_match_penalty, max_penalty_spread );
    
    /* initialize the threads */
    while( false == rawread_db_is_empty( rdb ) )
    {
        spawn_threads( &td_template );

        // update the maximum allowable readkey
        td_template.max_readkey += READS_STAT_UPDATE_STEP_SIZE;

        /* log the error data we're using for this round of mapping */
        log_error_data( td_template.error_data );

        /* update the error model from the new error data */
        update_error_model_from_error_data( error_model, td_template.error_data );
    }
    
    /* Print out performance information */
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Mapped (%i/%u) Partial Reads in %.2lf seconds ( %e/thread-hour )\n",
            *(td_template.mapped_cnt), rdb->readkey, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)*(td_template.mapped_cnt))*CLOCKS_PER_SEC*3600)/(stop-start)
        );

    free_td_template( &td_template );

    return;
}
