/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "statmap.h"
#include "find_candidate_mappings.h"
#include "quality.h"
#include "index_genome.h"
#include "mapped_location.h"
#include "error_correction.h"
#include "diploid_map_data.h"
#include "genome.h"
#include "rawread.h"
#include "read.h"
#include "mapped_read.h"
#include "pseudo_location.h"

const float untemplated_g_marginal_log_prb = -1.30103;

#define MAX_NUM_UNTEMPLATED_GS 3

float
subseq_penalty(
        struct read_subtemplate* rst,
        int offset,
        int subseq_len,
        struct penalty_array_t* penalty_array
    )
{
    float penalty = 0;

    /* loop over subsequence */
    int pos;
    for( pos = offset; pos < subseq_len; pos++ )
    {
        int bp = bp_code( rst->char_seq[pos] );
        /* take the match penalty for this bp */
        penalty += penalty_array->array[pos].penalties[bp][bp];
    }

    return penalty;
}

int
find_optimal_subseq_offset( 
    struct read_subtemplate* rst,
    struct penalty_array_t* penalty_array,

    int subseq_len,

    /* define region of underlying read to search for optimal subsequences */
    int region_start,
    int region_length
) {
    /* Make sure the read is at least as long as the subsequence */
    if( subseq_len > rst->length ) {
        fprintf( stderr, "============ %i \t\t %i \n", subseq_len, rst->length );
    }
    assert( subseq_len <= rst->length );

    /* Make sure the search region makes sense */
    assert( region_start >= 0 && region_start <= (rst->length - subseq_len) );
    assert( region_length <= rst->length );
    assert( (region_start + region_length - subseq_len) >= 0 );
    
    /*
       Remember: error_prb returns the inverse log probability of
       error: log10(1 - P(error)) for matches
    */
    int min_offset = 0;
    float max_so_far = -FLT_MAX;

    int i;
    /* each possible start bp in the subsequence */
    for( i = region_start;
         i < (region_start + region_length - subseq_len);
         i++ )
    {
        float ss_pen = subseq_penalty(rst, i, subseq_len, penalty_array);
        if( ss_pen > max_so_far ) {
            max_so_far = ss_pen;
            min_offset = i;
        }
    }

    return min_offset;
};

void
search_index(
        struct genome_data* genome,
        struct indexable_subtemplate* ist,

        float min_match_penalty,
        float max_penalty_spread,
        mapped_locations** results,

        enum bool only_find_unique_sequence
    )
{
    // reference to index
    struct index_t* index = genome->index;

    int subseq_length = index->seq_length;
    
    /* prepare the results container */
    init_mapped_locations( results, ist );
    
    /* Build bitpacked copies of the fwd and rev strand versions of this
     * indexable subtemplate */

    /* Store a copy of the read */
    /* This read has N's replaced with A's, and might be RC'd */
    char* sub_read = calloc(subseq_length + 1, sizeof(char));
    assert( sub_read != NULL );
    /* note that the NULL ending is pre-set from the calloc */
    memcpy( sub_read, ist->char_seq,
            sizeof(char)*(subseq_length) );

    /** Deal with the read on the fwd strand */
    /* Store the translated sequences here */
    LETTER_TYPE *fwd_seq;
    fwd_seq = translate_seq( sub_read, subseq_length );
    /* If we couldnt translate it */
    if( fwd_seq == NULL )
    {
        // fprintf(stderr, "Could Not Translate: %s\n", st->char_seq);
        return;
    }
    assert( fwd_seq != NULL );
    
    /** Deal with the read on the opposite strand */
    LETTER_TYPE *bkwd_seq;
    char* tmp_read = calloc(subseq_length + 1, sizeof(char));
    rev_complement_read( sub_read, tmp_read, subseq_length );
    bkwd_seq = translate_seq( tmp_read, subseq_length );
    assert( bkwd_seq != NULL );
    
    /* map the full read */
    find_matches_from_root(
            index, 

            min_match_penalty,
            max_penalty_spread,
            *results,

            genome,

            /* length of the reads */
            subseq_length,

            /* the fwd stranded sequence */
            fwd_seq, 
            /* the bkwd stranded sequence */
            bkwd_seq, 

            ist->fwd_penalties,
            ist->rev_penalties,

            only_find_unique_sequence
        );

    /* Cleanup memory */
    free( fwd_seq );
    free( bkwd_seq );

    free( sub_read );
    free( tmp_read );

    return;
};

static inline void 
recheck_location(
        struct genome_data* genome, 
        struct read_subtemplate* rst,
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
    int mapped_length = rst->length;
    
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
    mut_genome_seq = malloc(sizeof(char)*(rst->length+1));
    assert( mut_genome_seq != NULL ); 

    struct penalty_array_t* penalty_array;
                    
    if( BKWD == loc->rd_strnd )
    {
        rev_complement_read( genome_seq, mut_genome_seq, rst->length );
        penalty_array = rev_pa;
    } else {
        memcpy( mut_genome_seq, genome_seq, sizeof(char)*rst->length );
        mut_genome_seq[rst->length] = '\0';
        penalty_array = fwd_pa;
    }
    
    float rechecked_penalty = recheck_penalty( 
            mut_genome_seq, 
            rst->char_seq,
            mapped_length,
            penalty_array 
        );

    loc->penalty = rechecked_penalty + marginal_log_prb;
    
    free( mut_genome_seq );
}

static inline void
recheck_locations(
    struct genome_data* genome, 

    struct read_subtemplate* rst,
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
        recheck_location( genome, rst, mappings->mappings + k, fwd_pa, rev_pa );
                    
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

static inline void
build_candidate_mapping_from_mapped_location(
        struct genome_data* genome,
        struct read_subtemplate* rst,

        mapped_location* result, 
        mapped_locations* results,

        candidate_mappings* mappings,
        float max_penalty_spread
    )
{
    int subseq_len = results->probe->subseq_length;

    /****** Prepare the template candidate_mapping objects ***********/
    candidate_mapping cm 
        = init_candidate_mapping_from_template( 
            rst, max_penalty_spread 
        );

    /* Make sure the "subseq" is acutally shorter than the read */
    assert( subseq_len <= rst->length );

    /* set the strand */
    assert( result->strnd == FWD || result->strnd == BKWD );
    cm.rd_strnd = result->strnd;

    /* set metadata */
    cm.penalty = result->penalty;
    cm.subseq_offset = results->probe->subseq_offset;

    /* set location information */
    cm.chr = result->chr;

    /* We need to play with this a bit to account for index probes that are
     * shorter than the read length */
    int read_location =
        modify_mapped_read_location_for_index_probe_offset(
            result->loc,
            result->chr,
            result->strnd,
            results->probe->subseq_offset,
            results->probe->subseq_length,
            rst->length,
            genome
        );
    if( read_location < 0 ) // the read location was invalid; skip this mapped_locations_container
        return;

    cm.start_bp = read_location;

    add_candidate_mapping( mappings, &cm );
}

/* 
   Returns true if these candidate mappigns can be used to update the error 
   data. Basically, we just test for uniqueness. 
*/
static inline enum bool
can_be_used_to_update_error_data(
    struct genome_data* genome,
    candidate_mappings* mappings,
    struct read_subtemplate* rst
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
    int mapped_length = rst->length;
    
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
        
        /* this is guaranteed at the start of the function */
        assert( mappings->length == 2 );
        char* genome_seq_2 = find_seq_ptr( 
            genome, 
            mappings->mappings[1].chr, 
            mappings->mappings[1].start_bp, 
            rst->length
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
    if( !can_be_used_to_update_error_data( genome, mappings, rst ) )
        return;
    
    /* we need at least one valid mapping ( although this case should be handled
       above by can_be_used_to_update_error_data */
    if( mappings->length == 0 )
        return;
    
    // emphasize the array aspect with the + 0
    // but, since aal of the sequence is identical,
    // we only need to deal with this read
    candidate_mapping* loc = mappings->mappings + 0;         
    int mapped_length = rst->length;
    
    char* genome_seq = find_seq_ptr( 
        genome, 
        loc->chr, 
        loc->start_bp,
        mapped_length
    );            
        
    char* error_str = rst->error_str;

    /* get the read sequence - rev complement if on reverse strand */
    char* read_seq;
    if( loc->rd_strnd == BKWD )
    {
        read_seq = calloc( rst->length + 1, sizeof(char) );
        rev_complement_read( rst->char_seq,
                read_seq, rst->length);
    } else {
        read_seq = rst->char_seq;
    }
    
    update_error_data_record( 
        error_data_record, genome_seq, read_seq, error_str, mapped_length );

    /* free memory if we allocated it */
    if( loc->rd_strnd == BKWD )
        free( read_seq );
    
    return;
}

void
build_indexable_subtemplate(
        struct read_subtemplate* rst,
        struct indexable_subtemplates* ists,
        struct index_t* index,

        struct penalty_array_t* fwd_penalty_array,
        struct penalty_array_t* rev_penalty_array,

        // area of the read subtemplate to take an indexable subtemplate from
        int range_start,
        int range_length
    )
{
    int subseq_length = index->seq_length;

    /* Find the optimal subsequence offset for this read subtemplate */
    int subseq_offset = find_optimal_subseq_offset(
            rst,
            fwd_penalty_array, // TODO - different offset for rev comp?
            subseq_length,
            range_start,
            range_length
        );

    struct indexable_subtemplate* ist = NULL;
    init_indexable_subtemplate( &ist,
            subseq_length,
            subseq_offset,
            rst->char_seq,
            fwd_penalty_array,
            rev_penalty_array
        );

    // copy indexable subtemplate into set of indexable subtemplates
    add_indexable_subtemplate_to_indexable_subtemplates( ist, ists );

    // free working copy
    free_indexable_subtemplate( ist );
}

void
build_indexable_subtemplates_from_read_subtemplate(
        struct indexable_subtemplates* ists,
        struct read_subtemplate* rst,
        struct index_t* index,
        struct penalty_array_t* fwd_penalty_array,
        struct penalty_array_t* rev_penalty_array
    )
{
    /*
       TODO for now, we build 2 indexable subtemplates if the index sequence
       length is <= the read subtemplate length / 2. Otherwise build a single
       indexable subtemplate
     */
    int subseq_length = index->seq_length;

    if( subseq_length <= rst->length / 2 )
    {
        /* we can build 2 non-overlapping indexable subtemplates */
        build_indexable_subtemplate(
                rst, ists, index,
                fwd_penalty_array, rev_penalty_array,
                0, rst->length / 2
            );

        build_indexable_subtemplate(
                rst, ists, index,
                fwd_penalty_array, rev_penalty_array,
                rst->length / 2, rst->length
            );
    }
    else
    {
        /* build a single indexable subtemplate */
        build_indexable_subtemplate(
                rst, ists, index,
                fwd_penalty_array, rev_penalty_array,
                0, rst->length
            );
    }
}

void
search_index_for_indexable_subtemplates(
        struct indexable_subtemplates* ists,
        mapped_locations_container* ml_container,

        struct genome_data* genome,

        float min_match_penalty,
        float max_penalty_spread,

        enum bool only_collect_error_data
    )
{
    // search the index for each indexable subtemplate
    int ist_index;
    for( ist_index = 0; ist_index < ists->length; ist_index++ )
    {
        // reference to current indexable subtemplate
        struct indexable_subtemplate* ist = &(ists->container[ist_index]);

        /**** go to the index for mapping locations */
        mapped_locations *results = NULL;
        search_index(
                genome,
                ist,

                min_match_penalty,
                max_penalty_spread,

                &results,

                only_collect_error_data
            );

        // add these mapped locations to the mapped locations container
        add_mapped_locations_to_mapped_locations_container(
                results, ml_container );
    }
}

int
choose_base_mapped_locations_index(
        mapped_locations_container* mls_container
    )
{
    /* TODO for now, take the shortest set of mapped locations */
    int min_so_far = INT_MAX;
    int min_index = 0;

    int i;
    for( i = 0; i < mls_container->length; i++ )
    {
        mapped_locations* current = mls_container->container[i];

        if( current->length < min_so_far )
        {
            min_so_far = current->length;
            min_index = i;
        }
    }

    return min_index;
}

int
find_start_of_pseudo_mapped_locations_for_strand(
        mapped_locations* sorted_mapped_locs,
        enum STRAND strand
    )
{
    int start = -1; // Use -1 to signal no locations with given strand

    int i;
    for( i = 0; i < sorted_mapped_locs->length; i++ )
    {
        /* Assumes the sorted_mapped_locs are, indeed, sorted by
         * sort_mapped_locations_by_location */
        if( sorted_mapped_locs->locations[i].strnd == strand )
        {
            start = i;
            break;
        }
    }

    return start;
}

void
search_for_matches_in_pseudo_locations(
        mapped_locations* matching_subset,
        mapped_locations* potential_matches,
        mapped_location* key,
        struct genome_data* genome
    )
{
    /* Assume the pseudo chr is sorted to come before the rest of the chrs */
    assert( PSEUDO_LOC_CHR_INDEX == 0 );

    struct pseudo_locations_t *ps_locs = genome->index->ps_locs;

    /* Find the start of the set of pseudo mapped locations with strand
     * matching the key we are matching to */
    int pslocs_start_in_sorted_mapped_locations =
        find_start_of_pseudo_mapped_locations_for_strand(
            potential_matches,
            key->strnd );

    /* If there aren't any pseudo locations with the same strand as the key,
     * nothing to do here */
    if( pslocs_start_in_sorted_mapped_locations == -1 )
        return;

    /* Binary search each pseudo location's locations for matches to the key */
    int i;
    for( i = pslocs_start_in_sorted_mapped_locations;
         i < potential_matches->length;
         i++ )
    {
        /* We begin at the start of a group of pseudo mapped_location's. If we
         * encounter a non-pseudo mapped_location, we've left the block of sorted
         * pseudo mapped_locations and are done. */
        if( potential_matches->locations[i].chr != PSEUDO_LOC_CHR_INDEX )
            break;

        /* binary search each pseudo location for matches to key */
        int ps_loc_index = potential_matches->locations[i].loc;
        struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;

        /* The pseudo locations are INDEX_LOC_TYPE, so we construct an
         * INDEX_LOC_TYPE to use as the key to bsearch.
         * NOTE only chr and loc are considered in cmp_genome_location. */
        INDEX_LOC_TYPE key_iloc;
        key_iloc.chr = key->chr;
        key_iloc.loc = key->loc;

        INDEX_LOC_TYPE* match = bsearch( &key_iloc,
                ps_loc->locs,
                ps_loc->num,
                sizeof(INDEX_LOC_TYPE),
                (int(*)(const void*, const void*))cmp_genome_location
            );

        if( match != NULL )
        {
            /* build a mapped_location with the metadata from key and the
             * location from the matching INDEX_LOC_TYPE */
            mapped_location tmp;
            copy_mapped_location( &tmp, key );
            tmp.chr = match->chr;
            tmp.loc = match->loc;

            add_mapped_location( &tmp, matching_subset );
        }
    }
}

void
find_matching_mapped_locations(
        mapped_locations* matching_subset,
        mapped_locations* potential_matches,
        
        mapped_location* base,
        int base_loc_subseq_offset,

        struct genome_data* genome
    )
{
    /* any base location should have been expanded into real locations */
    assert( base->chr != PSEUDO_LOC_CHR_INDEX );

    /*
     * construct a key (mapped_location) to search for
     *
     * this is subtle - our criteria for matching are
     * 1) strand
     * 2) chromosome
     * 3) location
     *
     * and the locations in base_locs and potential_matches both have distinct
     * subseq_offsets.
     *
     * If i and j are the mapped_locations we're comparing, we
     * want i - i_offset = j - j_offset for a match. We build a key mapped
     * location with the loc = i - i_offset + j_offset = j in order to compare
     * directory into potential_matches.
     */

    mapped_location key;
    copy_mapped_location( &key, base );
    key.loc =
        base->loc 
        - base_loc_subseq_offset
        + potential_matches->probe->subseq_offset;

    /* match to the pseudo locations */
    search_for_matches_in_pseudo_locations(
            matching_subset, potential_matches, &key, genome );

    /* match to the remaining locations */
    /* since base is not a pseudo location, we can simply do a binary
     * search over the potential_matches. */
    mapped_location* match = bsearch( &key,
            potential_matches->locations,
            potential_matches->length,
            sizeof(mapped_location),
            (int(*)(const void*, const void*))cmp_mapped_locations_by_location
        );

    // TODO ? handle potential for multiple matches
    if( match != NULL )
        add_mapped_location( match, matching_subset );
}

void
build_candidate_mappings_for_ungapped_assay_type(
        struct genome_data* genome,
        struct read_subtemplate* rst,
        mapped_locations_container* matches,
        candidate_mappings* rst_mappings,
        float min_match_penalty,
        enum assay_type_t assay )
{
    assert( matches->length > 0 );
    mapped_locations* locs = matches->container[0];

    int i;
    for( i = 0; i < locs->length; i++ )
    {
        mapped_location* loc = locs->locations;

        build_candidate_mapping_from_mapped_location(
                genome,
                rst,

                loc,
                locs,

                rst_mappings,
                min_match_penalty
            );
    }
}

void
build_candidate_mappings_for_gapped_assay_type(
        struct genome_data* genome,
        struct read_subtemplate* rst,
        mapped_locations_container* matches,
        candidate_mappings* rst_mappings,
        float min_match_penalty,
        enum assay_type_t assay )
{
    /* STUB */
    fprintf("FATAL       :  Gapped assay types are not implemented yet.\n");
    assert( false );
    exit(-1);
}

void
build_candidate_mappings_from_matched_mapped_locations(
        struct genome_data* genome,
        struct read_subtemplate* rst,
        mapped_locations_container* matches,
        candidate_mappings* rst_mappings,
        float min_match_penalty,
        enum assay_type_t assay )
{
    /* The way candidate mappings are built here depends on the assay.
     *
     * If it is a non-gapped assay, then we just want to build 1 candidate
     * mapping from each matched set of candidate mappings.
     *
     * If it is a gapped assay, then we want to build candidate mappings from
     * each mapped location in the matched set
     */

    if( assay == STRANDED_RNASEQ )
    {
        /* Gapped assay */
        build_candidate_mappings_for_gapped_assay_type(
                genome, rst, matches, rst_mappings,
                min_match_penalty, assay );
    } else {
        /* Non-gapped assay */
        build_candidate_mappings_for_ungapped_assay_type(
                genome, rst, matches, rst_mappings,
                min_match_penalty, assay );
    }
}

mapped_locations*
mapped_locations_template( struct indexable_subtemplate* indexable_probe )
{
    /* construct an empty set of mapped_locations with the metadata
     * (original indexable_subtemplate) from template */
    mapped_locations* locs = NULL;
    init_mapped_locations( &locs, indexable_probe );

    return locs;
}

void
add_pseudo_loc_to_mapped_locations(
        INDEX_LOC_TYPE* gen_loc,
        mapped_locations* results,
        mapped_location* loc
    )
{
    mapped_location tmp_loc;
    copy_mapped_location( &tmp_loc, loc );

    tmp_loc.chr = gen_loc->chr;
    tmp_loc.loc = gen_loc->loc;

    // add_and_expand (diploid) ?
    add_mapped_location( &tmp_loc, results );

    return;
}

void
expand_pseudo_location_into_mapped_locations(
        mapped_location* loc,
        mapped_locations* results,
        struct genome_data* genome
    )
{
    assert( loc->chr == PSEUDO_LOC_CHR_INDEX );

    struct pseudo_locations_t* ps_locs = genome->index->ps_locs;

    int ps_loc_index = loc->loc;
    struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;

    /* add every location to the results list */
    INDEX_LOC_TYPE* gen_locs = ps_loc->locs;
    int i;
    for( i = 0; i < ps_loc->num; i++ )
    {
        add_pseudo_loc_to_mapped_locations(
                &( gen_locs[i] ),
                results,
                loc
            );
    }

    return;
}

void
sort_mapped_locations_in_container(
        mapped_locations_container* mls_container
    )
{
    int i;
    for( i = 0; i < mls_container->length; i++ )
    {
        sort_mapped_locations_by_location( mls_container->container[i] );
    }

    return;
}

enum bool
matching_mapped_locations_containers(
        mapped_locations_container* search_results,
        mapped_locations_container* matches
    )
{
    /* all of the indexable subtemplates should have matches */
    if( search_results->length != matches->length )
        return false;

    /* each set of mapped_locations in the matches should have the same length.
     * There is a correspondence across the rows of the mapped_locations for a
     * match for a specific indexable_subtemplate */

    int i;
    assert( matches->length > 0 );
    for( i = 1; i < matches->length; i++ )
    {
        /* just compare everything mapped_location's length to the first
         * mapped_locations */
        if( matches->container[0]->length != matches->container[i]->length )
            return false;
    }

    return true;
}

void
expand_base_mapped_locations(
        mapped_location* base_loc,
        mapped_locations* expanded_locs,
        struct genome_data* genome
    )
{
    /* we assume that expanded_locs has already been initalize */
    
    /* if base_loc is a pseudo location, expand it and add all of its
     * potential mapped_location's */
    if( base_loc->chr == PSEUDO_LOC_CHR_INDEX )
    {
        expand_pseudo_location_into_mapped_locations(
                base_loc, expanded_locs, genome );
    } else {
        // add the mapped_location as-is
        add_mapped_location( base_loc, expanded_locs );
    }

    /* sort in order to use binary search later */
    sort_mapped_locations_by_location( expanded_locs );

    return;
}

void
build_candidate_mappings_from_base_mapped_location(
        mapped_location* base,
        int base_locs_index,
        struct indexable_subtemplate* base_probe,
        mapped_locations_container* search_results,
        struct genome_data* genome,
        struct read_subtemplate* rst,
        candidate_mappings* rst_mappings,
        float min_match_penalty,
        enum assay_type_t assay )
{
    /*
     * Takes in a base location, which necessarily corresponds to a single 
     * indexable subtemplate. Then, we try and match this base location with
     * a location in each other indexable subtemplate.
     *
     * base_locs stores the mapped_locations that come from the same indexable
     * subtemplate as base_loc. We need 'base_locs' container because:
          - we need to determine what indexable subtemplate base_loc came from, 
            so that we dont try to join it to mapped locations that came from
            the same indexable subtemplate
          - we need the meta data associated with the indexable subtemplate
     */

    /* initialize the container of matching mapped_locations */
    mapped_locations_container* matches = NULL;
    init_mapped_locations_container( &matches );
    
    /* build a mapped_locations container for the current base location */
    mapped_locations* base_container = mapped_locations_template( base_probe );
    add_mapped_location( base, base_container );
    add_mapped_locations_to_mapped_locations_container(
            base_container, matches );

    /* search the other mapped_locations for matches to base */
    int i;
    for( i = 0; i < search_results->length; i++ )
    {
        /* skip locations that came from the same indexable subtemplate 
           as base */
        if( i == base_locs_index )
            continue;

        mapped_locations* current_locs = search_results->container[i];
        mapped_locations* matching_subset =
            mapped_locations_template( current_locs->probe );

        find_matching_mapped_locations(
                matching_subset,
                current_locs,
                base, 
                base_probe->subseq_offset,
                genome
            );

        /* if we found matches, add the subset to the set of matches.
         * Otherwise, optimize by terminating early */
        if( matching_subset->length > 0 )
        {
            add_mapped_locations_to_mapped_locations_container(
                    matching_subset, matches );
        } else {
            free_mapped_locations( matching_subset );
            break;
        }
    }

    /* if we were able to match across all of the indexable subtemplates, this
     * is a valid match for base */
    if( matches->length == search_results->length )
    {
        build_candidate_mappings_from_matched_mapped_locations(
                genome,
                rst,
                matches,
                rst_mappings,
                min_match_penalty,
                assay
            );
    }
    
    free_mapped_locations_container( matches );
}

void
build_candidate_mappings_from_search_results(
        candidate_mappings* rst_mappings,
        mapped_locations_container* search_results,
        struct read_subtemplate* rst,
        struct genome_data* genome,
        float min_match_penalty,
        enum assay_type_t assay )
{
    /* sort each mapped_locations in order to binary search later */
    sort_mapped_locations_in_container( search_results );

    /* pick a mapped_locations to use as the basis for matching */
    int base_locs_index = choose_base_mapped_locations_index( search_results );
    mapped_locations* base_locs = search_results->container[base_locs_index];
    
    /* consider each base location */
    int i, j;
    for( i = 0; i < base_locs->length; i++ )
    {
        mapped_location* base_loc = base_locs->locations + i;

        /*** if the base_loc is a pseudo location, expand it to consider all
         *** the possible locations */
        
        /* initialize the expanded set of locations with the base's metadata */
        mapped_locations* expanded_base =
            mapped_locations_template( base_locs->probe );
        expand_base_mapped_locations( base_loc, expanded_base, genome );

        /* match across each of the expanded locations from the base */
        for( j = 0; j < expanded_base->length; j++ )
        {
            mapped_location* base = expanded_base->locations + j;

            build_candidate_mappings_from_base_mapped_location(
                    base,
                    base_locs_index,
                    base_locs->probe,
                    search_results,
                    genome,
                    rst,
                    rst_mappings,
                    min_match_penalty,
                    assay );
        }

        free_mapped_locations( expanded_base );
    }
}

void
find_candidate_mappings_for_read_subtemplate(
        struct read_subtemplate* rst,
        candidate_mappings* rst_mappings,

        struct genome_data* genome,
        struct error_model_t* error_model,
        struct error_data_record_t* scratch_error_data_record,

        float min_match_penalty,
        float max_penalty_spread,
        
        enum bool only_collect_error_data,

        enum assay_type_t assay
    )
{
    /* build the penalty arrays for this read subtemplate */
    struct penalty_array_t fwd_penalty_array, rev_penalty_array;
    init_penalty_array( &fwd_penalty_array, rst->length );
    init_penalty_array( &rev_penalty_array, rst->length );

    build_penalty_array( &fwd_penalty_array,
            error_model, rst->error_str );
    build_reverse_penalty_array( &fwd_penalty_array,
            &rev_penalty_array );

    // build a set of indexable subtemplates from this read subtemplate
    struct indexable_subtemplates* ists = NULL;
    init_indexable_subtemplates( &ists );
    build_indexable_subtemplates_from_read_subtemplate(
            ists, rst, genome->index,
            &fwd_penalty_array, &rev_penalty_array
        );

    // store all of the mapped_locations for this read subtemplate
    mapped_locations_container* ml_container = NULL;
    init_mapped_locations_container( &ml_container );

    search_index_for_indexable_subtemplates(
            ists,
            ml_container,

            genome,

            min_match_penalty,
            max_penalty_spread,

            only_collect_error_data
        );
    
    /* build candidate mappings from matching subsets of the mapped_locations
     * for each indexable_subtemplate returned by the index search */
    build_candidate_mappings_from_search_results(
            rst_mappings, ml_container, rst,
            genome, min_match_penalty, assay
        );

    /* Note - mapped_locations_container contains references to memory
     * allocated in the indexable_subtemplates, so they must always be freed
     * simultaneously */
    free_indexable_subtemplates( ists );
    free_mapped_locations_container( ml_container );

    /****** Do the recheck ******/
    recheck_locations(
            genome, 

            rst, rst_mappings,

            min_match_penalty,
            max_penalty_spread,

            &fwd_penalty_array,
            &rev_penalty_array
        );

    /* update the thread local copy of error data (need the error data
     * and the subtemplate to do this) */
    update_error_data_record_from_candidate_mappings(
            genome,
            rst_mappings,
            rst,
            scratch_error_data_record
        );

    /* cleanup memory */
    free_penalty_array( &fwd_penalty_array );
    free_penalty_array( &rev_penalty_array );
}

void
find_candidate_mappings_for_read(
        struct read* r,
        candidate_mappings* read_mappings,

        struct genome_data* genome,
        struct error_model_t* error_model,
        struct error_data_record_t* scratch_error_data_record,

        float min_match_penalty,
        float max_penalty_spread,

        enum bool only_collect_error_data,

        enum assay_type_t assay )
{
    // for each subtemplate in the read
    int rst_index;
    for( rst_index=0; rst_index < r->num_subtemplates; rst_index++ )
    {
        // reference to current read subtemplate
        struct read_subtemplate* rst = &( r->subtemplates[rst_index] );

        /* initialize the candidate mappings container for this read
         * subtemplate */
        candidate_mappings* rst_mappings = NULL;
        init_candidate_mappings( &rst_mappings );

        find_candidate_mappings_for_read_subtemplate(
                rst,
                rst_mappings,

                genome,
                error_model,
                scratch_error_data_record,

                min_match_penalty,
                max_penalty_spread,
                
                only_collect_error_data,

                assay
            );

        /* append the candidate mappings from this read subtemplate to the set
         * of candidate mappings for this read */
        append_candidate_mappings( read_mappings, rst_mappings );

        free_candidate_mappings( rst_mappings );
    }
}

void
add_mapped_reads_from_candidate_mappings(
        struct mapped_reads_db* mpd_rds_db,
        candidate_mappings* mappings,
        readkey_t readkey
    )
{
    mapped_read_t* mpd_rd = NULL;

    /* The paired end reads joining code in build_mapped_read_from_candidate_mappings
     * requires the candidate_mappings to be sorted */
    sort_candidate_mappings( mappings );

    build_mapped_read_from_candidate_mappings(
            &mpd_rd,
            mappings,
            readkey
        );

    /* if we were actually able to build a mapped read with > 0 prb */
    if( NULL != mpd_rd ) {
        add_read_to_mapped_reads_db( mpd_rds_db, mpd_rd );
        free_mapped_read( mpd_rd );
    }
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

    //int thread_id = td->thread_id;

    struct genome_data* genome = td->genome;
    
    unsigned int* mapped_cnt = td->mapped_cnt;
    pthread_spinlock_t* mapped_cnt_lock = td->mapped_cnt_lock;

    struct rawread_db_t* rdb = td->rdb;
    
    struct mapped_reads_db* mpd_rds_db = td->mpd_rds_db;

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
    enum bool only_collect_error_data = td->only_collect_error_data;
    
    /* END parameter 'recreation' */

    assert( genome->index != NULL );

    /* how often we print out the mapping status */
    #define MAPPING_STATUS_GRANULARITY 100000

    /* create a thread local copy of the error data to avoid excess locking */
    struct error_data_record_t* scratch_error_data_record;
    init_error_data_record( &scratch_error_data_record, 
                            error_data->max_read_length, 
                            error_data->max_qual_score );

    /* The current read of interest */
    readkey_t readkey;
    struct read* r;
    int curr_read_index = 0;
    /* 
     * While there are still mappable reads in the read DB. All locking is done
     * in the get next read functions. 
     */
    while( EOF != get_next_read_from_rawread_db(
               rdb, &readkey, &r, td->max_readkey )  
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
        
        // Make sure this read has "enough" HQ bps before trying to map it
        if( filter_read( r, error_model ) )
        {
            continue; // skip the unmappable read
        }
        
        /* Initialize container for candidate mappings for this read */
        candidate_mappings* mappings = NULL;
        init_candidate_mappings( &mappings );

        find_candidate_mappings_for_read(
                r,
                mappings,
                genome,
                error_model,

                scratch_error_data_record,

                min_match_penalty,
                max_penalty_spread,

                only_collect_error_data,

                rdb->assay
            );

        // Free the read
        free_read( r );

        // build mapped reads from the set of candidate mappings for this read
        if( !only_collect_error_data )
        {
            // count the number of valid candidate mappings in order to update
            // the mapped read count
            int num_valid_mappings = 0;
            int i;
            for( i = 0; i < mappings->length; i++ )
            {
                if( (mappings->mappings + i)->recheck == VALID )
                    num_valid_mappings += 1;
            }

            /* update the mapped_cnt */
            if( num_valid_mappings > 0 )
            {
                pthread_spin_lock( mapped_cnt_lock );
                *mapped_cnt += num_valid_mappings;
                pthread_spin_unlock( mapped_cnt_lock );
            }

            /* build mapped reads from the returned candidate mappings and add
             * to the mapped reads db */
            add_mapped_reads_from_candidate_mappings(
                    mpd_rds_db,
                    mappings,
                    readkey
                );
        }

        curr_read_index += 1;

        /* cleanup memory */
        free_candidate_mappings( mappings );
    }

    /********* update the global error data *********/

    // add error_data to scratch_error_data
    merge_in_error_data_record( error_data, -1, scratch_error_data_record );
    
    // free local copy of error data
    free_error_data_record( scratch_error_data_record );

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
                       struct mapped_reads_db* mpd_rds_db,

                       struct error_model_t* error_model,
                       struct error_data_t* error_data,
                       
                       float min_match_penalty, float max_penalty_spread )
{
    int rc;
    pthread_mutexattr_t mta;
    rc = pthread_mutexattr_init(&mta);
    assert( rc == 0 );
    
    td_template->genome = genome;
    
    td_template->mapped_cnt = malloc( sizeof(unsigned int) );
    *(td_template->mapped_cnt) = 0;
    
    td_template->mapped_cnt_lock = malloc( sizeof(pthread_spinlock_t) );
    rc = pthread_spin_init( td_template->mapped_cnt_lock,
                            PTHREAD_PROCESS_PRIVATE );
    assert( rc == 0 );

    td_template->max_readkey = 0;

    td_template->rdb = rdb;
    
    td_template->mpd_rds_db = mpd_rds_db;
    
    td_template->min_match_penalty = min_match_penalty;
    td_template->max_penalty_spread = max_penalty_spread;

    td_template->error_data = error_data;
    
    td_template->error_model = error_model;

    td_template->thread_id = 0;
    
    /* Make sure we are finding all mappers for a normal index search */
    td_template->only_collect_error_data = false;
    
    return;
}

void 
free_td_template( struct single_map_thread_data* td_template )
{
    free( td_template->mapped_cnt );

    pthread_spin_destroy( td_template->mapped_cnt_lock );
    free( (void*) td_template->mapped_cnt_lock );
}

/* bootstrap an already initialized error model */
void
bootstrap_estimated_error_model( 
    struct genome_data* genome,
    
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mpd_rds_db, // TODO set to NULL for bootstrap?
    
    struct error_model_t* error_model
) 
{
    assert( error_model != NULL );

    clock_t start = clock();
    
    struct error_model_t* bootstrap_error_model;
    init_error_model( &bootstrap_error_model, MISMATCH );

    struct error_data_t* error_data;
    init_error_data( &error_data );
    
    /* This is the fraction of bases that can have a mismatch */
    #define MAX_NUM_MM_RATE 0.20
    #define MAX_MM_SPREAD_RATE 0.10
    
    int max_num_mm = -(int)(MAX_NUM_MM_RATE*genome->index->seq_length) - 1;
    int max_mm_spread = (int)(MAX_MM_SPREAD_RATE*genome->index->seq_length) + 1;
    
    printf( "NOTICE      :  Setting bootstrap mismatch rates to %i and %i\n",
            max_num_mm, max_mm_spread );
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, 
                      bootstrap_error_model, error_data, 
                      max_num_mm, max_mm_spread );
    
    /* 
       only use unique mappers for the initial bootstrap. This is just a small
       performance optimization, it prevents us from going too deeply into the 
       index as soon as we know that a mapping isn't unique.
    */
    td_template.only_collect_error_data = true;

    // Detyermine how many reads we should look through for the bootstrap
    #define NUM_READS_TO_BOOTSTRAP READS_STAT_UPDATE_STEP_SIZE
    td_template.max_readkey = NUM_READS_TO_BOOTSTRAP;

    // Add a new row to store error data in
    add_new_error_data_record( error_data, 0, td_template.max_readkey );

    spawn_threads( &td_template );
        
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Bootstrapped (%i/%u) Unique Reads in %.2lf seconds ( %e/thread-hour )\n",
            *(td_template.mapped_cnt), rdb->readkey, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)*(td_template.mapped_cnt))*CLOCKS_PER_SEC*3600)/(stop-start)
        );

    free_td_template( &td_template );
    
    rewind_rawread_db( rdb );
    
    update_error_model_from_error_data( error_model, td_template.error_data );
    
    FILE* error_data_ofp = fopen( BOOTSTRAP_ERROR_STATS_LOG, "w" );
    log_error_data( error_data_ofp, error_data );
    fclose( error_data_ofp );
    free_error_data( error_data );
    
    free_error_model( bootstrap_error_model );

    return;
}

void
find_all_candidate_mappings(
        struct genome_data* genome,

        struct rawread_db_t* rdb,
        struct mapped_reads_db* mpd_rds_db,

        struct error_model_t* error_model,

        float min_match_penalty,
        float max_penalty_spread
    )
{
    clock_t start = clock();

    struct error_data_t* error_data;
    init_error_data( &error_data );
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, error_model,
                      error_data, min_match_penalty, max_penalty_spread );
    
    /* initialize the threads */
    while( false == rawread_db_is_empty( rdb ) )
    {
        // Add a new row to store error data in
        add_new_error_data_record( 
            td_template.error_data, 
            td_template.max_readkey, 
            td_template.max_readkey+READS_STAT_UPDATE_STEP_SIZE-1  );

        // update the maximum allowable readkey
        td_template.max_readkey += READS_STAT_UPDATE_STEP_SIZE;
        
        spawn_threads( &td_template );
        
        /* update the error model from the new error data */
        update_error_model_from_error_data(error_model, td_template.error_data);
    }
    
    /* Print out performance information */
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Mapped (%i/%u) Partial Reads in %.2lf seconds ( %e/thread-hour )\n",
            *(td_template.mapped_cnt), rdb->readkey, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)*(td_template.mapped_cnt))*CLOCKS_PER_SEC*3600)/(stop-start)
        );

    FILE* error_data_ofp = fopen( ERROR_STATS_LOG, "w" );
    log_error_data( error_data_ofp, error_data );
    fclose( error_data_ofp );
    free_error_data( error_data );
    
    free_td_template( &td_template );

    return;
}
