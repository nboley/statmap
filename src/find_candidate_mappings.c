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
        struct index_t* index, 
        struct indexable_subtemplate* ist,

        float min_match_penalty,
        float max_penalty_spread,
        mapped_locations** results,

        enum bool only_collect_error_data
    )
{
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

            /* length of the reads */
            subseq_length,

            /* the fwd stranded sequence */
            fwd_seq, 
            /* the bkwd stranded sequence */
            bkwd_seq, 

            ist->fwd_penalties,
            ist->rev_penalties,

            only_collect_error_data
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

static void
add_candidate_mapping_from_haploid (
        struct read_subtemplate* rst,
        mapped_location*    result,
        mapped_locations*   results,
        candidate_mapping   cm,
        candidate_mappings* mappings,
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
            (result->location).loc,
            (result->location).chr,
            result->strnd, 
            results->probe->subseq_offset,
            results->probe->subseq_length,
            rst->length,
            genome
        );
    }
    if( read_location < 0 ) // the read location was invalid; skip this mapped_location
        return;
    cm.start_bp = read_location;

    /* add the candidate mapping */
    add_candidate_mapping( mappings, &cm );
}

static void
add_candidate_mapping_from_diploid (
        struct read_subtemplate* rst,
        mapped_location*    result,
        mapped_locations*   results,
        candidate_mapping   cm,
        candidate_mappings* mappings,
        struct genome_data* genome
    )
{
    /* add the paternal candidate mapping. */
    add_candidate_mapping_from_haploid(
        rst, result, results, cm, mappings, genome
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
            maternal_start,
            maternal_chr_index,
            result->strnd, 
            results->probe->subseq_offset,
            results->probe->subseq_length,
            rst->length,
            genome
        );
    }
    if( read_location < 0 ) // the read location was invalid; skip this mapped_location
        return;

    /* modify cm to be maternal complement of original paternal cm */
    cm.start_bp = read_location;
    cm.chr = maternal_chr_index;

    /* add maternal candidate mapping */
    add_candidate_mapping( mappings, &cm );
}

static inline void
build_candidate_mappings_from_mapped_location(
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
    candidate_mapping template_candidate_mapping 
        = init_candidate_mapping_from_template( 
            rst, max_penalty_spread 
        );

    /* Make sure the "subseq" is acutally shorter than the read */
    assert( subseq_len <= rst->length );

    /* set the strand */
    assert( result->strnd == FWD || result->strnd == BKWD );
    template_candidate_mapping.rd_strnd = result->strnd;

    /* set metadata */
    template_candidate_mapping.penalty = result->penalty;
    template_candidate_mapping.subseq_offset = results->probe->subseq_offset;

    /* Build, verify, and add candidate mappings depending on type of loc */
    /* TODO might not need to pass the subtemplate to these fns; it is only
     * used for rst->length, and that is also stored on the cm template */
    if( result->location.is_paternal && result->location.is_maternal )
        add_candidate_mapping_from_diploid(
            rst, result, results, template_candidate_mapping, mappings, genome
        );
    else
        add_candidate_mapping_from_haploid(
            rst, result, results, template_candidate_mapping, mappings, genome
        );
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
update_error_data_from_candidate_mappings(
    struct genome_data* genome,
    candidate_mappings* mappings,
    struct read_subtemplate* rst,
    struct error_data_t* error_data
)
{
    if( !can_be_used_to_update_error_data( genome, mappings, rst ) )
        return;
    
    /* we need at least one valid mapping ( although this case should be handled
       above by can_be_used_to_update_error_data */
    if( mappings->length == 0 )
        return;
    
    // emphasize the array aspect with the + 0
    // but, since the length is exactly 1, we know that
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

    update_error_data( 
        error_data, genome_seq, read_seq, error_str, mapped_length );

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
                genome->index, 
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

mapped_locations*
choose_base_mapped_locations(
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

    return mls_container->container[min_index];
}

int
find_strand_pivot_in_sorted_mapped_locations(
        mapped_locations* sorted_mapped_locs
    )
{
    /* Since mapped locations are sorted by
     * 1) strand
     * 2) chromosome
     * 3) start bp
     *
     * we can find the index where sorted block of FWD strand mapped_locations
     * ends and the sorted block of BKWD strand mapped_locations begins
     */

    int pivot = -1; // If there are no BKWD strand reads, pivot will be -1

    int i;
    for( i = 0; i < sorted_mapped_locs->length; i++ )
    {
        if( sorted_mapped_locs->locations[i].strnd == BKWD )
        {
            pivot = i;
            break;
        }
    }

    return pivot;
}

void
search_for_matches_in_pseudo_locations(
        mapped_locations* matching_subset,
        mapped_locations* potential_matches,
        mapped_location* key,
        struct genome_data* genome
    )
{
    struct pseudo_locations_t *ps_locs = genome->index->ps_locs;

    /*
     * potential_matches is a sorted mapped_locations. Since mapped_locations
     * are sorted by strand, chromosome, and bp_start, and pseudo_locations are
     * identified by having their chromosome == PSEUDO_LOC_CHR_INDEX, there are
     * potentially two contiguous blocks of contiguous pseudo locations in
     * potential_matches.
     *
     * Furthermore, the specific key we are searching for has a strand. So we
     * want to identify which subset of pseudo locations could contain matches
     * for the given key, and then search inside of it.
     */

    int strand_pivot =
        find_strand_pivot_in_sorted_mapped_locations( potential_matches );

    /* Pseudo locations all have chromsome index PSEUDO_LOC_CHR_INDEX (0).
     * Assuming this does not change, there will be a block of all of the
     * pseudo locations at the beginning of each of the FWD and BKWD strand
     * mapped_locations segments */

    /* find the start and end of the region of potential_matches that contains
     * pseudo locations matching the strand of key */
    assert( FWD < BKWD );

    int pslocs_start, pslocs_end;

    if( key->strnd == FWD )
    {
        pslocs_start = 0; // FWD < BKWD
        if( strand_pivot > -1 )
        {
            /* if there were any BKWD stranded mapped_locations */
            pslocs_end = strand_pivot;
        } else {
            pslocs_end = potential_matches->length;
        }
    }
    else if( key->strnd == BKWD )
    {
        pslocs_end = potential_matches->length;

        if( strand_pivot > -1 )
        {
            pslocs_start = strand_pivot;
        } else {
            /* set start == end so we don't iterate over anything */
            pslocs_start = potential_matches->length;
        }
    } else {
        fprintf( stderr, "FATAL       :  Unrecognized strand %i on mapped_location.\n", key->strnd );
        assert(false);
        exit( -1 );
    }

    int i;
    for( i = pslocs_start; i < pslocs_end; i++ )
    {
        /* binary search each pseudo location for matches to key */
        int ps_loc_index = potential_matches->locations[i].location.loc;
        struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;

        /* The pseudo locations are just GENOME_LOC_TYPEs, so we need to
         * extract the GENOME_LOC_TYPE in order to use bsearch */
        GENOME_LOC_TYPE key_loc = key->location;

        GENOME_LOC_TYPE* match = bsearch( &key_loc,
                ps_loc->locs,
                ps_locs->num,
                sizeof(GENOME_LOC_TYPE),
                (int(*)(const void*, const void*))cmp_genome_location
            );

        // TODO ? handle potential for multiple matches
        if( match != NULL )
        {
            /* reconstruct mapped_location from key and the pseudo_location's
             * GENOME_LOC_TYPE */
            mapped_location tmp;
            copy_mapped_location( &tmp, key );
            tmp.location = *match;

            add_mapped_location( &tmp, matching_subset );
        }
    }
}

void
find_matching_mapped_locations(
        mapped_locations* matching_subset,
        mapped_locations* potential_matches,
        mapped_locations* base,

        struct genome_data* genome
    )
{
    int bm;
    /* for every location in the current location's matching subset */
    for( bm = 0; bm < base->length; bm++ )
    {
        // reference to the current base location we're considering
        mapped_location* base_loc = base->locations + bm;
        /* if any base locatino was a pseudo location, it should have been
         * expanded */
        assert( base_loc->location.chr != PSEUDO_LOC_CHR_INDEX );

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
        copy_mapped_location( &key, base_loc );
        key.location.loc =
            base_loc->location.loc 
            - base->probe->subseq_offset
            + potential_matches->probe->subseq_offset;

        /* match to the pseudo locations */
        search_for_matches_in_pseudo_locations(
                matching_subset, potential_matches, &key, genome );

        /* match to the remaining locations */
        /* since base_loc is not a pseudo location, we can simply do a binary
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
}

void
build_candidate_mappings_from_matched_mapped_locations(
        struct genome_data* genome,
        struct read_subtemplate* rst,
        mapped_locations_container* matches,
        candidate_mappings* rst_mappings,
        float min_match_penalty
    )
{
    /* for now, build a candidate mapping from the first mapped_location in
     * the first mapped_locations in the matches */

    /* for now, bulid a candidate mapping from all the mapped_locations
     * in the first set */
    assert( matches->length > 0 );
    mapped_locations* locs = matches->container[0];

    int i;
    for( i = 0; i < locs->length; i++ )
    {
        mapped_location* loc = &( locs->locations[i] );

        build_candidate_mappings_from_mapped_location(
                genome,
                rst,

                loc,
                locs,

                rst_mappings,
                min_match_penalty
            );
    }
}

mapped_locations*
mapped_locations_template(
        mapped_locations* template
    )
{
    /* construct an empty set of mapped_locations with the metadata from
     * template */
    mapped_locations* locs = NULL;
    init_mapped_locations( &locs, template->probe );

    return locs;
}

void
add_pseudo_loc_to_mapped_locations(
        GENOME_LOC_TYPE* gen_loc,
        mapped_locations* results,
        mapped_location* loc
    )
{
    mapped_location tmp_loc = *loc;
    tmp_loc.location = *gen_loc; // TODO check read start
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
    assert( loc->location.chr == PSEUDO_LOC_CHR_INDEX );

    struct pseudo_locations_t* ps_locs = genome->index->ps_locs;

    int ps_loc_index = loc->location.loc;
    struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;

    /* add every location to the results list */
    GENOME_LOC_TYPE* gen_locs = ps_loc->locs;
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

mapped_locations*
init_mapped_locations_container_for_matches(
        mapped_locations_container** matches,
        mapped_locations* base_locs,
        mapped_location* base_loc,

        struct genome_data* genome
    )
{
    init_mapped_locations_container( matches );

    /* initialize the set of matches with the base location's metadata */
    mapped_locations* matches_base =
        mapped_locations_template( base_locs );

    /* handle pseudo locations for the base location */
    /* if base_loc is a pseudo location, we expand it here and add all of its
     * potential mapped locations */
    if( base_loc->location.chr == PSEUDO_LOC_CHR_INDEX )
    {
        expand_pseudo_location_into_mapped_locations(
                base_loc, matches_base, genome );
    } else {
        // add the mapped_location as-is
        add_mapped_location( base_loc, matches_base );
    }

    // sort so we can use binary search when searching for matches_base
    sort_mapped_locations_by_location( matches_base );

    add_mapped_locations_to_mapped_locations_container(
            matches_base, *matches );

    return matches_base;
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

void
build_candidate_mappings_from_mapped_locations_container(
        candidate_mappings* rst_mappings,
        mapped_locations_container* search_results,
        struct read_subtemplate* rst,

        struct genome_data* genome,

        float min_match_penalty
    )
{
    /* sort so we can binary search in the mapped_locations later */
    sort_mapped_locations_in_container( search_results );

    /* pick a mapped_locations to use as the basis for matching */
    mapped_locations* base_locs = choose_base_mapped_locations( search_results );

    /* consider each location in base_locs as a candidate for building a set of
     * matching mapped locations */
    int i, j;
    for( i = 0; i < base_locs->length; i++ )
    {
        mapped_location* base_loc = base_locs->locations + i;
        
        /* store the matches as subsets of each set of mapped_locations in the
         * search_results mapped_locations_container */
        mapped_locations_container* matches = NULL;
        /* If base_loc is a pseudo location, it is expanded into base. O.w.,
         * base just contains base_loc */
        mapped_locations* base_for_matches = 
            init_mapped_locations_container_for_matches(
                    &matches, base_locs, base_loc, genome );
        
        /* search the other mapped_locations for matches */
        for( j = 0; j < search_results->length; j++ )
        {
            mapped_locations* current_locs = search_results->container[j];
            if( current_locs == base_locs )
                continue;

            /*
             * store the matching subset of mapped_locations from the current
             * base location and the current locations
             */
            mapped_locations* matching_subset = 
                mapped_locations_template( current_locs );

            find_matching_mapped_locations(
                    matching_subset,
                    current_locs,
                    base_for_matches,
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
        
        /* if we were able to match across all of the indexable subtemplates,
         * this is valid */
        if( matches->length == search_results->length )
        {
            build_candidate_mappings_from_matched_mapped_locations(
                    genome,
                    rst,
                    matches,
                    rst_mappings,
                    min_match_penalty
                );
        }
        
        free_mapped_locations_container( matches );
    }
}

void
find_candidate_mappings_for_read_subtemplate(
        struct read_subtemplate* rst,
        candidate_mappings* rst_mappings,

        struct genome_data* genome,
        struct error_model_t* error_model,
        struct error_data_t* thread_error_data,

        float min_match_penalty,
        float max_penalty_spread,
        
        enum bool only_collect_error_data
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

    /* build candidate mappings from each set of mapped locations in the mapped
     * locations container */
    build_candidate_mappings_from_mapped_locations_container(
            rst_mappings,
            ml_container,
            rst,

            genome,

            min_match_penalty
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
    update_error_data_from_candidate_mappings(
            genome,
            rst_mappings, rst,
            thread_error_data
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
        struct error_data_t* thread_error_data,

        float min_match_penalty,
        float max_penalty_spread,

        enum bool only_collect_error_data
    )
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
                thread_error_data,

                min_match_penalty,
                max_penalty_spread,

                only_collect_error_data
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
    struct mapped_read_t* mpd_rd;

    /* The paired end reads joining code in build_mapped_read_from_candidate_mappings
     * requires the candidate_mappings to be sorted */
    sort_candidate_mappings( mappings );

    build_mapped_read_from_candidate_mappings(
            &mpd_rd,
            mappings,
            readkey
        );

    /* if we were actual able to build a mapped read with > 0 prb */
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
    struct error_data_t* thread_error_data;
    init_error_data( &thread_error_data );

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

                thread_error_data,

                min_match_penalty,
                max_penalty_spread,

                only_collect_error_data
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

    // add thread_error_data to (global) error_data
    pthread_mutex_lock( error_data->mutex );
    add_error_data( error_data, thread_error_data );
    pthread_mutex_unlock( error_data->mutex );
    
    // free local copy of error data
    free_error_data( thread_error_data );
    
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
    
    struct error_model_t* bootstrap_error_model;
    init_error_model( &bootstrap_error_model, MISMATCH );
    
    /* put the search arguments into a structure */
    #define MAX_NUM_MM 5 // XXX should this be negative?
    #define MAX_MM_SPREAD 3
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, 
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

    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, error_model,
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
