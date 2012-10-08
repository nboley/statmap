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
        int bp = rst->char_seq[pos];
        if( bp == 'N' || bp == 'n' )
        {
            penalty += N_penalty;
        } else {
            int bpc = bp_code( bp );
            /* take the match penalty for this bp */
            penalty += penalty_array->array[pos].penalties[bpc][bpc];
        }
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
        struct index_search_params* search_params,
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
    replace_ns_inplace( sub_read, subseq_length );

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

            search_params,
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
    
    /* we need at least one valid mapping ( although this case should be handled
       above by can_be_used_to_update_error_data */
    if( mappings->length == 0 )
        return;
    
    // emphasize the array aspect with the + 0
    // but, since all of the sequence is identical,
    // we only need to deal with this read
    candidate_mapping* cm = mappings->mappings + 0;         
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

    update_error_data_record( 
        error_data_record, genome_seq, read_seq, error_str, mapped_length );

    /* free memory if we allocated it */
    if( cm->rd_strnd == BKWD )
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

    if( rst->length - softclip_len < subseq_length )
    {
        fprintf( stderr,
                 "FATAL       :  Probe lengths must be at least %i basepairs short, to account for the specified --soft-clip-length (-S)\n", softclip_len );
        assert( false );
        exit(1);
    }

    /* Finding the optimal subsequence offset is desirable for ungapped assays,
     * but is not best for gapped assays (where we want the ists to be as far
     * apart as possible to avoid overlapping with potential introns).
     *
     * For now, we manually set the ists to the first and last subsequences of
     * the read of length subseq_length, which will work well for both types of
     * assays. */

    int subseq_offset = 0;

    /* Our strategy for choosing indexable subtemplates depends on the type of
     * assay. For gapped assays (RNA_SEQ), we wish to maximize the distance
     * between probes in order to optimize intron finding. For ungapped assays,
     * we want to use the highest quality subsequence in the read for the index
     * search. */
    if( _assay_type == RNA_SEQ ) // Gapped assay
    {
        /* FIXME soft clipping for gapped assays? */
        if( range_start == 0 ) {
            /* subseq_offset = softclip_len ? will that mess up the matching
             * code later? */
            subseq_offset = 0;
        } else {
            subseq_offset = range_length - subseq_length;
        }
    } else if( _assay_type == CAGE ) {
        /* FIXME this will work for now but the value saved in config.data
         * (which is currently unused) will be wrong if we make the adjustment
         * here. */
        if( softclip_len < MAX_NUM_UNTEMPLATED_GS ) {
            fprintf( stderr,
                     "WARNING     :  Specified soft clip length (%i) for CAGE assay type is less than the maximum number of untemplated G's (%i)\n",
                     softclip_len, MAX_NUM_UNTEMPLATED_GS );
            fprintf( stderr,
                     "NOTICE      :  Setting soft clip length to maxmimum number of untemplated G's (minimum required for mapping CAGE)\n" );

            softclip_len = MAX_NUM_UNTEMPLATED_GS;
        }
    } else {
        subseq_offset = find_optimal_subseq_offset(
                rst,
                fwd_penalty_array,
                subseq_length,
                range_start,
                range_length
            );
    }

    /* The probe must be offset by at least the soft clip length; if there is
     * an optimal subseq offset greater than the soft clip length, that will be
     * used instead. */
    subseq_offset = MAX( subseq_offset, softclip_len );

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
        mapped_locations** search_results,

        struct genome_data* genome,

        struct index_search_params* index_search_params,

        enum bool only_collect_error_data
    )
{
    int ist_index;
    for( ist_index = 0; ist_index < ists->length; ist_index++ )
    {
        // reference to current indexable subtemplate
        struct indexable_subtemplate* ist = ists->container + ist_index;

        /**** go to the index for mapping locations */
        mapped_locations *results = NULL;

        search_index(
                genome,
                ist,
                index_search_params + ist_index,
                &results,
                only_collect_error_data
            );

        /* add the results for this indexable subtemplate to the array of
         * all the search results */
        search_results[ist_index] = results;
    }
}

int
choose_mapped_locations_base_index( int search_results_length )
{
    assert( search_results_length > 0 );

    /* TODO for now, use the first set of mapped locations, which corresponds
     * to the search results for the first (5'-most) indexable subtemplate */
    return 0;
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

/* TODO generic search, so we can use it for both the pseudo locations
 * and the regular locations */

/*
 * Our criteria for matching are
 * 1) strand
 * 2) chromosome
 * 3) location
 */

/* TODO this should consider strand, chr, and loc
 * can i drop in cmp_mapped_locations_by_location somehow?*/
int
search_for_matching_mapped_locations(
        int base_pos,
        int base_offset,
        mapped_locations* match_candidates )
{
    /* Search for a location with a start equal to base_pos - base_offset
     * (the true start of the base location ) */
    int base_start = base_pos - base_offset;

    /* Binary search */
    int low = 0;
    int high = match_candidates->length;
    while( low < high ) {
        int mid = low + ((high-low) / 2);

        /* Compute the true start of the location currently being considered
         * by the search */
        int current_loc_start = match_candidates->locations[mid].loc 
                              - match_candidates->probe->subseq_offset;

        if( current_loc_start < base_start ) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }

    /* make sure the binary search is working */
    assert( low <= high );
    assert( low <= match_candidates->length );
    assert( low >= 0 );

    /* This is the "deferred detection of equality" variant of the binary
     * search algorithm (http://en.wikipedia.org/wiki/Binary_search_algorithm).
     *
     * It has the useful property that if there are multiple matching keys, it
     * returns the index of the first matching key in the region of matching
     * keys. */

    /* Check for a match */
    int found_loc_start = match_candidates->locations[low].loc 
                        - match_candidates->probe->subseq_offset;
    if( (low == high) && found_loc_start == base_start )
    {
        /* the index of the matching location in match_candidates */
        return low;
    }

    /* Otherwise, we did not find a perfect match. However, we now have a
     * idea of the "closest" location in the match_candidates from the
     * failed binary search.
     *
     * Starting with low, do a linear search until we get to the first value
     * that is greater than the base_start (it could be within the range of
     * base_start + max_reference_insert_len) */

    int i;
    for( i = low; i < match_candidates->length; i++ )
    {
        int current_loc_start = match_candidates->locations[low].loc 
                              - match_candidates->probe->subseq_offset;

        if( current_loc_start > base_start )
        {
            /* the index of the first "greater" location in match_candidates */
            return i;
        }
    }

    /* Otherwise, no potential matches were found for base. Return -1 to indicate
     * no index of a potential match */
    return -1;
}

void
build_candidate_mapping_cigar_string_from_match(
        candidate_mapping* cm,
        struct ml_match* match,
        struct read_subtemplate* rst )
{
    /* Index of the current entry in the cigar string to update. This is
     * incremented every time there is a sequence with a reference gap */
    int cigar_index = 0;

    int i;
    for( i = 0; i < match->len; i++ )
    {
        /* If there is a gap in the reference, represent it with an N op. 
         * We check the value of i because obviously there cannot be an intron
         * before the first mapped location.*/
        if( i > 0 )
        {
            /* The reference gap is computed differently based on whether the
             * candidate mapping is fwd or rev stranded */
            int ref_gap;

            if( cm->rd_strnd == FWD )
            {
                ref_gap = (match->locations[i]->loc - match->subseq_offsets[i])
                        - (match->locations[i-1]->loc - match->subseq_offsets[i-1]);
            } else {
                ref_gap = (match->locations[i-1]->loc + match->subseq_offsets[i-1])
                        - (match->locations[i]->loc + match->subseq_offsets[i]);
            }

            if( ref_gap > 0 )
            {
                /* If there is a reference gap, add an N op with the known
                 * length of the intron, and then add two surrounding markers
                 * to show that we don't (yet) know the final length of the
                 * adjacent M ops. */
                cm->cigar[cigar_index+1].op = 'U';
                cm->cigar[cigar_index+1].len = 0;

                cm->cigar[cigar_index+2].op = 'N';
                cm->cigar[cigar_index+2].len = ref_gap;

                cm->cigar[cigar_index+3].op = 'U';
                cm->cigar[cigar_index+3].len = 0;
                cigar_index += 4;
            }
        }

        /* For each location, add an M op for the region of aligned sequence */
        cm->cigar[cigar_index].op = 'M';
        /* if this is the first indexable subtemplate, we assume that this
         * matches from the beginning of the read to the start of the indexable
         * subtemplate, so we add the subsequence offset for the first ist. */
        if( i == 0 ) {
            cm->cigar[cigar_index].len += match->subseq_offsets[i];
        }

        /* cigar[cigar_index].len was initialized to zero in 
         * init_candidate_mapping_from_read_subtemplate */
        cm->cigar[cigar_index].len += match->subseq_lengths[i];
    }

    /* Finally, we assume the last indexable subtemplate in the match matches
     * to the end of the read subtemplate. We add the additional length between
     * the end of the last indexable subtemplate and the length of the entire
     * read */
    cm->cigar[cigar_index].len += rst->length
        - (match->subseq_offsets[i-1] + match->subseq_lengths[i-1]);

    cm->cigar_len = cigar_index + 1;
    assert( cm->cigar_len <= MAX_CIGAR_STRING_ENTRIES );

    /* If this wasn't a gapped candidate mapping, hard-code the full read
     * length. TODO hack */
    if( cm->cigar_len == 1 )
    {
        cm->cigar[0].op = 'M';
        cm->cigar[0].len = cm->mapped_length;
    }

    return;
}

float
compute_candidate_mapping_penalty_from_match(
        struct ml_match* match )
{
    float cum_penalty = 0;

    int i;
    for( i = 0; i < match->len; i++ ) {
        /* the product of the marginal (log) probabilites */
        cum_penalty += match->locations[i]->penalty;
    }

    return cum_penalty;
}

void
build_candidate_mapping_from_match(
        struct ml_match* match,
        candidate_mappings* mappings,
        struct read_subtemplate* rst,
        struct genome_data* genome )
{
    /* build a candidate mapping from the base location, then build a cigar
     * string representing all the locations in the match */

    candidate_mapping cm = init_candidate_mapping();

    assert( match->len > 0 );
    mapped_location* base = match->locations[0];

    /* set the strand */
    assert( base->strnd == FWD || base->strnd == BKWD ); // XXX correct?
    cm.rd_strnd = base->strnd;

    /* set metadata */
    cm.penalty = compute_candidate_mapping_penalty_from_match( match );

    /* set location information */
    cm.chr = base->chr;

    /* set (partial) READ_TYPE information
     * (rd_type.pos is set in update_read_type_pos) */
    /* TODO this is unnecessary if we follow through with this
     * reworking on candidate_mapping */
    //cm.rd_type.follows_ref_gap = follows_ref_gap;

    /* We need to play with this a bit to account for index probes that are
     * shorter than the read length */
    cm.mapped_length = rst->length - softclip_len;
    int read_location = modify_mapped_read_location_for_index_probe_offset(
                base->loc, base->chr, base->strnd,
                match->subseq_offsets[0], match->subseq_lengths[0],
                cm.mapped_length, genome );

    /* HACK */
    if( match->len > 1 && base->strnd == BKWD ) {
        read_location = match->locations[match->len-1]->loc;
    }

    if( read_location < 0 ) // the read location was invalid; skip this matched set
        return;

    cm.start_bp = read_location;

    /* build the cigar string from the mapped_locations */
    build_candidate_mapping_cigar_string_from_match( &cm, match, rst );

    //print_candidate_mapping( &cm ); // DEBUG

    add_candidate_mapping( mappings, &cm );
}

void
condense_candidate_mapping_cigar_string(
        candidate_mapping* cm )
{
    /* Merge the final U lengths into their adjacent match counterparts and
     * remove the U ops from the cigar string */

    /* TODO for now, this is hard-coded */
    assert( cm->cigar_len == 5 );
    assert( cm->cigar[1].op == 'U' && cm->cigar[3].op == 'U' );

    /* TODO - remove all code related to 'U' - since we now examine all
     * possibilities within the index probes as well, it is unnecessary and
     * complicates things */
    int first_match_len = cm->cigar[1].len;
    int gap_len = cm->cigar[2].len;
    int second_match_len = cm->cigar[3].len;

    /* Clear the original cigar string array */
    memset( cm->cigar, 0,
            sizeof( struct CIGAR_ENTRY ) * MAX_CIGAR_STRING_ENTRIES );

    /* Build the new, condensed cigar string array */
    cm->cigar[0].op = 'M';
    cm->cigar[0].len = first_match_len;

    cm->cigar[1].op = 'N';
    cm->cigar[1].len = gap_len;

    cm->cigar[2].op = 'M';
    cm->cigar[2].len = second_match_len;

    cm->cigar_len = 3;

    return;
}

void
build_candidate_mappings_from_matches(
        struct genome_data* genome,
        struct read_subtemplate* rst,
        struct ml_matches* matches,
        candidate_mappings* rst_mappings )
{
    int i;
    for( i = 0; i < matches->length; i++ )
    {
        struct ml_match* match = matches->matches[i];

        build_candidate_mapping_from_match(
                match, rst_mappings, rst, genome );
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
add_matches_from_pseudo_locations_to_stack(
        mapped_locations* candidate_locs,
        struct ml_match* match,
        int match_index,
        int prev_matched_location_start,
        struct ml_match_stack* stack,
        struct genome_data* genome,
        struct mapping_params* mapping_params )
{
    /* Get a reference to the base mapped location for the strand */
    mapped_location* base = match->locations[0];

    /* Assume the pseudo chr is sorted to come before the rest of the chrs */
    assert( PSEUDO_LOC_CHR_INDEX == 0 );

    struct pseudo_locations_t *ps_locs = genome->index->ps_locs;

    /* Find the start of the set of pseudo mapped locations with the same
     * strand as the location we are matching to */
    int pslocs_start_in_sorted_mapped_locations =
        find_start_of_pseudo_mapped_locations_for_strand(
            candidate_locs, base->strnd );

    /* If there are no pseudo locations in the set of candidate locations for
     * matching, nothing to do here. */
    if( pslocs_start_in_sorted_mapped_locations < 0 )
        return;

    /* Iterate over the sets of pseduo locations */
    int i;
    for( i = pslocs_start_in_sorted_mapped_locations;
         i < candidate_locs->length;
         i++ )
    {
        mapped_location* candidate_loc = candidate_locs->locations + i;

        /* Once we're in the right section of strand, the locations are sorted
         * by chromosome. Since the pseudo chromosome is sorted to come before
         * the rest of the chromosomes, we know that if the current location's
         * chromosome is not the pseudo chromosome, there aren't any more
         * pseudo locations and we are done. */
        if( candidate_loc->chr > PSEUDO_LOC_CHR_INDEX )
            break;

        /* check the cumulative penalty for the match. If it is less
         * than the minimum, skip this match */
        float cum_penalty = match->cum_penalty + candidate_loc->penalty;
        if( cum_penalty < mapping_params->recheck_min_match_penalty )
            continue;

        int ps_loc_index = candidate_loc->loc;
        struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;

        /* Consider each of the pseudo locations. This is a linear search, but
         * the pseudo locations are sorted by location so we could use a binary
         * search to make it faster */
        int j;
        for( j = 0; j < ps_loc->num; j++ )
        {
            /* get the location */
            INDEX_LOC_TYPE* iloc = ps_loc->locs + j;

            /* check that the chromosome matches (strand already checked) */
            if( iloc->chr != base->chr )
                continue;

            /* Get the gap between this candidate for the match and the rest of the
             * partial match (depends on strand) */
            int candidate_ref_gap;
            if( base->strnd == FWD )
            {
                candidate_ref_gap = iloc->loc
                    - candidate_locs->probe->subseq_offset
                    - prev_matched_location_start;
            } else {
                candidate_ref_gap = iloc->loc
                    + candidate_locs->probe->subseq_offset
                    - prev_matched_location_start;
            }

            int candidate_cum_ref_gap = match->cum_ref_gap + candidate_ref_gap;

            if( candidate_cum_ref_gap < 0 )
                continue;

            if( candidate_cum_ref_gap <= max_reference_insert_len )
            {
                /* construct a mapped location */
                mapped_location* tmp = malloc( sizeof( mapped_location ));
                copy_mapped_location( tmp, candidate_loc );
                tmp->chr = iloc->chr;
                tmp->loc = iloc->loc;

                /* Since we're building this mapped_location dynamically to
                 * represent a pseudo location, it won't be free with the rest
                 * of the mapped_locations pointers from search_results.
                 * For now, we set an explicit flag to distinguish these
                 * different types of mapped locations. */
                tmp->free_with_match = true;

                struct ml_match* new_match = copy_ml_match( match );
                add_location_to_ml_match( tmp,
                        new_match,
                        candidate_locs->probe->subseq_length,
                        candidate_locs->probe->subseq_offset,
                        match_index,
                        candidate_cum_ref_gap );

                ml_match_stack_push( stack, new_match );
            } else {
                break;
            }
        }
    }

    return;
}

void
add_matches_from_locations_to_stack(
        mapped_locations* candidate_locs,
        struct ml_match* match,
        int match_index,
        int prev_matched_location_start,
        struct ml_match_stack* stack,
        struct mapping_params* mapping_params )
{
    mapped_location* base = match->locations[0];

    /* try to continue building match with each location in the
     * candidate_locs. Since each set of mapped_locations is sorted by start
     * location, we can terminate as soon as matching criteria fails (every
     * following attempt will also fail) */
    int i;
    for( i = 0; i < candidate_locs->length; i++ )
    {
        mapped_location* candidate_loc = candidate_locs->locations + i;

        if( candidate_loc->strnd != base->strnd ||
            candidate_loc->chr != base->chr )
            continue;

        /* Get the gap between this candidate for the match and the rest of the
         * partial match (depends on strand) */
        int candidate_ref_gap;
        if( base->strnd == FWD )
        {
            candidate_ref_gap = candidate_loc->loc
                - candidate_locs->probe->subseq_offset
                - prev_matched_location_start;
        } else {
            candidate_ref_gap = prev_matched_location_start
                - ( candidate_loc->loc + candidate_locs->probe->subseq_offset );
        }

        int candidate_cum_ref_gap = match->cum_ref_gap + candidate_ref_gap;

        if( candidate_cum_ref_gap < 0 )
            continue;

        if( candidate_cum_ref_gap <= max_reference_insert_len )
        {
            /* check the cumulative penalty for the match. If it is less
             * than the minimum, skip this match */
            float cum_penalty = match->cum_penalty + candidate_loc->penalty;
            if( cum_penalty < mapping_params->recheck_min_match_penalty ) {
                continue;
            }

            /* build a new match with the current location, and push it onto
             * the stack */
            struct ml_match* new_match = copy_ml_match( match );

            add_location_to_ml_match( candidate_loc,
                    new_match,
                    candidate_locs->probe->subseq_length,
                    candidate_locs->probe->subseq_offset,
                    match_index,
                    candidate_cum_ref_gap );

            ml_match_stack_push( stack, new_match );
        } else {
            /* adding this mapped_location to the match would cause it to have
             * total intron length greater than the maximum accepted size.
             * Since the mapped_locations in search_results are sorted by
             * start location, we then know that none of the following
             * mapped_locations can build a valid match, and we can optimize
             * by terminating early. */
            break;
        }
    }

    return;
}

int
find_index_of_next_indexable_subtemplate_to_match(
        struct ml_match* match )
{
    int match_index = 0;
    /*
     * find the index of the next set of index subtemplate search results
     * to consider. Since match is a partially built set of matched locations,
     * the index is the first entry in it's mapped_locations array that is NULL.
     */
    int i;
    for( i = 0; i < match->len; i++ )
    {
        if( match->locations[i] == NULL )
        {
            match_index = i;
            break;
        }
    }
    /* Since we initialize match with the base matched location, this index
     * should always be greater than zero */
    assert( match_index > 0 );

    return match_index;
}

void
add_potential_matches_to_stack(
        struct ml_match* match,
        struct ml_match_stack* stack,
        mapped_locations** search_results,
        struct genome_data* genome,
        struct mapping_params* mapping_params )
{
    int match_index =
        find_index_of_next_indexable_subtemplate_to_match( match );

    mapped_locations* candidate_locs = search_results[match_index];

    /* Compute the true start (loc - subseq_offset) of the last mapped
     * location in the match */
    int prev_matched_location_start;
    mapped_location* base = match->locations[0];
    if( base->strnd == FWD )
    {
        prev_matched_location_start = match->locations[match_index - 1]->loc
                                    - match->subseq_offsets[match_index - 1];
    } else {
        prev_matched_location_start = match->locations[match_index - 1]->loc
                                    + match->subseq_offsets[match_index - 1];
    }

    add_matches_from_pseudo_locations_to_stack(
            candidate_locs,
            match,
            match_index,
            prev_matched_location_start,
            stack,
            genome,
            mapping_params
        );

    add_matches_from_locations_to_stack(
            candidate_locs,
            match,
            match_index,
            prev_matched_location_start,
            stack,
            mapping_params
        );

    return;
}

void
find_matching_mapped_locations(
        struct ml_match* base_match, 
        struct ml_matches* matches,
        mapped_locations** search_results,
        struct genome_data* genome,
        struct mapping_params* mapping_params )
{
    /* Initialize the stack of partially completed matches */
    struct ml_match_stack* stack = NULL;
    init_ml_match_stack( &stack );

    /* Push the base_match on to set up the algorithm */
    ml_match_stack_push( stack, base_match );

    while( !ml_match_stack_is_empty(stack) )
    {
        /* pop a partially completed match object */
        struct ml_match* curr_match = ml_match_stack_pop( stack );

        /* if this is a completed, valid match object (i.e. it contains a
         * matched mapped location from each indexable subtemplate */
        if( ml_match_is_valid( curr_match ) )
        {
            /* then add it to the list of found matches */
            copy_ml_match_into_matches( curr_match, matches );
        } else {
            /* if this match object is incomplete, build partial matches using
             * the next set of indexable subtemplates, and add the new
             * potential matches to the stack */
            add_potential_matches_to_stack(
                    curr_match,
                    stack,
                    search_results,
                    genome,
                    mapping_params
                );
        }

        free_ml_match( curr_match, false );
    }

    free_ml_match_stack( stack );
}

void
build_candidate_mappings_from_base_mapped_location(
        mapped_location* base,
        int base_locs_index,
        struct indexable_subtemplate* base_probe,

        mapped_locations** search_results,
        int search_results_length,

        struct genome_data* genome,
        struct read_subtemplate* rst,
        candidate_mappings* rst_mappings,
        
        struct mapping_params* mapping_params )
{
    /* For now, we assume that we always start building from 5' -> 3' */
    assert( base_locs_index == 0 );

    /* Initialize the match object for all matches from this base location */
    struct ml_match* base_match = NULL;
    init_ml_match( &base_match, search_results_length );
    add_location_to_ml_match( base, base_match,
            base_probe->subseq_length, base_probe->subseq_offset,
            base_locs_index, 0 );

    /* Initialize container for complete, valid matches to this base location */
    struct ml_matches* matches = NULL;
    init_ml_matches( &matches );

    find_matching_mapped_locations(
            base_match, 
            matches,
            search_results,
            genome,
            mapping_params );

    build_candidate_mappings_from_matches(
            genome,
            rst,
            matches,
            rst_mappings
        );
    
    /* cleanup memory */
    free_ml_matches( matches, true );
}

void
sort_search_results(
        mapped_locations** search_results,
        int search_results_length )
{
    int i;
    for(i = 0; i < search_results_length; i++ )
    {
        sort_mapped_locations_by_location( search_results[i] );
    }

    return;
}

void
build_candidate_mappings_from_search_results(
        candidate_mappings* rst_mappings,
        mapped_locations** search_results,
        int search_results_length,
        struct read_subtemplate* rst,
        struct genome_data* genome,
        struct mapping_params* mapping_params )
{
    /* sort each mapped_locations in order to binary search later */
    sort_search_results( search_results, search_results_length );

    /* pick a mapped_locations to use as the basis for matching */
    int base_locs_index = choose_mapped_locations_base_index(
            search_results_length );
    mapped_locations* base_locs = search_results[base_locs_index];
    
    /* consider each base location */
    int i, j;
    for( i = 0; i < base_locs->length; i++ )
    {
        mapped_location* base_loc = base_locs->locations + i;

        /* If the base_loc is a diploid or pseudo location, build a set of all
         * possible expansions to consider for matching */
        mapped_locations* expanded_base =
            mapped_locations_template( base_locs->probe );

        expand_base_mapped_locations( base_loc, expanded_base, genome );

        /* match across each of the expanded locations */
        for( j = 0; j < expanded_base->length; j++ )
        {
            mapped_location* base = expanded_base->locations + j;

            build_candidate_mappings_from_base_mapped_location(
                    base, base_locs_index, base_locs->probe,
                    search_results, search_results_length,
                    genome,
                    rst, rst_mappings,
                    mapping_params
                );
        }

        free_mapped_locations( expanded_base );
    }
}

void
free_search_results(
        mapped_locations** search_results,
        int search_results_length
    )
{
    /* free each of the mapped_locations stored in the array */
    int i;
    for( i = 0; i < search_results_length; i++ )
    {
        free_mapped_locations( search_results[i] );
    }

    /* free the array of pointers */
    free( search_results );

    return;
}

int
count_gaps_in_candidate_mapping(
        candidate_mapping* cm )
{
    /* Count the number of gaps in a candidate mapping. Each gap is represented
     * by the same sequence of characters in the cigar string -
     * MUNUM. We can therefore count the number of M's - 1 to get the number of
     * gaps. */
    int num_m = 0;

    int i;
    for( i = 0; i < cm->cigar_len; i++ )
    {
        if( cm->cigar[i].op == 'M' )
            num_m++;
    }

    assert( num_m > 0 );
    return num_m - 1;
}

void
build_gapped_candidate_mappings_for_candidate_mapping(
        candidate_mapping* cm,
        candidate_mappings* gapped_mappings,
        int num_intron_configurations,
        struct read_subtemplate* rst,
        struct genome_data* genome,
        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa )
{
    /* Build the initial proposed gapped candidate mapping and compute its
     * penalty (all of the remaining sequence is appened to the first index
     * probe) */
    candidate_mapping gapped_cm = *cm;

    /* Update the tmp mapping's cigar string with the proposed probe lengths */
    assert( gapped_cm.cigar_len == 5 );
    assert( gapped_cm.cigar[1].op == 'U' && gapped_cm.cigar[3].op == 'U' );
    gapped_cm.cigar[1].len = num_intron_configurations;
    gapped_cm.cigar[3].len = 0;

    /* Build a copy of the sequence that was actually mapped to use in
     * comparing to the genome */
    char* mapped_seq = calloc( rst->length+1, sizeof(char) );
    assert( mapped_seq != NULL );
    struct penalty_array_t* penalty_array;

    if( gapped_cm.rd_strnd == FWD ) {
        memcpy( mapped_seq, rst->char_seq, sizeof(char)*rst->length );
        penalty_array = fwd_pa;
    } else {
        rev_complement_read( rst->char_seq, mapped_seq, rst->length );
        penalty_array = rev_pa;
    }

    /* The full length of the fragment, including any gaps */
    int gapped_length = gapped_cm.mapped_length + gapped_cm.cigar[2].len;

    /* Get a pointer to the genome sequence this read maps to. This is the
     * whole genome sequence it covers, including the intron. */
    char* genome_seq = find_seq_ptr( genome,
            gapped_cm.chr,
            gapped_cm.start_bp,
            gapped_length
        );

    /* Start by computing the penalty of a candidate mapping where the exon is
     * the first num_intron_configurations bp's in the read.
     * 
     * This corresponds to 
     * 
     *          rst->length
     *              |
     *      -----------------
     *      |  p1  |   e1   |    intron     |  p2  |
     *      -------
     *         |
     *      indexed_seq_len
     *
     *  where e2 is not shown because it has zero length.
     *
     * */
    float exon1_penalty = recheck_penalty(
            genome_seq,
            mapped_seq,
            penalty_array->array,
            gapped_cm.cigar[1].len
        );

    float exon2_penalty = 0;

    /* Add the initial proposed candidate mapping */
    /* TODO could optimize here by not adding anything below the minimum
     * match penalty. could potentially optimize further by breaking once
     * penalties start climbing again (assumes approximately normal
     * distribution around the intron). */
    gapped_cm.penalty = cm->penalty + exon1_penalty + exon2_penalty;
    add_candidate_mapping( gapped_mappings, &gapped_cm );

    int i;
    for( i = 1; i < num_intron_configurations; i++ )
    {
        gapped_cm.cigar[1].len -= 1;
        gapped_cm.cigar[3].len += 1;

        /* Compute the penalty of the shifted bps */

        /* We are subtracting one base pair from the end of the first exon, and
         * adding one to the start of the second exon */
        float exon1_bp_penalty = recheck_penalty(
                genome_seq + gapped_cm.cigar[1].len,
                mapped_seq + gapped_cm.cigar[1].len,
                penalty_array->array + gapped_cm.cigar[1].len,
                1
            );

        float exon2_bp_penalty = recheck_penalty(
                genome_seq + gapped_length - gapped_cm.cigar[3].len,
                mapped_seq + num_intron_configurations - gapped_cm.cigar[3].len,
                penalty_array->array + num_intron_configurations - gapped_cm.cigar[3].len,
                1
            );

        exon1_penalty -= exon1_bp_penalty;
        exon2_penalty += exon2_bp_penalty;

        /* Update the cm penalty (sum of log probs is product of probs)
         * and add the candidate mapping to the list of gapped mappings */
        gapped_cm.penalty = cm->penalty + exon1_penalty + exon2_penalty;
        add_candidate_mapping( gapped_mappings, &gapped_cm );
    }

    free( mapped_seq );

    return;
}

enum bool
find_potential_gapped_candidate_mappings(
        candidate_mapping* cm,
        candidate_mappings* gapped_mappings,
        struct read_subtemplate* rst,
        struct indexable_subtemplates* ists,
        struct genome_data* genome,
        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa,
        struct mapping_params* mapping_params )
{
    int num_gaps = count_gaps_in_candidate_mapping( cm );

    /* If there is no gap in this candidate mapping, then it must not have
     * spanned an intron. Skip it. */
    if( num_gaps == 0 )
        return false;

    /* for now, only allow one gap (one intron in the read subtemplate) */
    assert( num_gaps == 1 );
    assert( ists->length == 2 );

    int num_intron_configurations = rst->length;

    /* Allocate memory to store the potential gapped candidate mappings */
    candidate_mappings* potential_gapped_mappings;
    init_candidate_mappings( &potential_gapped_mappings );

    build_gapped_candidate_mappings_for_candidate_mapping(
            cm,
            potential_gapped_mappings,
            num_intron_configurations,
            rst,
            genome,
            fwd_pa, rev_pa
        );

    // DEBUG
    // print out the penalties of the potential gapped mappings so we can see
    // if we're getting the kind of distribution we expect
    /*
    fprintf(stderr, "===============================================\n");
    int pi;
    for( pi = 0; pi < potential_gapped_mappings->length; pi++ )
    {
        fprintf( stderr, "%f\n", potential_gapped_mappings->mappings[pi].penalty );
    }
    */

    /* Unpack the search parameters */
    float min_match_penalty = mapping_params->recheck_min_match_penalty;

    /* If a potential gapped mapping has penalty above the minimum, add it to
     * the final set of gapped mappings. */
    int i;
    for( i = 0; i < potential_gapped_mappings->length; i++ )
    {
        candidate_mapping* current_mapping
            = potential_gapped_mappings->mappings + i;

        if( current_mapping->penalty >= min_match_penalty )
        {
            condense_candidate_mapping_cigar_string( current_mapping );
            add_candidate_mapping( gapped_mappings, current_mapping );
        }
    }

    free_candidate_mappings( potential_gapped_mappings );

    return true;
}

void
build_gapped_candidate_mappings_for_read_subtemplate(
        candidate_mappings* mappings,
        candidate_mappings* gapped_mappings,
        struct read_subtemplate* rst,
        struct indexable_subtemplates* ists,
        struct genome_data* genome,
        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa,
        struct mapping_params* mapping_params )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        candidate_mapping* cm = mappings->mappings + i;

        enum bool mapping_was_gapped
            = find_potential_gapped_candidate_mappings(
                cm,
                gapped_mappings,
                rst,
                ists,
                genome,
                fwd_pa, rev_pa,
                mapping_params
            );

        if( !mapping_was_gapped )
        {
            /* Then this candidate mapping did not have a gap, and should be
             * added to the set of mappings to return as-is */
            add_candidate_mapping( gapped_mappings, cm );
        }
    }

    return;
}

void
find_candidate_mappings_for_read_subtemplate(
        struct read_subtemplate* rst,
        candidate_mappings* final_mappings,
        struct genome_data* genome,
        
        struct error_data_record_t* scratch_error_data_record,
        enum bool only_collect_error_data,

        struct mapping_params* mapping_params
    )
{
    int rst_pos = rst->pos_in_template.pos;
    
    // build a set of indexable subtemplates from this read subtemplate
    struct indexable_subtemplates* ists = NULL;
    init_indexable_subtemplates( &ists );
    build_indexable_subtemplates_from_read_subtemplate(
            ists, rst, genome->index,
            mapping_params->fwd_penalty_arrays[rst_pos],
            mapping_params->rev_penalty_arrays[rst_pos]
        );
    
    /* Stores the results of the index search for each indexable subtemplate */
    int search_results_length = ists->length;
    mapped_locations** search_results = malloc(
            sizeof(mapped_locations*) * search_results_length );

    /* initialize search parameters for the index probes */
    struct index_search_params* index_search_params = NULL;
    init_index_search_params( &index_search_params, ists, mapping_params );

    search_index_for_indexable_subtemplates(
            ists,
            search_results,
            genome,
            index_search_params,
            only_collect_error_data
        );

    free( index_search_params );
    
    /* build candidate mappings from matching subsets of the mapped_locations
     * for each indexable_subtemplate returned by the index search */
    candidate_mappings* mappings = NULL;
    init_candidate_mappings( &mappings );

    build_candidate_mappings_from_search_results(
            mappings, search_results, search_results_length,
            rst, genome, mapping_params );

    /* NOTE search_results contains references to memory allocated in the
     * indexable_subtemplates. search_results must be freed before ists */
    free_search_results( search_results, search_results_length );

    /* Return the set of gapped mappings - if mappings were ungapped, they were
     * included in this set of candidate mappings as-is (so for an ungapped
     * assay, gapped_mappings and mappings are identical). */
    build_gapped_candidate_mappings_for_read_subtemplate(
            mappings,
            final_mappings,
            rst,
            ists,
            genome,
            mapping_params->fwd_penalty_arrays[rst_pos],
            mapping_params->rev_penalty_arrays[rst_pos],
            mapping_params
        );

    free_candidate_mappings( mappings );
    free_indexable_subtemplates( ists );

    /* update the thread local copy of error data (need the error data
     * and the subtemplate to do this) */
    update_error_data_record_from_candidate_mappings(
            genome,
            final_mappings,
            rst,
            scratch_error_data_record
        );
    
    return;
}

void
update_read_type_pos(
        candidate_mappings* mappings,
        int rst_index )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        mappings->mappings[i].rd_type.pos = rst_index;
    }
}

void
find_candidate_mappings_for_read(
        struct read* r,
        candidate_mappings* read_mappings,
        struct genome_data* genome,

        struct error_model_t* error_model,
        struct error_data_record_t* scratch_error_data_record,
        enum bool only_collect_error_data,

        struct mapping_params* mapping_params
    )
{
    assert( NULL != error_model );
    
    int rst_index;
    for( rst_index=0; rst_index < r->num_subtemplates; rst_index++ )
    {
        // reference to current read subtemplate
        struct read_subtemplate* rst = r->subtemplates + rst_index;

        /* initialize the candidate mappings container for this read
         * subtemplate */
        candidate_mappings* rst_mappings = NULL;
        init_candidate_mappings( &rst_mappings );

        find_candidate_mappings_for_read_subtemplate(
                rst,
                rst_mappings,
                genome,
                
                scratch_error_data_record,
                only_collect_error_data,

                mapping_params
            );

        /* Update pos in READ_TYPE with the index of the underlying read
         * subtemplate for these candidate mappings */
        update_read_type_pos( rst_mappings, rst_index );

        /* append the candidate mappings from this read subtemplate to the set
         * of candidate mappings for this read */
        append_candidate_mappings( read_mappings, rst_mappings );

        free_candidate_mappings( rst_mappings );
    }
}

void
add_candidate_mappings_for_untemplated_gs(
        candidate_mappings* assay_corrected_mappings,
        candidate_mapping* cm,
        struct read* r
    )
{
    /** Add additional mappings if there is a G (that could be untemplated) in
     * the read sequence **/
    struct read_subtemplate* rst = r->subtemplates + cm->rd_type.pos;

    int i;
    for( i = 0; i < MAX_NUM_UNTEMPLATED_GS; i++ )
    {
        /* Get the BP we're looking at (after any soft clipped bases) */
        char bp = rst->char_seq[ cm->trimmed_length + i ];

        if( bp == 'G' || bp == 'g' )
        {
            candidate_mapping cm_copy = *cm;

            /* Soft clip the possibly untemplated G */
            cm_copy.mapped_length -= i+1;
            cm_copy.trimmed_length += i+1;

            if( cm_copy.rd_strnd == FWD )
            {
                cm_copy.start_bp += i+1;
                cm_copy.cigar[0].len -= i+1;
            } else {
                cm_copy.cigar[cm_copy.cigar_len-1].len -= i+1;
            }

            cm_copy.num_untemplated_gs += i+1;

            add_candidate_mapping( assay_corrected_mappings, &cm_copy );
        } else {
            break;
        }
    }

    return;
}

void
make_cage_specific_corrections(
        candidate_mappings* mappings,
        candidate_mappings* assay_corrected_mappings,
        struct read* r
    )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        candidate_mapping* cm = mappings->mappings + i;

        /* Add initial candidate mapping with no untemplated G's */
        candidate_mapping corrected_cm = *cm;

        /* Make sure we soft clipped enough bases when searching the index to be
         * able to make this adjustment */
        assert( corrected_cm.trimmed_length >= MAX_NUM_UNTEMPLATED_GS );

        /* Check that the candidate mapping isn't so close to a contig boundary
         * that this adjustment would cause it to cross over */
        if( corrected_cm.start_bp - MAX_NUM_UNTEMPLATED_GS < 0 ) {
            /* FIXME what's the right thing to do here? (skipping it for now) */
            /* Run the same algorithm with as many bp's as we can spare? */
            continue;
        }

        corrected_cm.mapped_length += MAX_NUM_UNTEMPLATED_GS;
        corrected_cm.trimmed_length -= MAX_NUM_UNTEMPLATED_GS;

        /* Adjust the entries in the cigar string to account for soft
         * clipping adjacent to the matches. */
        if( corrected_cm.rd_strnd == FWD )
        {
            /* Also adjust the start position of reads that map to the forward
             * strand */
            corrected_cm.start_bp -= MAX_NUM_UNTEMPLATED_GS;

            assert( corrected_cm.cigar[0].op == 'M' );
            corrected_cm.cigar[0].len += MAX_NUM_UNTEMPLATED_GS;

        } else { // corrected_cm.rd_strnd == BKWD
            assert( corrected_cm.cigar[corrected_cm.cigar_len-1].op == 'M' );
            corrected_cm.cigar[corrected_cm.cigar_len-1].len += MAX_NUM_UNTEMPLATED_GS;
        }

        add_candidate_mapping( assay_corrected_mappings, &corrected_cm );

        /* This candidate mapping is now the basis for adding additional
         * candidate mappings if it possible that they could have untemplated
         * G's */
        add_candidate_mappings_for_untemplated_gs( assay_corrected_mappings,
                &corrected_cm, r );
    }

    return;
}

void
make_assay_specific_corrections(
        candidate_mappings* mappings,
        candidate_mappings* assay_corrected_mappings,
        struct read* r
    )
{
    /* Only handle CAGE (untemplated G's) for now */
    if( _assay_type == CAGE )
    {
        make_cage_specific_corrections( mappings, assay_corrected_mappings, r );
    } else {
        /* Copy the original set of mappings to the corrected set (without
         * alteration, since this assay type does not require corrections) */
        copy_candidate_mappings( assay_corrected_mappings, mappings );
    }
        
    return;
}

mapped_read_t*
build_mapped_read_from_candidate_mappings(
        candidate_mappings* mappings,
        struct genome_data* genome,
        struct read* r,
        struct error_model_t* error_model,
        struct fragment_length_dist_t* fl_dist,

        struct mapping_params* mapping_params
    )
{            
    int joined_mappings_len = 0;
    candidate_mapping** joined_mappings = NULL;
    float* joined_mapping_penalties = NULL;
    mapped_read_t* rd = NULL;

    /* Build additional candidate mappings to make corrections for known
     * problems in the assay (e.g. untemplated G's in CAGE).
     *
     * For now, we assume these corrections all take places on the level of
     * a read subtemplate / candidate mapping.
     * */
    candidate_mappings* assay_corrected_mappings = NULL;
    init_candidate_mappings( &assay_corrected_mappings );
    make_assay_specific_corrections( mappings, assay_corrected_mappings, r );
    /* the original set of candidate mappings is freed in the calling fn */
    
    join_candidate_mappings( assay_corrected_mappings,
                             &joined_mappings,
                             &joined_mapping_penalties,
                             &joined_mappings_len );
            
    filter_joined_candidate_mappings( &joined_mappings,
                                      &joined_mapping_penalties,
                                      &joined_mappings_len,
                                              
                                      genome,
                                      r,
                                      error_model,
                                      fl_dist,
                                              
                                      mapping_params );

    rd = build_mapped_read_from_joined_candidate_mappings( r->read_id,
            joined_mappings, joined_mappings_len, joined_mapping_penalties );
    
    free_candidate_mappings( assay_corrected_mappings );
    
    free( joined_mappings );
    free( joined_mapping_penalties );
    
    return rd;
}

static inline void
increment_counter_with_lock( unsigned int *counter, pthread_spinlock_t* lock )
{
    pthread_spin_lock( lock );
    *counter += 1;
    pthread_spin_unlock( lock );
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
    
    struct genome_data* genome = td->genome;
    
    unsigned int* mapped_cnt = td->mapped_cnt;
    pthread_spinlock_t* mapped_cnt_lock = td->mapped_cnt_lock;

    struct rawread_db_t* rdb = td->rdb;
    
    struct mapped_reads_db* mpd_rds_db = td->mpd_rds_db;
    
    /* Store observed error data from the current thread's 
       execution in scratch */
    struct error_model_t* error_model = td->error_model;
    
    /* store statistics about mapping quality here  */
    struct error_data_t* error_data = td->error_data;

    /* if we only want error data, then there is not reason to find antyhing 
       except unique reads. */
    enum bool only_collect_error_data = td->only_collect_error_data;

    struct mapping_metaparams* metaparams = td->metaparams;
    
    /* END parameter 'recreation' */

    assert( genome->index != NULL );
    
    /* create a thread local copy of the error data to avoid excess locking */
    struct error_data_record_t* scratch_error_data_record;
    init_error_data_record( &scratch_error_data_record, 
                            error_data->max_read_length, 
                            error_data->max_qual_score );
    
    /* The current read of interest */
    struct read* r;
    int curr_read_index = 0;
    /* 
     * While there are still mappable reads in the read DB. All locking is done
     * in the get next read functions. 
     */
    while( EOF != get_next_read_from_rawread_db(
               rdb, &r, td->max_readkey )  
         ) 
    {
        /* We dont memory lock mapped_cnt because it's read only and we dont 
           really care if it's wrong 
         */
        if( r->read_id > 0 && 0 == r->read_id%MAPPING_STATUS_GRANULARITY )
        {
            fprintf(stderr, "DEBUG       :  Mapped %u reads, %i successfully\n", 
                    r->read_id, *mapped_cnt);
        }

        struct mapping_params* mapping_params = NULL;
        init_mapping_params_for_read(&mapping_params, r, metaparams, error_model);
        
        // Make sure this read has "enough" HQ bps before trying to map it
        if( filter_read( r, error_model ) )
        {
            add_unmappable_read_to_mapped_reads_db( r, mpd_rds_db );
            free_read( r );
            free_mapping_params( mapping_params );
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
                only_collect_error_data,
                mapping_params
            );
        
        /* unless we're only collecting error data, build candidate mappings */
        if( !only_collect_error_data )
        {
            mapped_read_t* mapped_read = 
                build_mapped_read_from_candidate_mappings(
                    mappings,
                    genome,
                    r,
                    error_model,
                    mpd_rds_db->fl_dist,
                    mapping_params
                );
            
            if( mapped_read != NULL )
            {
                /* mapped count is the number of reads that successfully mapped
                 * (not the number of mappings) */
                increment_counter_with_lock( mapped_cnt, mapped_cnt_lock );

                /* the read has at least one mapping - add it to the mapped
                 * reads database */
                add_read_to_mapped_reads_db( mpd_rds_db, mapped_read );
                free_mapped_read( mapped_read );
            } else {
                /* the read was declared mappable, but did not map */
                add_nonmapping_read_to_mapped_reads_db( r, mpd_rds_db );
            }
        } else {
            /* Update the number of bootstrapped reads (also uses mapped_cnt) */
            if( mappings->length > 0 )
                increment_counter_with_lock( mapped_cnt, mapped_cnt_lock );
        }

        curr_read_index += 1;

        /* cleanup memory */
        free_mapping_params( mapping_params );
        free_candidate_mappings( mappings );
        free_read( r );
    }

    /********* update the global error data *********/

    // add error_data to scratch_error_data
    merge_in_error_data_record( error_data, -1, scratch_error_data_record );
    
    // free local copy of error data
    free_error_data_record( scratch_error_data_record );

    return NULL;
}

void
spawn_find_candidate_mappings_threads( struct single_map_thread_data* td_template )
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
    size_t num_reads = 0;
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
        num_reads += (size_t) status;
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

                       struct mapping_metaparams* metaparams,
                       struct error_model_t* error_model,
                       struct error_data_t* error_data
    )
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
    
    td_template->metaparams = metaparams;
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
        struct mapping_metaparams* mapping_metaparams,
        struct error_model_t* error_model
    ) 
{
    assert( error_model != NULL );

    clock_t start = clock();
    
    struct error_model_t* bootstrap_error_model;
    init_error_model( &bootstrap_error_model, MISMATCH );

    struct error_data_t* error_data;
    init_error_data( &error_data );
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, 
                      mapping_metaparams, bootstrap_error_model, error_data   );
    
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

    spawn_find_candidate_mappings_threads( &td_template );
        
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
        
        struct mapping_metaparams* mapping_metaparams,
        struct error_model_t* error_model
    )
{
    clock_t start = clock();

    struct error_data_t* error_data;
    init_error_data( &error_data );
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, mpd_rds_db, 
                      mapping_metaparams, error_model, error_data );

    /* initialize the threads */
    while( false == rawread_db_is_empty( rdb ) )
    {
        // Add a new row to store error data in
        add_new_error_data_record( 
            td_template.error_data, 
            td_template.max_readkey, 
            td_template.max_readkey+READS_STAT_UPDATE_STEP_SIZE-1  );

        // update the maximum allowable readkey
        // update this dynamically
        td_template.max_readkey += READS_STAT_UPDATE_STEP_SIZE;
        
        spawn_find_candidate_mappings_threads( &td_template );
        
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
