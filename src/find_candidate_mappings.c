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
#include <time.h>
#include <unistd.h> // for sysconf

#include "statmap.h"
#include "log.h"
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
#include "fragment_length.h"

/*******************************************************************************
 *
 * Psuedo location expansion code
 *
 ******************************************************************************/

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

mapped_locations*
expand_base_mapped_locations(
        mapped_location* base_loc,
        struct indexable_subtemplate* index_probe,
        struct genome_data* genome
    )
{
    mapped_locations* expanded_locs = NULL;
    init_mapped_locations( &expanded_locs, index_probe );
    
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

    return expanded_locs;
}


/*******************************************************************************
 * END Psuedo location expansion code
 ******************************************************************************/

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
    for( i = 0; i < match->matched; i++ )
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
                ref_gap = (match->locations[i].loc-match->subseq_offsets[i])
                       - (match->locations[i-1].loc-match->subseq_offsets[i-1]);
            } else {
                ref_gap = (match->locations[i-1].loc+match->subseq_offsets[i-1])
                        - (match->locations[i].loc+match->subseq_offsets[i]);
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

        /* cigar[cigar_index].len was initialized to zero when candidate
         * mapping was initialized */
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
    for( i = 0; i < match->matched; i++ ) {
        /* the product of the marginal (log) probabilites */
        cum_penalty += match->locations[i].penalty;
    }

    return cum_penalty;
}

enum bool
build_candidate_mapping_from_match(
    candidate_mapping* cm,
    struct ml_match* match,
    struct read_subtemplate* rst,
    struct genome_data* genome )
{
    mapped_location* base = match->locations + 0;
    cm->chr = base->chr;
    cm->rd_strnd = base->strnd;

    /* the length of the sequence that was mapped */
    cm->trimmed_length = softclip_len;
    cm->pos_in_template = -1;
    cm->trimmed_length = 0;
    cm->num_untemplated_gs = 0;
    cm->mapped_length = rst->length - cm->trimmed_length;

    /* The first location in the match is the first location in the mapping,
     * relative to the strand. We want to normalize to the start in the 5'
     * genome. */
    int norm_start_loc_index;

    assert( cm->rd_strnd == FWD || cm->rd_strnd == BKWD );
    if( cm->rd_strnd == FWD ) {
        /* If this is a forward mapped read, then the start relative to the
         * forward strand is in the first index probe's mapped location, since
         * index probes are chosen from 5' to 3' across the read */
        norm_start_loc_index = 0;
    } else // strnd == BKWD
    {
        /* If this a reverse mapped read, then the start relative to the
         * forward strand is the *end* of the read in the 3' genome - so we
         * want to use the last mapped_location */
        norm_start_loc_index = match->matched - 1;
    }

    int read_location
        = modify_mapped_read_location_for_index_probe_offset(
                match->locations[norm_start_loc_index].loc,
                match->locations[norm_start_loc_index].chr,
                match->locations[norm_start_loc_index].strnd,
                match->subseq_offsets[norm_start_loc_index],
                match->subseq_lengths[norm_start_loc_index],
                rst->length,
                genome
            );
    
    // the read location was invalid; skip this matched set
    if( read_location < 0 ) 
    {
        return false;
    }
    cm->start_bp = read_location;
    
    build_candidate_mapping_cigar_string_from_match( cm, match, rst );

    cm->penalty = calc_candidate_mapping_penalty( cm, rst, genome );

    cm->pos_in_template = rst->pos_in_template.pos;
    
    return true;
}

enum bool
build_ungapped_candidate_mapping_from_mapped_location(
    mapped_location* ml, 
    candidate_mapping* mapping,
    struct read_subtemplate* rst,
    struct indexable_subtemplate* ist,
    struct genome_data* genome) 
{
    struct ml_match* match;
    init_ml_match( &match, 1 );
    add_location_to_ml_match(ml, match, ist->subseq_length, ist->subseq_offset);
    int rv = build_candidate_mapping_from_match(
        mapping, match, rst, genome );
    free_ml_match(match);
    return rv;
};


int
compare_index_probes(mapped_location* loc1, 
                     int loc1_probe_offset,
                     mapped_location* loc2,
                     int loc2_probe_offset)
{
    if( loc1->strnd != loc2->strnd )
        return loc1->strnd - loc2->strnd;
    
    if( loc1->chr != loc2->chr )
        return loc1->chr - loc2->chr;
    
    assert( loc1->strnd == loc2->strnd );
    if(loc1->strnd == FWD ) 
    {
        return loc1->loc - loc1_probe_offset - loc2->loc + loc2_probe_offset;
    } else {
        return loc1->loc + loc1_probe_offset - loc2->loc - loc2_probe_offset;
    }
}

int
build_candidate_mappings_from_search_results(
        mapped_locations** search_results,
        candidate_mappings* mappings,
        struct read_subtemplate* rst,
        struct genome_data* genome )
{
    init_candidate_mappings( mappings );

    /* sort each mapped_locations in order to use optimized merge algorithm */
    sort_search_results( search_results );
    
    int i;
    
    /* initiaize the current probe indices */
    int num_probes;
    for(num_probes = 0; search_results[num_probes] != NULL; num_probes++ );
    
    int* curr_loc_indices = calloc(sizeof(int), num_probes);
    int* probe_offsets = calloc(sizeof(int), num_probes);
    int* probe_lengths = calloc(sizeof(int), num_probes);
    for(i = 0; i < num_probes; i++ ) {
        probe_offsets[i] = search_results[i]->probe->subseq_offset;
        probe_lengths[i] = search_results[i]->probe->subseq_length;
        curr_loc_indices[i]=0; ;
    }
    
    struct ml_match* mlm;
    init_ml_match( &mlm, MAX_NUM_INDEX_PROBES);

    while( true )
    {
        /* reset the ml match - this is jsut to avoid unnecessary 
           memory allocations */
        reset_ml_match(mlm);

        /*** increment the probe index ***/
        /* find the probe location with the smallest start position */
        int curr_probe_index = -1;
        for( i = 0; i < num_probes; i++ ) {
            /*Move past pseudo chromosomes */
            for(; curr_loc_indices[i] < search_results[i]->length
                    && ( search_results[i]->locations[curr_loc_indices[i]].chr 
                         == PSEUDO_LOC_CHR_INDEX );
             curr_loc_indices[i]++ );
            
            /* skip probes that don't have any remaing mapped locations */
            if( search_results[i]->length == curr_loc_indices[i] )
                continue;
            
            /* if we dont yet have a current probe or the current probe is 
               larger than i, upadte the current probe index */
            if( curr_probe_index < 0
                || 0 < compare_index_probes(
                    search_results[curr_probe_index]->locations 
                        + curr_loc_indices[curr_probe_index], 
                    probe_offsets[curr_probe_index],
                    search_results[i]->locations 
                        + curr_loc_indices[i], 
                    probe_offsets[i]))
            { curr_probe_index = i; }
        }
        
        /* if no probes are left, then we are done */
        if( curr_probe_index == -1 ) break;
                                        
        /* initialize the match structure, and add the first match, and 
           increment the current match point */
        add_location_to_ml_match(
            search_results[curr_probe_index]->locations 
                + curr_loc_indices[curr_probe_index], 
            mlm, 
            probe_lengths[curr_probe_index], 
            probe_offsets[curr_probe_index]);
        
        /* find matching probes */
        for( i = 0; i < num_probes; i++ ) 
        {
            /* if this is the current index, then it's already been added 
               to the match list */
            if( i == curr_probe_index ) continue;
            if( search_results[i]->length == curr_loc_indices[i] )
                continue;
            
            /* this must be greater than or equal to the base, because 
               of the previous step in which we chose the smallest */
            assert(0 >= compare_index_probes(
                       search_results[curr_probe_index]->locations 
                           + curr_loc_indices[curr_probe_index], 
                       probe_offsets[curr_probe_index],
                       search_results[i]->locations 
                           + curr_loc_indices[i], 
                       probe_offsets[i] ) );

            /* if this is a match, add it to the match lsit and 
               then incremenet the pointer */
            if(0 == compare_index_probes(
                       search_results[curr_probe_index]->locations 
                           + curr_loc_indices[curr_probe_index], 
                       probe_offsets[curr_probe_index],
                       search_results[i]->locations 
                           + curr_loc_indices[i], 
                       probe_offsets[i] ) )
            {
                add_location_to_ml_match(
                    search_results[i]->locations 
                        + curr_loc_indices[i], 
                    mlm, 
                    probe_lengths[i], 
                    probe_offsets[i]);
                /*increment the probe ndex, moving past any pseudo chromosomes*/
                for( curr_loc_indices[i] = curr_loc_indices[i] + 1;
                     curr_loc_indices[i] < search_results[i]->length
                       && (search_results[i]->locations[curr_loc_indices[i]].chr
                           == PSEUDO_LOC_CHR_INDEX );
                     curr_loc_indices[i]++ ); 
            }
        }

        /* increment the index for the base match */
        curr_loc_indices[curr_probe_index] += 1;
        
        candidate_mapping mapping;
        if( true == build_candidate_mapping_from_match(
                &mapping, mlm, rst, genome )  ) 
        {
            add_candidate_mapping(mappings, &mapping);
        }        
    }
    
    free_ml_match(mlm);
    
    return 0;
}


int
find_candidate_mappings_for_read_subtemplate(
        struct read_subtemplate* rst,
        candidate_mappings* mappings,
        struct genome_data* genome,
        
        struct mapping_params* mapping_params
    )
{
    mapped_locations **search_results;
    
    int rv = search_index_for_read_subtemplate( 
        rst, mapping_params, &search_results, genome, false);
    if( rv != 0 ) return rv;
    
    rv = build_candidate_mappings_from_search_results( 
        search_results, mappings, rst, genome);
    if( rv != 0 ) return rv;
    
    free_search_results(search_results);
    /* Return the set of gapped mappings - if mappings were ungapped, they were
     * included in this set of candidate mappings as-is (so for an ungapped
     * assay, gapped_mappings and mappings are identical). */
    /*
    build_gapped_candidate_mappings_for_read_subtemplate(
            mappings,
            final_mappings,
            rst,
            ists,
            genome,
            mapping_params
        );

    free_candidate_mappings( mappings );
    */
        
    return 0;
}

void
update_pos_in_template(
        candidate_mappings* mappings,
        int rst_index )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        mappings->mappings[i].pos_in_template = rst_index;
    }
}

int
paired_cms_frag_len(candidate_mapping* pair_1_mapping, 
                    candidate_mapping* pair_2_mapping)
{
    /* make sure this fragment length is allowed */
    int frag_start = MIN( pair_1_mapping->start_bp, 
                          pair_1_mapping->start_bp );
    int frag_stop = MAX( 
        pair_1_mapping->start_bp + pair_1_mapping->mapped_length, 
        pair_2_mapping->start_bp + pair_2_mapping->mapped_length );
    return frag_stop - frag_start;
}
/* returns 0 on success */
int
find_candidate_mappings_for_read(
        struct read* r,
        mapped_read_t** mapped_read,
        
        struct genome_data* genome,

        struct mapping_params* mapping_params,
        struct fragment_length_dist_t* fl_dist
    )
{
    int rv = -1;
    /* make sure that the reads are either single ended or paired */
    assert( r->num_subtemplates >= 1 && r->num_subtemplates <= 2 );
    candidate_mappings* rst_mappings[2];
    rst_mappings[0] = alloca(sizeof(candidate_mappings));
    init_candidate_mappings(rst_mappings[0]);
    rst_mappings[1] = alloca(sizeof(candidate_mappings));
    init_candidate_mappings(rst_mappings[1]);
    
    /* first find the candidate mappings for each read subtemplate */
    int i;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        // reference to current read subtemplate
        struct read_subtemplate* rst = r->subtemplates + i;
        
        /* initialize the candidate mappings container for this read
         * subtemplate */
        rv = find_candidate_mappings_for_read_subtemplate(
                rst,
                rst_mappings[i],
                genome,
                mapping_params
            );
                        
        if(rv != 0){ return rv; };
    }
    
    /* next pair them, making sure to keep track of the best match penalty */
    int r1_i, r2_i, num_joined_cms;
    float max_penalty = 2*mapping_params->recheck_min_match_penalty ;
    num_joined_cms = 0;
    struct joined_candidate_mappings joined_cms[MAX_NUM_CAND_MAPPINGS+1];
    for( r1_i = 0; r1_i < rst_mappings[0]->length; r1_i++ )
    {
        candidate_mapping* pair_1_mapping = &(rst_mappings[0]->mappings[r1_i]);        
        /* if single end, then there is no need to find a pair, so continue */
        if( r->num_subtemplates == 1 )
        {
            joined_cms[num_joined_cms] = (struct joined_candidate_mappings){
                pair_1_mapping, NULL, -1, pair_1_mapping->penalty };

            max_penalty = MAX(
                joined_cms[num_joined_cms].log_penalty, max_penalty);
            num_joined_cms += 1;
            continue;
        }
        
        for( r2_i = 0; r2_i < rst_mappings[1]->length; r2_i++ )
        {
            candidate_mapping* pair_2_mapping = &(
                rst_mappings[1]->mappings[r2_i]);

            /* If the chrs mismatch, since these are sorted we know there is no 
             * need to continue. */
            if( pair_2_mapping->chr > pair_1_mapping->chr )
                break;

            if( pair_2_mapping->chr < pair_1_mapping->chr )
                continue;

            assert( pair_1_mapping->chr == pair_2_mapping->chr );
            
            int frag_len = paired_cms_frag_len(pair_1_mapping, pair_2_mapping);
            if( frag_len > r->prior.max_fragment_length )
                continue;

            float fl_penalty = get_fl_log_prb( fl_dist, frag_len );
            
            float penalty = \
                pair_1_mapping->penalty + pair_2_mapping->penalty + fl_penalty;
            /*if we already know that the penalty is too small, then skip this*/
            if (penalty + mapping_params->recheck_max_penalty_spread 
                    < max_penalty )
                continue;
            max_penalty = MAX(max_penalty, penalty);
            
            /* add the joined mapping */
            joined_cms[num_joined_cms] = (struct joined_candidate_mappings){
                pair_1_mapping, pair_2_mapping, frag_len, penalty };
            num_joined_cms += 1;
            
            /* Update the array indices */
            if( num_joined_cms > MAX_NUM_CAND_MAPPINGS )
                return TOO_MANY_CANDIDATE_MAPPINGS_ERROR;
        }
    }

    if( 0 == num_joined_cms ) 
        return NO_JOINED_CANDIDATE_MAPPINGS;
    
    /* filter out matches that havea  penalty that is too small */
    int new_first_index = 0;
    for(i=0; i<num_joined_cms; i++)
    {
        if( max_penalty - mapping_params->recheck_max_penalty_spread 
            < joined_cms[i].log_penalty )
        {
            joined_cms[new_first_index] = joined_cms[i];
            new_first_index++;
        }
    }
    num_joined_cms = new_first_index;

    if( 0 == num_joined_cms ) 
        return NO_UNFILTERED_CANDIDATE_MAPPINGS;
    
    /* finally, build a mpped read from these */
    *mapped_read = 
        build_mapped_read_from_joined_candidate_mappings( 
            r->read_id, joined_cms, num_joined_cms);
    
    return rv;
}

void
add_candidate_mappings_for_untemplated_gs(
        candidate_mappings* assay_corrected_mappings,
        candidate_mapping* cm,
        struct read* r,
        struct genome_data* genome
    )
{
    /** Add additional mappings if there is a G (that could be untemplated) in
     * the read sequence **/
    struct read_subtemplate* rst = r->subtemplates + cm->pos_in_template;

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

            /* Add this soft clipped variation on the base candidate mapping if
             * it is a valid mapping */
            if( NULL != find_seq_ptr( genome, cm_copy.chr, cm_copy.start_bp,
                        cm_copy.mapped_length ) )
            {
                add_candidate_mapping( assay_corrected_mappings, &cm_copy );
            }
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
        struct read* r,
        struct genome_data* genome
    )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        /* Start by adjusting the mapping to be an initial candidate mapping
         * with no untemplated G's */
        candidate_mapping corrected_cm = mappings->mappings[i];

        /* Make sure we soft clipped enough bases when searching the index to be
         * able to make this adjustment */
        assert( corrected_cm.trimmed_length >= MAX_NUM_UNTEMPLATED_GS );
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

        /* If a soft clipped read maps near the boundary of a contig, it is
         * possible that this adjusted candidate mapping will be beyond the
         * boundary and thus invalid. However, soft clipped variations on the
         * mapping may be valid - so we don't add this to the set of corrected
         * mappings, but do use it as a basis to build possible valid soft
         * clipped mappings. */
        if( NULL != find_seq_ptr( genome, corrected_cm.chr,
                    corrected_cm.start_bp, corrected_cm.mapped_length ) )
        {
            add_candidate_mapping( assay_corrected_mappings, &corrected_cm );
        }

        /* This candidate mapping is now the basis for adding additional
         * candidate mappings if it possible that they could have untemplated
         * G's */
        add_candidate_mappings_for_untemplated_gs( assay_corrected_mappings,
                &corrected_cm, r, genome );
    }

    return;
}

void
make_assay_specific_corrections(
        candidate_mappings* mappings,
        candidate_mappings* assay_corrected_mappings,
        struct read* r,
        struct genome_data* genome
    )
{
    /* Only handle CAGE (untemplated G's) for now */
    if( _assay_type == CAGE )
    {
        make_cage_specific_corrections( 
            mappings, assay_corrected_mappings, r, genome );
    }
    return;
}

#if 0
int
build_mapped_read_from_candidate_mappings(
        candidate_mappings* mappings,
        mapped_read_t** rd,
        struct genome_data* genome,
        struct read* r,
        struct fragment_length_dist_t* fl_dist,

        struct mapping_params* mapping_params
    )
{            
    int rv = 0;    
    if( mappings->length == 0 ) {
        *rd = NULL;
        return NO_CANDIDATE_MAPPINGS;
    }
    

    /* Build additional candidate mappings to make corrections for known
     * problems in the assay (e.g. untemplated G's in CAGE).
     *
     * For now, we assume these corrections all take places on the level of
     * a read subtemplate / candidate mapping.
     * 
     
     */
    // XXX Dead code path
    //bool mappings_were_updated = make_assay_specific_corrections( 
    //    mappings, assay_corrected_mappings, r, genome );
    
    /* the original set of candidate mappings is freed in the calling fn */

    int joined_mappings_len = 0;
    candidate_mapping* joined_mappings[3*(MAX_NUM_CAND_MAPPINGS+1)];
    float joined_mapping_penalties[MAX_NUM_CAND_MAPPINGS+1];    
    join_candidate_mappings( mappings,
                             r,
                             joined_mappings,
                             joined_mapping_penalties,
                             &joined_mappings_len );
    
    if( joined_mappings_len > MAX_NUM_CAND_MAPPINGS ) {
        *rd = NULL;
        return TOO_MANY_CANDIDATE_MAPPINGS_ERROR;
    }

    if( joined_mappings_len == 0 ) {
        *rd = NULL;
        return NO_JOINED_CANDIDATE_MAPPINGS;
    }
    
    rv = filter_joined_candidate_mappings( joined_mappings,
                                           joined_mapping_penalties,
                                           joined_mappings_len,
                                           
                                           genome,
                                           r,
                                           fl_dist,
                                           
                                           mapping_params );

    if( rv != 0 ) {
        *rd = NULL;
        return rv;
    }

    *rd = build_mapped_read_from_joined_candidate_mappings( 
        r->read_id,
        joined_mappings, joined_mappings_len, joined_mapping_penalties );
    
    return rv;
}
#endif

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
    
    struct mapping_metaparams* metaparams = td->metaparams;
    
    /* END parameter 'recreation' */

    assert( genome->index != NULL );
    
    /* The current read of interest */
    struct read* r;
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
            statmap_log( LOG_INFO, "Mapped %u reads, %i successfully",  
                         r->read_id, *mapped_cnt );
        }

        if( r->read_id > 0 && 0 == (r->read_id%1000) )
        {
            statmap_log( LOG_DEBUG, "Mapped %u reads, %i successfully",  
                         r->read_id, *mapped_cnt );
        }

        struct mapping_params* mapping_params
            = init_mapping_params_for_read( r, metaparams, error_model );
        cache_penalty_arrays_in_read_subtemplates(r, mapping_params);
        
        /* Initialize container for candidate mappings for this read */        
        mapped_read_t* mapped_read;
        int rv = find_candidate_mappings_for_read(
                    r,
                    &mapped_read,
                    genome,
                    mapping_params,
                    mpd_rds_db->fl_dist
            );
        
        if( rv != 0 )
        {
            add_unmappable_read_to_mapped_reads_db( r, rv, mpd_rds_db );
            free_read( r );
            free_mapping_params( mapping_params );
            continue; // skip the unmappable read            
        }
        
        assert( mapped_read != NULL );
        /* mapped count is the number of reads that successfully mapped
         * (not the number of mappings) */
        increment_counter_with_lock( mapped_cnt, mapped_cnt_lock );

        /* the read has at least one mapping - add it to the mapped
         * reads database */
        add_read_to_mapped_reads_db( mpd_rds_db, mapped_read );
        free_mapped_read( mapped_read );
        
        /* cleanup memory */
        free_mapping_params( mapping_params );
        free_read( r );
    }
    
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
            statmap_log( LOG_FATAL, "Return code from pthread_create() is %d", rc );
        }
    }
    
    /* Free attribute and wait for the other threads */    
    size_t num_reads = 0;
    for(t=0; t < num_threads; t++) {
        rc = pthread_join(thread[t], &status);
        pthread_attr_destroy(attrs+t);
        if (rc) {
            statmap_log( LOG_FATAL, "Return code from pthread_join() is %d", rc );
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

    td_template->max_readkey = -1;

    td_template->rdb = rdb;
    
    td_template->mpd_rds_db = mpd_rds_db;
    
    td_template->metaparams = metaparams;
    td_template->reads_min_match_penalty = 0;
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

void* 
collect_error_data( void* params ) 
{
    /* 
     * recreate the struct parameters for readability
     * this should be optimized out 
     */

    struct single_map_thread_data* td = params;
    
    struct genome_data* genome = td->genome;
    
    struct rawread_db_t* rdb = td->rdb;
    
    struct error_model_t* error_model = td->error_model;
    
    /* store statistics about mapping quality here  */
    struct error_data_t* error_data = td->error_data;
    
    struct mapping_metaparams* metaparams = td->metaparams;
    
    /* END parameter 'recreation' */

    assert( genome->index != NULL );
    
    /* create a thread local copy of the error data to avoid excess locking */
    struct error_data_t* scratch_error_data;
    init_error_data( &scratch_error_data );
    
    /* The current read of interest */
    struct read* r;
    /* 
     * While there are still mappable reads in the read DB. All locking is done
     * in the get next read functions. 
     */
    while( EOF != get_next_read_from_rawread_db(
               rdb, &r, td->max_readkey )  
         ) 
    {
        if( rdb->readkey%1000 == 0 ) 
        {
            statmap_log(LOG_DEBUG, "Bootstrapped %i reads", rdb->readkey);
        }
        
        /* We dont memory lock mapped_cnt because it's read only and we dont 
           really care if it's wrong 
         */
        struct mapping_params* mapping_params
            = init_mapping_params_for_read( r, metaparams, error_model );
        cache_penalty_arrays_in_read_subtemplates(r, mapping_params);
        
        int i;
        for( i = 0; i < r->num_subtemplates; i++ )
        {
            // reference to current read subtemplate
            struct read_subtemplate* rst = r->subtemplates + i;
            
            mapped_locations** search_results;
            int rv = search_index_for_read_subtemplate(
                rst, mapping_params, &search_results, genome, true );
            if( rv != 0 ) continue;
            
            update_error_data_from_index_search_results(
                rst,
                search_results,
                genome,
                scratch_error_data
            );
            
            free_search_results(search_results);
        }
        
        free_mapping_params( mapping_params );
        free_read( r );
    }

    /********* update the global error data *********/

    // add error_data to scratch_error_data
    merge_in_error_data( error_data, scratch_error_data );
    
    // free local copy of error data
    free_error_data( scratch_error_data );

    return NULL;
}

/* bootstrap an already initialized error model */
void
bootstrap_estimated_error_model( 
        struct genome_data* genome,
        struct rawread_db_t* rdb,
        struct mapping_metaparams* mapping_metaparams,
        struct error_model_t* error_model
    ) 
{
    assert( error_model != NULL );
    
    int rc = 0;
    void* status;
    
    pthread_t thread[num_threads];
    pthread_attr_t attrs[num_threads];
    
    struct single_map_thread_data tds[num_threads];

    clock_t start = clock();
    
    struct error_model_t* bootstrap_error_model;
    init_error_model( &bootstrap_error_model, MISMATCH );

    struct error_data_t* error_data;
    init_error_data( &error_data );
    
    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    init_td_template( &td_template, genome, rdb, NULL, 
                      mapping_metaparams, bootstrap_error_model, error_data   );
        
    // Determine how many reads we should look through for the bootstrap
    td_template.max_readkey = NUM_READS_TO_BOOTSTRAP;
    
    if( num_threads == 1 )
    {
        collect_error_data(&td_template);
    } else {        
        int t;
        for( t = 0; t < num_threads; t++ )
        {  
            memcpy( tds+t,  &td_template, sizeof(struct single_map_thread_data) );
            tds[t].thread_id = t;
        
            pthread_attr_init(attrs + t);
            pthread_attr_setdetachstate(attrs + t, PTHREAD_CREATE_JOINABLE);
            pthread_create( thread + t, 
                            attrs + t, 
                            collect_error_data, 
                            (void *)(tds + t) );
            if (rc) {
                statmap_log( LOG_FATAL, "Return code from pthread_create() is %d", rc );
            }
        }
    
        /* Free attribute and wait for the other threads */    
        for(t=0; t < num_threads; t++) {
            rc = pthread_join(thread[t], &status);
            pthread_attr_destroy(attrs+t);
            if (rc) {
                statmap_log( LOG_FATAL, 
                             "Return code from pthread_join() is %d", rc );
            }
        }
    }
    free_td_template( &td_template );
    
    clock_t stop = clock();
    int num_bootstrapped_reads = 0;
    int i;
    for(i=0; i<error_data->num_records; i++) 
        num_bootstrapped_reads += error_data->records[i]->num_unique_reads;
    statmap_log( LOG_NOTICE,
            "Bootstrapped (%i/%u) Unique Reads in %.2lf seconds ( %e/thread-hour )",
            num_bootstrapped_reads, rdb->readkey,
            ((float)(stop-start))/CLOCKS_PER_SEC,
            ((float)(num_bootstrapped_reads*CLOCKS_PER_SEC*3600)/(stop-start))
        );
    
    update_error_model_from_error_data( error_model, error_data );
    
    FILE* error_data_ofp = fopen( BOOTSTRAP_ERROR_STATS_LOG, "w" );
    log_error_data( error_data_ofp, error_data );
    fclose( error_data_ofp );

    free_error_data( error_data );
    
    free_error_model( bootstrap_error_model );
    
    rewind_rawread_db( rdb );

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
        // update the maximum allowable readkey
        // update this dynamically
        td_template.max_readkey += READS_STAT_UPDATE_STEP_SIZE;
        spawn_find_candidate_mappings_threads( &td_template );
        
        /* update the error model from the new error data */
        #ifdef INCREMENTLY_UPDATE_ERROR_MODEL
        update_error_model_from_error_data(error_model, td_template.error_data);
        #endif
    }
    
    /* Print out performance information */
    clock_t stop = clock();
    statmap_log( LOG_NOTICE,
            "Mapped (%i/%u) Partial Reads in %.2lf seconds ( %e/thread-hour )",
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
