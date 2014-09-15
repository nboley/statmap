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
                ref_gap = (match->locations[i].loc - match->subseq_offsets[i])
                        - (match->locations[i-1].loc - match->subseq_offsets[i-1]);
            } else {
                ref_gap = (match->locations[i-1].loc + match->subseq_offsets[i-1])
                        - (match->locations[i].loc + match->subseq_offsets[i]);
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
    for( i = 0; i < match->len; i++ ) {
        /* the product of the marginal (log) probabilites */
        cum_penalty += match->locations[i].penalty;
    }

    return cum_penalty;
}

candidate_mapping*
build_candidate_mapping_from_match(
        struct ml_match* match,
        struct read_subtemplate* rst,
        struct genome_data* genome )
{
    /* Make sure this is a valid, completed match */
    assert( match->len > 0 && match->matched == match->len );

    candidate_mapping* cm;
    cm = calloc( sizeof(candidate_mapping), 1 );

    mapped_location* base = match->locations + 0;
    cm->chr = base->chr;
    cm->rd_strnd = base->strnd;

    /* the length of the sequence that was mapped */
    cm->trimmed_length = softclip_len;
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
        norm_start_loc_index = match->len - 1;
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

    if( read_location < 0 ) // the read location was invalid; skip this matched set
        return NULL;
    cm->start_bp = read_location;
    
    build_candidate_mapping_cigar_string_from_match( cm, match, rst );

    cm->penalty = calc_candidate_mapping_penalty( cm, rst, genome );

    return cm;
}

candidate_mapping* 
build_ungapped_candidate_mapping_from_mapped_location(
    mapped_location* ml, 
    struct read_subtemplate* rst,
    struct indexable_subtemplate* ist,
    struct genome_data* genome) 
{
    struct ml_match* match = NULL;
    init_ml_match( &match, 1 );
    add_location_to_ml_match( ml, match, ist->subseq_length, ist->subseq_offset);
    candidate_mapping* mapping = build_candidate_mapping_from_match(
        match, rst, genome );
    free_ml_match(match);
    return mapping;
};

candidate_mappings*
build_candidate_mappings_from_search_results(
        mapped_locations** search_results,
        struct read_subtemplate* rst,
        struct genome_data* genome,
        struct mapping_params* mapping_params )
{
    // silence compiler warning
    assert( NULL != mapping_params );
    
    candidate_mappings* mappings = NULL;
    init_candidate_mappings( &mappings );

    /* sort each mapped_locations in order to use optimized merge algorithm */
    sort_search_results( search_results );
    
    /* consider each base location */
    int index_probe_i, i, j;
    for( index_probe_i = 0; 
         search_results[index_probe_i] != NULL; 
         index_probe_i++ )
    {
        /* Always use the mapped locations from the first indexable subtemplate
         * as the basis for building matches across the whole read subtemplate.
         * We always build matches from 5' -> 3' */
        mapped_locations* locs = search_results[index_probe_i];
        
        for( i = 0; i < locs->length; i++ )
        {
            mapped_location* loc = locs->locations + i;

            /* If the loc is a pseudo location, build a set of all its
             * possible expansions to consider for matching. Otherwise, returns
             * the original location */
            mapped_locations* expanded_locs = expand_base_mapped_locations(
                loc, locs->probe, genome );

            /* match across each of the expanded locations */
            for( j = 0; j < expanded_locs->length; j++ )
            {
                struct ml_match* match = NULL;
                init_ml_match( &match, 1 );
                add_location_to_ml_match( expanded_locs->locations + j, 
                                          match, 
                                          expanded_locs->probe->subseq_length,
                                          expanded_locs->probe->subseq_offset);
                candidate_mapping* mapping = 
                    build_candidate_mapping_from_match(match, rst, genome );
                if( NULL != mapping )add_candidate_mapping(mappings, mapping);
                free_ml_match(match);
            }
        
            free_mapped_locations( expanded_locs );
        }
    }

    return mappings;
}


int
find_candidate_mappings_for_read_subtemplate(
        struct read_subtemplate* rst,
        candidate_mappings** mappings,
        struct genome_data* genome,
        
        struct mapping_params* mapping_params
    )
{
    mapped_locations **search_results;
    *mappings = NULL;
    
    int rv = search_index_for_read_subtemplate( 
        rst, mapping_params, &search_results, genome, false);
    if( rv != 0 ) return rv;
    
    *mappings = build_candidate_mappings_from_search_results( 
        search_results, rst, genome, mapping_params);

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

/* returns 0 on success */
int
find_candidate_mappings_for_read(
        struct read* r,
        candidate_mappings* read_mappings,
        struct genome_data* genome,

        struct error_model_t* error_model,

        struct mapping_params* mapping_params
    )
{
    assert( NULL != error_model );
    int rv = -1;

    int i;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        // reference to current read subtemplate
        struct read_subtemplate* rst = r->subtemplates + i;

        /* initialize the candidate mappings container for this read
         * subtemplate */
        candidate_mappings* rst_mappings = NULL;
        rv = find_candidate_mappings_for_read_subtemplate(
                rst,
                &rst_mappings,
                genome,
                mapping_params
            );
        if( rv != 0 ) {
            if( NULL != rst_mappings )
                free_candidate_mappings( rst_mappings );
            
            free_candidate_mappings( read_mappings );
            return rv;
        }
        
        /* Update pos in READ_TYPE with the index of the underlying read
         * subtemplate for these candidate mappings */
        update_pos_in_template( rst_mappings, i );

        /* append the candidate mappings from this read subtemplate to the set
         * of candidate mappings for this read */
        append_candidate_mappings( read_mappings, rst_mappings );

        free_candidate_mappings( rst_mappings );
    }
    
    if( read_mappings->length > MAX_NUM_CAND_MAPPINGS ) {
        free_candidate_mappings(read_mappings);
        return TOO_MANY_CANDIDATE_MAPPINGS;
    }

    return 0;
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
        make_cage_specific_corrections( mappings, assay_corrected_mappings, r,
                genome );
    } else {
        /* Copy the original set of mappings to the corrected set (without
         * alteration, since this assay type does not require corrections) */

        /* FIXME can just return pointer to original mappings, eliminating
         * unnecessary memcpy. However, freeing the correct memory then becomes
         * a bit tricky. Revisit. */
        copy_candidate_mappings( assay_corrected_mappings, mappings );
    }
        
    return;
}

mapped_read_t*
build_mapped_read_from_candidate_mappings(
        candidate_mappings* mappings,
        struct genome_data* genome,
        struct read* r,
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
    make_assay_specific_corrections( 
        mappings, assay_corrected_mappings, r, genome );
    /* the original set of candidate mappings is freed in the calling fn */
    
    join_candidate_mappings( assay_corrected_mappings,
                             r,
                             &joined_mappings,
                             &joined_mapping_penalties,
                             &joined_mappings_len );
    
    if( joined_mappings_len > MAX_NUM_CAND_MAPPINGS ) {
        statmap_log( LOG_DEBUG, 
                     "Skipping read %i: too many candidate mappings ( %i )",
                     r->read_id, joined_mappings_len  );
        rd = NULL;
        goto cleanup;
    }
    
    filter_joined_candidate_mappings( &joined_mappings,
                                      &joined_mapping_penalties,
                                      &joined_mappings_len,
                                              
                                      genome,
                                      r,
                                      fl_dist,
                                              
                                      mapping_params );

    rd = build_mapped_read_from_joined_candidate_mappings( r->read_id,
            joined_mappings, joined_mappings_len, joined_mapping_penalties );
    
cleanup:
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
    
    struct mapping_metaparams* metaparams = td->metaparams;
    float reads_min_match_penalty = td->reads_min_match_penalty;
    
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
        #ifdef PROFILE_CANDIDATE_MAPPING
        statmap_log( LOG_DEBUG, "begin read_id %i", r->read_id );

        /* Log CPU time used by the current thread in processing this candidate mapping */
        int err;
        struct timespec start, stop;
        double elapsed;

        assert( sysconf(_POSIX_THREAD_CPUTIME) );
        err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        #endif

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
            = init_mapping_params_for_read( r, metaparams, error_model, 
                reads_min_match_penalty );
        cache_penalty_arrays_in_read_subtemplates(r, mapping_params);
        
        // Make sure this read has "enough" HQ bps before trying to map it
        if( filter_read( r, mapping_params, genome ) )
        {
            add_unmappable_read_to_mapped_reads_db( r, mpd_rds_db );
            free_read( r );
            free_mapping_params( mapping_params );
            continue; // skip the unmappable read
        }

        /* Initialize container for candidate mappings for this read */
        candidate_mappings* mappings = NULL;
        init_candidate_mappings( &mappings );
        
        int rv = find_candidate_mappings_for_read(
                    r,
                    mappings,
                    genome,
                    error_model,
                    mapping_params
            );
        
        if( rv != 0 )
        {
            add_unmappable_read_to_mapped_reads_db( r, mpd_rds_db );
            free_read( r );
            free_mapping_params( mapping_params );
            continue; // skip the unmappable read            
        }
        
        mapped_read_t* mapped_read = 
            build_mapped_read_from_candidate_mappings(
                mappings,
                genome,
                r,
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
        
        #ifdef PROFILE_CANDIDATE_MAPPING
        err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stop);
        elapsed = (stop.tv_sec - start.tv_sec);
        elapsed += (stop.tv_nsec - start.tv_nsec) / 1000000000.0;

        statmap_log(LOG_DEBUG, "find_candidate_mappings_time %f", elapsed);
        statmap_log(LOG_DEBUG, "end read_id %i", r->read_id);
        #endif

        /* cleanup memory */
        free_mapping_params( mapping_params );
        free_candidate_mappings( mappings );
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
    float reads_min_match_penalty = td->reads_min_match_penalty;
    
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
            = init_mapping_params_for_read( r, metaparams, error_model, 
                reads_min_match_penalty );
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
            
            // cleanup
            //free_search_results(search_results);
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
                statmap_log( LOG_FATAL, "Return code from pthread_join() is %d", rc );
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
        
        if( mapping_metaparams->error_model_type == ESTIMATED )
        {
            /* Compute the min match penalty for this block of reads that will
             * map the desired percentage of reads given in metaparameters */
            float reads_min_match_penalty
                = compute_min_match_penalty_for_reads( 
                    rdb, error_model, 
                    MAX(0, td_template.max_readkey - rdb->readkey),
                    mapping_metaparams->error_model_params[0] );

            statmap_log( LOG_INFO,
                "Computed min_match_penalty %f for reads [%i, %i]",
                reads_min_match_penalty,
                td_template.max_readkey - READS_STAT_UPDATE_STEP_SIZE,
                td_template.max_readkey );

            /* Save in the mapping metaparameters */
            td_template.reads_min_match_penalty = reads_min_match_penalty;
        }

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
