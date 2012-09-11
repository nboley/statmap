#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include <sys/stat.h>
/* mmap() is defined in this header */
#include <sys/mman.h> 

#include "genome.h"
#include "fragment_length.h"
#include "mapped_read.h"
#include "candidate_mapping.h"
#include "find_candidate_mappings.h"
#include "pseudo_location.h"

// this is needed for the wiggle writing code
#include "iterative_mapping.h"
#include "trace.h"

/*************************************************************************
 *
 *  Mapped Reads
 * 
 *  Reads that have been joined, but unlike mapped reads proper have not
 *  had the 'extra' read information attached.
 *
 */

char*
skip_read_id_nodes_in_mapped_read( char* ptr )
{
    /* Given a pointer to the start of a mapped_read_t, this returns a pointer
     * to the start of the mapped_read_locations (after any read_id_nodes) */

    /* loop over read_id_nodes */
    while(true)
    {
        read_id_node* curr_node = (read_id_node *) ptr;
        /* Save this read_id_node's are_more flag */
        int are_more = curr_node->are_more;

        /* Move the pointer */
        ptr += sizeof( read_id_node );

        /* If this was the last read id node, break */
        if( !are_more ) {
            break;
        }
    }

    return ptr;
}

char*
skip_mapped_read_sublocations_in_mapped_read_location( char* ptr )
{
    /* Loop over each mapped_read_sublocation
     * (we assume there is at least one) */
    while(true)
    {
        mapped_read_sublocation* curr_subloc = (mapped_read_sublocation*) ptr;

        /* If both next_subread flags are false, then there are no more
         * mapped_read_sublocations for this mapped_read_location */
        enum bool more_sublocs = true;
        if( !( curr_subloc->next_subread_is_gapped ||
               curr_subloc->next_subread_is_ungapped ) )
        {
            more_sublocs = false;
        }

        /* Skip this sub location */
        ptr += sizeof( mapped_read_sublocation );

        if( !more_sublocs )
            break;
    }

    return ptr;
}

char*
skip_mapped_read_location_in_mapped_read_locations( char* ptr )
{
    ptr += sizeof( mapped_read_location_prologue );
    ptr = skip_mapped_read_sublocations_in_mapped_read_location( ptr );

    return ptr;
}

char*
skip_mapped_read_locations_in_mapped_read_t( char* ptr )
{
    /* Loop over each mapped_read_location
     * (we assume there is at least one) */
    while(true)
    {
        /* Save the value of are_more for this mapped_read_location */
        mapped_read_location_prologue* curr_loc_prologue =
            (mapped_read_location_prologue*) ptr;

        int more_locs = curr_loc_prologue->are_more;

        ptr = skip_mapped_read_location_in_mapped_read_locations( ptr );

        if( !more_locs )
            break;
    }

    return ptr;
}

size_t
get_size_of_mapped_read_location(
        mapped_read_location* loc )
{
    assert( loc != NULL );

    /* Get a pointer we can use to iterate bytewise over the mapped_read_location */
    char* ptr = (char*) loc;

    /* Skip the prologue */
    ptr += sizeof( mapped_read_location_prologue );
    /* Skip the sublocations */
    ptr = skip_mapped_read_sublocations_in_mapped_read_location( ptr );

    return (size_t) (ptr - (char*) loc);
}

size_t
get_size_of_mapped_read( mapped_read_t* rd )
{
    /*
     * ASSUMPTIONS
     *
     * This code assumes that each mapped_read_t has at least
     * 1) 1 read_id_node
     * 2) 1 mapped_read_location, which has at least
     *     a) 1 mapped_read_sublocation
     */

    char* ptr = (char*) rd;
    ptr = skip_read_id_nodes_in_mapped_read( ptr );
    ptr = skip_mapped_read_locations_in_mapped_read_t( ptr );

    return (size_t)(ptr - (char*) rd);
}

void
free_mapped_read( mapped_read_t* rd )
{
    if( rd == NULL )
        return;
    
    free( rd );
    return;
}

size_t 
write_mapped_read_to_file( mapped_read_t* read, FILE* of  )
{
    size_t num_written = 0;
    size_t num_allocated = get_size_of_mapped_read( read );

    num_written = fwrite( read, sizeof(char), num_allocated, of );
    if( num_written != num_allocated )
        return -num_written;

    return 0;
}

MPD_RD_ID_T
get_read_id_from_mapped_read( mapped_read_t* rd )
{
    /* for now, assume each mapped_read_t has only one read_id_node.
     * Just take the first read id node and return the stored read id */
    char* rd_ptr = (char*) rd;
    
    /* cast to read_id_node to dereference */
    read_id_node *node = (read_id_node*) rd_ptr;
    return node->read_id;
}

/*************************************************************************
 *
 *  Mapped Read Index
 *
 */

void
index_mapped_read( mapped_read_t* rd,
                   mapped_read_index* index )
{
    const int REALLOC_BLOCK_SIZE = 1000;

    /* initialize array of pointers to mapped_read_location's */
    int num_allocated_locations = REALLOC_BLOCK_SIZE;
    index->mappings = malloc( num_allocated_locations *
                              sizeof(mapped_read_location*) );

    int num_mapped_locations = 0;

    /* Locate the start of the mapped read locations in this mapped read. */
    char* loc_ptr = skip_read_id_nodes_in_mapped_read( (char*) rd );

    /* index and count the mapped_read_location(s) */
    while(true)
    {
        /* Note - this assumes there is at least one mapped_read_location for
         * this mapped_read_t */

        /* resize index */
        if( num_mapped_locations + 1 == num_allocated_locations )
        {
            num_allocated_locations += REALLOC_BLOCK_SIZE;
            index->mappings = realloc( index->mappings,
                                       num_allocated_locations*
                                       sizeof(mapped_read_location*));
        }

        index->mappings[num_mapped_locations] = (mapped_read_location*) loc_ptr;
        num_mapped_locations += 1;

        /* cast read pointer to mapped_read_location_prologue */
        mapped_read_location_prologue* prologue =
            (mapped_read_location_prologue *) loc_ptr;

        /* if this was the last mapped_read_location, we're done counting */
        if( !( prologue->are_more ) )
        {
            break;
        }

        /* move rd_ptr to the start of the next mapped_read_location */
        loc_ptr = skip_mapped_read_location_in_mapped_read_locations( loc_ptr );
    }

    /* reclaim any wasted memory */
    index->mappings = realloc( index->mappings,
                               num_mapped_locations*
                               sizeof(mapped_read_location*) );
    /* save number of mapped locations */
    index->num_mappings = num_mapped_locations;

    return;
}

void
init_mapped_read_index( mapped_read_index** index,
                        mapped_read_t* rd )
{
    *index = malloc( sizeof( mapped_read_index ));

    (*index)->rd = rd;
    (*index)->read_id = get_read_id_from_mapped_read( rd );
    (*index)->num_mappings = 0;
    (*index)->mappings = NULL;

    /* Index the locations in this mapped read */
    index_mapped_read( rd, *index );

    return;
}

void
free_mapped_read_index( mapped_read_index* index )
{
    if( index == NULL ) return;

    if( index->mappings != NULL ) {
        free( index->mappings );
    }

    free( index );

    return;
}

/* Longstanding TODO -- expand joined_mappings dynamically for performance */

void
join_candidate_mappings_for_single_end( candidate_mappings* mappings, 
                                        candidate_mapping*** joined_mappings, 
                                        float** penalties,
                                        int* joined_mappings_len )
{
    /* for single end reads, just copy the candidate mappings over, adding
     * NULL delimiters */
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        *joined_mappings_len += 1;

        /* add new joined mapping to the joined_mappings
         * multiply by 2 for the NULL pointers separating entries */
        *joined_mappings = realloc( *joined_mappings,
                sizeof(candidate_mapping*) * (*joined_mappings_len * 2) );

        /* second to last pointer is the pointer to the current candidate
         * mapping */
        (*joined_mappings)[(*joined_mappings_len*2) - 2] = mappings->mappings + i;
        /* last pointer is NULL pointer indicating the end of this set of
         * joined candidate mappings */
        (*joined_mappings)[(*joined_mappings_len*2) - 1] = NULL;

        /* add the penalty for this new joined_mapping */
        *penalties = realloc( *penalties,
                sizeof(float) * *joined_mappings_len );
        (*penalties)[*joined_mappings_len - 1] = mappings->mappings[i].penalty;
    }

    return;
}

void
join_candidate_mappings_for_paired_end( candidate_mappings* mappings, 
                                        candidate_mapping*** joined_mappings, 
                                        float** penalties,
                                        int* joined_mappings_len )
{
    int pair_1_start = 0;
    int pair_2_start = -1;

    /* find the start of the second pair candidate mappings */
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        candidate_mapping* mapping = mappings->mappings + i;

        if( mapping->rd_type.pos == 1 )
        {
            pair_2_start = i;
            break;
        }
    }
    /* TODO instead of asserting, just don't build a mapped read if there are
     * no paired end 2 reads */
    assert( pair_2_start > 0 );

    /* sanity check */
    assert( pair_2_start > pair_1_start );

    int j;
    for( i = pair_1_start; i < pair_2_start; i++ )
    {
        for( j = pair_2_start; j < mappings->length; j++ )
        {
            candidate_mapping* pair_1_mapping = mappings->mappings + i;
            candidate_mapping* pair_2_mapping = mappings->mappings + j;

            /* If the chrs mismatch, since these are sorted we know there is no 
             * need to continue. */
            if( pair_2_mapping->chr > pair_1_mapping->chr )
                break;

            if( pair_2_mapping->chr < pair_1_mapping->chr )
                continue;

            assert( pair_1_mapping->chr == pair_2_mapping->chr );

            *joined_mappings_len += 1;

            /* allocate memory for new set of joined candidate mappings
             * multiply by 3 for 2 candidate mappings and 1 NULL pointer
             * separating the entries */
            *joined_mappings = realloc( *joined_mappings,
                    sizeof(candidate_mapping*) * (*joined_mappings_len * 3));

            /* third to last pointer is 1st candidate mapping */
            (*joined_mappings)[(*joined_mappings_len*3) - 3] = pair_1_mapping;

            /* second to last pointer is 2nd candidate mapping */
            (*joined_mappings)[(*joined_mappings_len*3) - 2] = pair_2_mapping;

            /* last pointer is NULL pointer indicating the end of this set of
             * joined candidate mappings */
            (*joined_mappings)[(*joined_mappings_len*3) - 1] = NULL;

            /* add the penalty for this (joined) candidate mapping. Simply
             * add the penalties from both candidate mappings (since these are
             * log probabilities, adding is equivalent to the product of the
             * marginal probabilities) */
            *penalties = realloc( *penalties,
                    sizeof(float) * *joined_mappings_len );
            (*penalties)[*joined_mappings_len - 1]
                = pair_1_mapping->penalty
                + pair_2_mapping->penalty;
        }
    }
}

void
join_candidate_mappings( candidate_mappings* mappings, 
                         candidate_mapping*** joined_mappings, 
                         float** penalties,
                         int* joined_mappings_len )
{
    /* joined_mappings is a list of pointers to candidate_mappings. Each
     * "joined mapping" is a list of pointers to candidate_mappings separated
     * by NULL poiners */

    /* sort the candidate_mappings for efficiency */
    sort_candidate_mappings( mappings );

    /* TODO special case single end and paired end reads joining, for now */

    enum bool paired_end = false;
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        candidate_mapping* mapping = mappings->mappings + i;

        if( mapping->rd_type.pos == 1 )
        {
            paired_end = true;
            break;
        }
    }

    if( paired_end )
    {
        join_candidate_mappings_for_paired_end( mappings,
                                                joined_mappings,
                                                penalties,
                                                joined_mappings_len );
    } else {
        join_candidate_mappings_for_single_end( mappings,
                                                joined_mappings,
                                                penalties,
                                                joined_mappings_len );
    }

    return;
}

void
recheck_candidate_mapping(
        candidate_mapping* mapping,

        struct read_subtemplate* rst,
        struct genome_data* genome,

        struct error_model_t* error_model )
{
    /* build penalty arrays from the underlying read subtemplate */
    struct penalty_array_t fwd_pa, rev_pa;
    build_penalty_array( &fwd_pa, rst->length,
            error_model, rst->error_str );
    build_reverse_penalty_array( &rev_pa, &fwd_pa );

    float rechecked_penalty = 0;

    int seq_pos = 0;
    int genome_pos = mapping->start_bp;

    /* loop over the entries in this candidate mapping's cigar string */
    int i;
    for( i = 0; i < mapping->cigar_len; i++ )
    {
        struct CIGAR_ENTRY current_entry = mapping->cigar[i];

        if( current_entry.op == 'M' )
        {
            /* This indicates a portion of matched sequence. Compare it to the
             * genome to recompute the penalty across the entire sequence */
            char* genome_seq = find_seq_ptr(
                    genome, mapping->chr, genome_pos, current_entry.len );

            /* if the genome_seq pointer is null, the sequence isn't valid
               for some reason ( ie, it falls off the end of the chromosome). 
               In this case, mark the location as invalid and continue TODO how? */
            if( NULL == genome_seq )
            {
                // TODO signal invalid
                return;
            }

            /* Build the sequence to compare and pick a penalty array depending
             * on the strand of the mapping
             *
             * TODO is this necesarily correct? Is it possible/logical for
             * different indexable subtemplates to have different strands? If
             * that were the case, we would just take the strand of the base
             * location for the cm's strand, and this would end up being incorrect. 
             */

            // TODO do we need the +1 here, since we know the length?
            char* real_seq = calloc( current_entry.len+1, sizeof(char) );
            assert( real_seq != NULL );
            struct penalty_array_t* penalty_array;

            if( BKWD == mapping->rd_strnd )
            {
                rev_complement_read( genome_seq, real_seq, current_entry.len );
                penalty_array = &rev_pa;
            } else {
                memcpy( real_seq, genome_seq, sizeof(char)*current_entry.len );
                penalty_array = &fwd_pa;
            }

            rechecked_penalty += recheck_penalty(
                    real_seq,
                    rst->char_seq + seq_pos,
                    current_entry.len,
                    penalty_array
                );

            seq_pos += current_entry.len;

            free( real_seq );

        } else if ( current_entry.op == 'N' ) {
            /* Move the genome_pos pointer forward to skip the intron indicated
             * by this cigar string entry */
            genome_pos += current_entry.len;
        } else {
            assert( false );
        }
    }

    mapping->penalty = rechecked_penalty;

    /* cleanup memory */
    free_penalty_array( &fwd_pa );
    free_penalty_array( &rev_pa );

    return;
}

int
get_fl_for_joined_candidate_mappings(
        candidate_mapping** group_start )
{
    /* There must be at least one candidate mapping in the joined group */
    assert( (*group_start) != NULL );

    /* Loop over the candidate mappings in the group until we get to the last
     * one */
    candidate_mapping** current_mapping = group_start;
    while( *(current_mapping + 1) != NULL )
        current_mapping++;

    /* The fragment length for the joined candidate mappings is the difference
     * between the start of the first cm and the start of the last cm, plus the
     * length given by the cigar string of the last candidate mapping. The
     * difference between the starts of the first and last cm's also includes
     * any gaps represented by the cigar strings of intervening
     * candidate_mapping's */

    /* TODO this might not be valid for all assay types. Verify */
    candidate_mapping* first_read = *group_start;
    candidate_mapping* last_read = *current_mapping;

    int frag_len = 0;
    int i;
    if( first_read->start_bp < last_read->start_bp )
    {
        frag_len = last_read->start_bp - first_read->start_bp;

        /* Add the length encoded in the last candidate mapping's cigar string */
        for( i = 0; i < last_read->cigar_len; i++ )
        {
            frag_len += last_read->cigar[i].len;
        }
    } else {
        frag_len = first_read->start_bp - last_read->start_bp;

        /* Add the length encoded in the first candidate mapping's cigar string */
        for( i = 0; i < first_read->cigar_len; i++ )
        {
            frag_len += first_read->cigar[i].len;
        }
    }

    return frag_len;
}

void
recheck_joined_candidate_mappings(
        candidate_mapping** group_start,
        float* penalty,
        
        struct genome_data* genome,
        struct read* r,
        struct error_model_t* error_model,
        struct fragment_length_dist_t* fl_dist )
{
    float rechecked_group_penalty = 0;

    candidate_mapping** current_mapping = group_start;
    while( (*current_mapping) != NULL )
    {
        /* each candidate_mapping in a joined group corresponds to a read
         * subtemplate. Get a reference to this cm's read subtemplate */
        struct read_subtemplate* rst
            = r->subtemplates + (*current_mapping)->rd_type.pos;

        recheck_candidate_mapping( *current_mapping, rst, genome, error_model );

        /* Add the log probability (equivalent to product) */
        rechecked_group_penalty += (*current_mapping)->penalty;

        /* Move pointer to the next candidate mapping in the group */
        current_mapping++;
    }

    /* Get the probability of the given fragment length from the fragment
     * length distribution */
    int fl = get_fl_for_joined_candidate_mappings( group_start );
    float fl_penalty = get_fl_log_prb( fl_dist, fl );

    /* Save the rechecked group penalty in penalty */
    *penalty = rechecked_group_penalty + fl_penalty;
}

void
advance_pointer_to_start_of_next_joined_candidate_mappings(
        candidate_mapping*** current_mapping,
        int current_set_index,
        int num_sets )
{
    /* skip the rest of the candidate_mappings in the current group */
    while( **current_mapping != NULL )
        (*current_mapping)++;

    /* unless we know we're in the last set of joined candidate mappings,
     * move to the pointer to the start of the next set. This check avoids
     * accessing uninitialized memory after the last set. */
    if( current_set_index < num_sets - 1 )
    {
        /* skip any NULL markers until we get to the start of the next set of
         * candidate_mappings */
        while( **current_mapping == NULL )
            (*current_mapping)++;
    }
}

void
remove_candidate_mapping_group(
        candidate_mapping** current_mapping )
{
    /* loop over the candidate mappings in this group (delimited by a NULL
     * pointer), and setting their pointers to NULL to remove them from the list
     * of joined candidate_mappings */

    while( *current_mapping != NULL )
    {
        /* remove the reference to the invalid candidate_mapping */
        *current_mapping = NULL;
        /* advance the pointer to the next candidate mapping pointer */
        current_mapping++;
    }
}

void
filter_joined_candidate_mappings( candidate_mapping*** joined_mappings, 
                                  float** penalties,
                                  int* joined_mappings_len,
                                  
                                  struct genome_data* genome,
                                  struct read* r,
                                  struct error_model_t* error_model,
                                  struct fragment_length_dist_t* fl_dist,

                                  float min_match_penalty,
                                  float max_penalty_spread )
{
    float max_penalty = min_match_penalty;

    /* Pointer to the start of the current set of joined candidate_mappings */
    candidate_mapping** current_mapping = *joined_mappings;

    int i;
    /* recheck the locations, updating penalties, and find the max penalty */
    for( i = 0; i < *joined_mappings_len; i++ )
    {
        recheck_joined_candidate_mappings( current_mapping, (*penalties) + i,
                genome, r, error_model, fl_dist );

        /* we may need to update the max penalty */
        if( (*penalties)[i] > max_penalty ) {
            max_penalty = (*penalties)[i];
        }

        advance_pointer_to_start_of_next_joined_candidate_mappings(
                &current_mapping, i, *joined_mappings_len );
    }

    /* Reset the pointer to the start of the joined mappings */
    current_mapping = (*joined_mappings);

    /* To filter mappings, we simply replace their pointers with NULL pointers
     * and update the joined_mappings_len appropriately */
    int filtered_mappings_len = *joined_mappings_len;

    // int i declared earlier
    for( i = 0; i < *joined_mappings_len; i++ )
    {
        enum bool filter_current_group = false;

        /* 
         * We always need to do this because of the way the search queue
         * handles poor branches. If our branch prediction fails we could
         * add a low probability read, and then go back and find a very 
         * good read which would invalidate the previous. Then, the read
         * may not belong, but it isn't removed in the index searching method.
         */

        /* we check max_penalty_spread > -0.00001 to make sure that the
           likelihood ratio threshold is actually active */

        /* I set the safe bit to 1e-6, which is correct for most floats */
        if( max_penalty_spread > 1e-6 )
        {
            assert( max_penalty_spread >= 0.0 ); // ?
            if( (*penalties)[i] < ( max_penalty - max_penalty_spread ) )
            {
                filter_current_group = true;
            }
        }

        /* Make sure that the penalty isn't too low */
        if( (*penalties)[i] < min_match_penalty )
        {
            filter_current_group = true;
        }

        /* if either test failed, filter the current group */
        if( filter_current_group )
        {
            remove_candidate_mapping_group( current_mapping );
            filtered_mappings_len -= 1;
        }

        advance_pointer_to_start_of_next_joined_candidate_mappings(
                &current_mapping, i, *joined_mappings_len );
    }

    *joined_mappings_len = filtered_mappings_len;

    return;
}

size_t
calculate_mapped_read_space_from_joined_candidate_mappings(
        candidate_mapping** joined_mappings,
        int joined_mappings_len )
{
    size_t size = 0;

    /* loop over the joined_mappings, counting the amount of space needed for
     * mapped_read_location's, mapped_read_sublocation's, etc. */
    candidate_mapping** current_mapping = joined_mappings;

    int i;
    for( i = 0; i < joined_mappings_len; i++ )
    {
        /* Every set of joined candidate mappings corresponds to
         * a mapped_read_location. We add the size of the mapped read location
         * prologue here */

        size += sizeof( mapped_read_location_prologue );

        /* for each candidate mapping in the current set of candidate mappings,
         * we will add a sublocation. */
        while( *current_mapping != NULL )
        {
            size += sizeof( mapped_read_sublocation );
            current_mapping++;
        }

        if( i < joined_mappings_len - 1 )
        {
            /* skip any NULL markers until we get to the start of the next set of
             * candidate_mappings */
            while( *current_mapping == NULL )
                current_mapping++;
        }
    }

    return size;
}

void
init_new_mapped_read_from_single_read_id( mapped_read_t** rd,
                                          MPD_RD_ID_T read_id,
                                          size_t mapped_read_size )
{
    /* Allocate memory for the mapped read */
    //*rd = malloc( mapped_read_size );
    /* Using calloc avoids an error in valgrind - 
     * "Syscall param write(buf) points to uninitialised byte(s)" */
    *rd = calloc( mapped_read_size, 1 );

    /* Add the read_id_node. For now, we assume there is only one */
    read_id_node* node = (read_id_node*) *rd;
    node->read_id = read_id;
    node->are_more = 0;
}

void
populate_mapped_read_locations_from_joined_candidate_mappings(
        mapped_read_t** rd,
        candidate_mapping** joined_mappings,
        int joined_mappings_len )
{
    /* loop over the joined candidate mappings, adding mapped locations as we
     * go */

    /* skip the read_id_nodes at the start of the mapped read to get to the
     * start of the first mapped_read_location */
    char* rd_ptr = skip_read_id_nodes_in_mapped_read( (char*) *rd );

    candidate_mapping** current_mapping = joined_mappings;

    /* move the current_mapping pointer to the start of the first set of joined
     * candidate_mappings */
    while( *current_mapping == NULL )
        current_mapping++;

    int i;
    for( i = 0; i < joined_mappings_len; i++ )
    {
        /* Add the mapped location prologue */
        mapped_read_location_prologue* prologue = 
            (mapped_read_location_prologue*) rd_ptr;

        /* TODO for now, take the chr and strand from the first candidate
         * mapping arbitrarily */
        prologue->chr = (*current_mapping)->chr;

        if( (*current_mapping)->rd_strnd == FWD ) {
            prologue->strand = 0;
        } else if ( (*current_mapping)->rd_strnd == BKWD ) {
            prologue->strand = 1;
        } else {
            assert( (*current_mapping)->rd_strnd == FWD ||
                    (*current_mapping)->rd_strnd == BKWD );
        }

        if( i == joined_mappings_len - 1 )
        {
            /* If this is the last set of joined mappings, this is the last
             * mapped_read_location we will add. Set are_more to 0 */
            prologue->are_more = 0;
        } else {
            /* There are more mapped_read_location's, set are_more to 1 */
            prologue->are_more = 1;
        }

        rd_ptr += sizeof( mapped_read_location_prologue );

        /* Sum the penalties of the candidate mappings ( since the penalties,
         * are log probabilites, we can just add (multiply) them, and do the
         * conversion back to [0,1] probability space once at the end */
        float cum_penalty = 0;

        while( *current_mapping != NULL )
        {
            /* Add a sublocation for each candidate mapping in the set of
             * joined candidate mappings */
            mapped_read_sublocation* subloc = 
                (mapped_read_sublocation*) rd_ptr;

            subloc->start_pos = (*current_mapping)->start_bp;
            subloc->length = (*current_mapping)->rd_len;

            if( (*current_mapping)->rd_strnd == FWD )
            {
                subloc->rev_comp = 0;
            } else if( (*current_mapping)->rd_strnd == BKWD ) {
                subloc->rev_comp = 1;
            } else {
                assert( (*current_mapping)->rd_strnd == FWD ||
                        (*current_mapping)->rd_strnd == BKWD );
            }

            /* TODO unused for now */
            subloc->is_full_contig = 0;

            // look ahead to figure out if the next subread is gapped or not.
            candidate_mapping* next_mapping = *(current_mapping + 1);
            if( next_mapping == NULL )
            {
                /* there is no next subread */
                subloc->next_subread_is_gapped = 0;
                subloc->next_subread_is_ungapped = 0;
            } else {
                if( next_mapping->rd_type.follows_ref_gap ) {
                    subloc->next_subread_is_gapped = 1;
                    subloc->next_subread_is_ungapped = 0;
                } else {
                    subloc->next_subread_is_gapped = 0;
                    subloc->next_subread_is_ungapped = 1;
                }
            }

            /* Add the penalty to the total */
            cum_penalty += (*current_mapping)->penalty;

            /* Advance pointers */
            rd_ptr += sizeof( mapped_read_sublocation );
            current_mapping++;
        }

        /* Convert the sum of the log penalties to a probability in standard
         * [0,1] probability space */
        prologue->seq_error = pow( 10, cum_penalty );
        assert( prologue->seq_error >= 0 && prologue->seq_error <= 1 );

        if( i < joined_mappings_len - 1 )
        {
            /* skip any NULL markers until we get to the start of the next set of
             * candidate_mappings */
            while( *current_mapping == NULL )
                current_mapping++;
        }
    }
}

mapped_read_t*
build_mapped_read_from_joined_candidate_mappings(
        MPD_RD_ID_T read_id,
        candidate_mapping** joined_mappings,
        int joined_mappings_len )
{
    if( joined_mappings_len == 0 )
        return NULL;
    
    mapped_read_t* rd = NULL;
    
    /* TODO for now, assume there is only one read_id_node */
    size_t mapped_read_size = sizeof( read_id_node );
    mapped_read_size +=
        calculate_mapped_read_space_from_joined_candidate_mappings(
                joined_mappings, joined_mappings_len );
    
    init_new_mapped_read_from_single_read_id( &rd, read_id, mapped_read_size );
    
    populate_mapped_read_locations_from_joined_candidate_mappings(
            &rd, joined_mappings, joined_mappings_len );
    
    return rd;
}

/*****************************************************************************
 *
 * Mapped Reads DB
 *
 ***************************************************************************/

static void
init_mapped_reads_db( 
    struct mapped_reads_db** rdb, char* fname, const char* mode )
{
    *rdb = malloc(sizeof(struct mapped_reads_db));

    /* Generic mutex attributes - used for all mutexes */
    /* we use a mutex because these operations are IO bound, plus it eliminates
     * the spurious helgrind warnings that we got when using a spinlock */
    pthread_mutexattr_t mta;
    pthread_mutexattr_init(&mta);

    /***** Mapped reads *****/

    (*rdb)->mapped_fp = fopen( fname, mode );
    if( (*rdb)->mapped_fp == NULL )
    {
        perror("FATAL       :  Could not open mapped reads file");
        assert( false );
        exit(-1);
    }
    
    /* number of mapped reads in the DB */
    (*rdb)->num_mapped_reads = 0;
    
    (*rdb)->mapped_mutex = malloc( sizeof(pthread_mutex_t) );
    pthread_mutex_init( (*rdb)->mapped_mutex, &mta );

    /***** Unmappable reads *****/

    char fname_buffer[255];

    sprintf( fname_buffer, "%s.unmappable", fname );
    (*rdb)->unmappable_fp = fopen( fname_buffer, mode );
    if( (*rdb)->unmappable_fp == NULL )
    {
        perror("FATAL       :  Could not open unmappable reads file");
        assert(false);
        exit(-1);
    }

    (*rdb)->num_unmappable_reads = 0;

    (*rdb)->unmappable_mutex = malloc( sizeof( pthread_mutex_t ) );
    pthread_mutex_init( (*rdb)->unmappable_mutex, &mta );

    /***** Nonmapping reads *****/

    sprintf( fname_buffer, "%s.nonmapping", fname );
    (*rdb)->nonmapping_fp = fopen( fname_buffer, mode );
    if( (*rdb)->nonmapping_fp == NULL )
    {
        perror("FATAL       :  Could not open nonmapping reads file");
        assert(false);
        exit(-1);
    }

    (*rdb)->num_nonmapping_reads = 0;

    (*rdb)->nonmapping_mutex = malloc( sizeof( pthread_mutex_t ));
    pthread_mutex_init( (*rdb)->nonmapping_mutex, &mta );
    
    /* Initialize the mode to 0, this will be set by
       the mode specific init function */
    (*rdb)->mode = 0;

    /* mmapped data */
    (*rdb)->mmapped_data = NULL;    
    (*rdb)->mmapped_data_size = 0;

    /* index */
    (*rdb)->index = NULL;
 
    /* fl dist */
    (*rdb)->fl_dist = NULL;

    /* store the number of times that we have iterated through the read
       db in the iterative mapping code. This is very hacky, and probably 
       shouldnt be here... TODO - remove this */
    (*rdb)->num_succ_iterations = NULL;

    /* the current read location ( either the number of reads written in 'w' 
       mode, or the read we're on in 'r' mode */
    (*rdb)->current_read = 0;

    return;
}

void
open_mapped_reads_db_for_reading(
        struct mapped_reads_db** rdb,
        char* fname
    )
{
    init_mapped_reads_db( rdb, fname, "r+" );

    mmap_mapped_reads_db( *rdb );

    // read num_mapped_reads from start of file
    (*rdb)->num_mapped_reads = *((MPD_RD_ID_T*) (*rdb)->mmapped_data);

    index_mapped_reads_db( *rdb );

    (*rdb)->mode = 'r';
}

void
open_mapped_reads_db_for_writing(
        struct mapped_reads_db** rdb,
        char* fname
    )
{
    init_mapped_reads_db( rdb, fname, "w+" );

    /* write placeholder for the number of mapped reads in the mapped reads db,
     * this will be updated when we close the mapped read db. */
    MPD_RD_ID_T placeholder = 0;
    fwrite( &placeholder, sizeof(MPD_RD_ID_T), 1, (*rdb)->mapped_fp );
    
    (*rdb)->mode = 'w';
}

void
close_reading_specific_portions_of_mapped_reads_db( struct mapped_reads_db* rdb)
{
    munmap_mapped_reads_db( rdb );

    if( rdb->index != NULL ) {
        free( rdb->index );
        rdb->index = NULL;
    }

    return;
}

void
close_writing_specific_portions_of_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* update the number of reads that we have written since the 
       file was opened. */
    fseek( rdb->mapped_fp, 0, SEEK_SET );
    fwrite( &(rdb->num_mapped_reads), sizeof(MPD_RD_ID_T), 1, rdb->mapped_fp );
    
    /* close the file pointers */
    fclose( rdb->mapped_fp );
    fclose( rdb->unmappable_fp );
    fclose( rdb->nonmapping_fp );

    return;
}

void
close_mapped_reads_db( struct mapped_reads_db** rdb )
{
    if( NULL == *rdb )
        return;

    if( (*rdb)->mode == 'r' )
    {
        close_reading_specific_portions_of_mapped_reads_db( *rdb );
    } else if( (*rdb)->mode == 'w' ) {
        close_writing_specific_portions_of_mapped_reads_db( *rdb );
    } else {
        fprintf( stderr,
                "FATAL       :  Unrecognized mode '%c' for open mapped reads db.\n",
                (*rdb)->mode );
        assert( false );
        exit( -1 );
    }
    
    /* clean up the mutexes */
    pthread_mutex_destroy( (*rdb)->mapped_mutex );
    free( (*rdb)->mapped_mutex );

    pthread_mutex_destroy( (*rdb)->unmappable_mutex);
    free( (*rdb)->unmappable_mutex );

    pthread_mutex_destroy( (*rdb)->nonmapping_mutex );
    free( (*rdb)->nonmapping_mutex );
    
    free( *rdb );
    *rdb = NULL;
    
    return;
}

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    mapped_read_t* rd)
{
    if ( rdb->mode == 'r' )
    {
        fprintf(stderr, "ERROR       :  Mapped Reads DB is read-only - cannot add read.\n");
        assert( false );
        exit( -1 );
    }

    int error;

    pthread_mutex_lock( rdb->mapped_mutex );
    error = write_mapped_read_to_file( rd, rdb->mapped_fp );
    rdb->num_mapped_reads += 1;
    pthread_mutex_unlock( rdb->mapped_mutex );

    if( error < 0 )
    {
        fprintf(stderr, "FATAL       :  Error writing to packed mapped reads db.\n");
        assert( false );
        exit( -1 );
    }
    
    return;
}

void
add_unmappable_read_to_mapped_reads_db(
        struct read* r,
        struct mapped_reads_db* db )
{
    /* lock the mutex for the corresponding file pointer */
    pthread_mutex_lock( db->unmappable_mutex );

    /* just write out the read id on one line */
    fprintf( db->unmappable_fp, "%d\n", r->read_id );
    db->num_unmappable_reads += 1;

    pthread_mutex_unlock( db->unmappable_mutex );
}

void
add_nonmapping_read_to_mapped_reads_db(
        struct read* r,
        struct mapped_reads_db* db )
{
    /* lock the mutex for the corresponding file pointer */
    pthread_mutex_lock( db->nonmapping_mutex );

    /* just write out the read id on one line */
    fprintf( db->nonmapping_fp, "%d\n", r->read_id );
    db->num_nonmapping_reads += 1;

    pthread_mutex_unlock( db->nonmapping_mutex );
}

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if( rdb->mode != 'r' )
    {
        fprintf(stderr, "FATAL       :  Can only rewind mapped reads db that is open for reading ( mode 'r' ).\n");
        assert( false );
        exit( -1 );
    }

    /* since the rdb is mmapped, we just need to reset the current read id */
    rdb->current_read = 0;
    
    return;
}

int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* rdb, 
    mapped_read_t** rd
)
{
    /* Make sure the db is open for reading */
    if( rdb->mode != 'r' )
    {
        fprintf(stderr, "FATAL       :  Cannot get read from mapped reads db unless it is open for reading.\n" );
        assert( false );
        exit( -1 );
    }

    /** Get the next read **/
    pthread_mutex_lock( rdb->mapped_mutex );
    /* if we have read every read */
    if( rdb->current_read == rdb->num_mapped_reads )
    {
        pthread_mutex_unlock( rdb->mapped_mutex );
        *rd = NULL;
        return EOF;
    }
    
    MPD_RD_ID_T current_read_id = rdb->current_read;
    rdb->current_read += 1;
    pthread_mutex_unlock( rdb->mapped_mutex );

    /* Set mapped_read_t to be a pointer into the mmapped mapped reads db */
    *rd = rdb->index[current_read_id].ptr;

    return 0;
}

void
reset_read_cond_probs( struct cond_prbs_db_t* cond_prbs_db,
                       mapped_read_t* rd,
                       struct mapped_reads_db* mpd_rds_db )
{
    struct fragment_length_dist_t* fl_dist = mpd_rds_db->fl_dist;
    
    /* build an index for this mapped_read */
    mapped_read_index* rd_index = NULL;
    init_mapped_read_index( &rd_index, rd );

    if( 0 == rd_index->num_mappings )
        return;
    
    float *prbs = calloc( rd_index->num_mappings, sizeof(float) );
    
    /* prevent divide by zero */
    double prb_sum = ML_PRB_MIN;
    MPD_RD_ID_T i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        mapped_read_location* loc = rd_index->mappings + i;

        float cond_prob = get_seq_error_from_mapped_read_location( loc );
        /* WARNING - this needs to be reconsidered. Simply removing the flag
         * makes this code incorrect. */
#if 0
        if( get_flag_from_mapped_read_location( loc )&IS_PAIRED )
            cond_prob *= get_fl_prb( fl_dist, get_fl_from_mapped_read_location( loc ) );
#endif
        prbs[i] = cond_prob;
        prb_sum += cond_prob;
    }
    assert( rd_index->num_mappings == 0 || prb_sum > ML_PRB_MIN );

    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        set_cond_prb( cond_prbs_db, rd_index->read_id, i, prbs[i]/prb_sum );
        
        assert( (prbs[i]/prb_sum < 1.001) && (prbs[i]/prb_sum) >= -0.001 );
    }
    
    free( prbs );
    free_mapped_read_index( rd_index );
    
    return;
}

void
reset_all_read_cond_probs( struct mapped_reads_db* rdb,
                           struct cond_prbs_db_t* cond_prbs_db )
{
    rewind_mapped_reads_db( rdb );
    mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) ) 
    {
        reset_read_cond_probs( cond_prbs_db, r, rdb );
    }
}


/* use this for wiggles */
void
update_traces_from_read_densities( 
    struct mapped_reads_db* rdb,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces
)
{    
    zero_traces( traces );

    mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) )     
    {
        mapped_read_index* rd_index;
        init_mapped_read_index( &rd_index, r );

        /* Update the trace from this mapping */
        MPD_RD_ID_T i;
        double cond_prob_sum = 0;
        for( i = 0; i < rd_index->num_mappings; i++ )
        {
            MRL_CHR_TYPE chr_index 
                = get_chr_from_mapped_read_location( rd_index->mappings + i );
            MRL_START_POS_TYPE start
                = get_start_from_mapped_read_location( rd_index->mappings + i );
            MRL_START_POS_TYPE stop
                = get_stop_from_mapped_read_location( rd_index->mappings + i );

            float cond_prob = get_cond_prb( cond_prbs_db, rd_index->read_id, i );
            cond_prob_sum += cond_prob;
            
            assert( cond_prob >= -0.0001 );
            assert( stop >= start );            
            assert( chr_index < traces->num_chrs );
            assert( traces->chr_lengths[chr_index] >= stop );

            unsigned int j = 0;
            for( j = start; j < stop; j++ )
            {
                /* update the trace */
                traces->traces[0][chr_index][j] 
                    += (1.0/(stop-start))*cond_prob;
            }
        }

        free_mapped_read_index( rd_index );
    }
    
    return;
}

void
mmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* get the file descriptor for the file we wish to mmap */
    int fdin = fileno( rdb->mapped_fp );
    
    /* make sure the entire file has been written to disk */
    fflush( rdb->mapped_fp );

    /* check that the file is not empty before trying to mmap */
    fseek(rdb->mapped_fp, 0L, SEEK_END);   // seek to end
    long fp_size = ftell(rdb->mapped_fp);  // tell() to get size
    if( fp_size == 0 )
    {
        fprintf( stderr,
                 "FATAL       :  Cannot mmap empty mapped reads db.\n" );
        exit( 1 );
    }
    rewind( rdb->mapped_fp ); // reset fp

    /* find the size of the opened file */
    struct stat buf;
    fstat(fdin, &buf);
    rdb->mmapped_data_size = buf.st_size;
    
    #ifdef MALLOC_READS_DB
    fseek( rdb->mapped_fp, 0, SEEK_SET );

    fprintf( stderr, 
             "NOTICE        : Allocating %zu bytes for the mapped reads db.", 
             buf.st_size );
    rdb->mmapped_data = malloc( buf.st_size );
    if( NULL == rdb->mmapped_data )
    {
        fprintf( stderr, 
                 "FATAL       : Failed to allocate %zu bytes for the mapped reads.\n",
                 (size_t) buf.st_size );
        exit( 1 );
    }
    
    fread( rdb->mmapped_data, buf.st_size, 1, rdb->mapped_fp );
           
    #else
    /* mmap the file */
    rdb->mmapped_data
        = mmap( NULL, rdb->mmapped_data_size,  
                PROT_READ|PROT_WRITE, 
		MAP_POPULATE|MAP_SHARED, fdin, (off_t) 0 );

    if( rdb->mmapped_data == (void*) -1 )
    {
        char* buffer;
        buffer = malloc( sizeof(char)*500 );
        sprintf(buffer, "Can not mmap the fdescriptor '%i'", fdin );
        perror( buffer );
        assert( false );
        exit( -1 );
    }
    #endif
    
    return;
}

void
munmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if( rdb->mmapped_data == NULL )
        return;

    #ifdef MALLOC_READS_DB
    free( rdb->mmapped_data );
    #else 
    int error = munmap( rdb->mmapped_data, rdb->mmapped_data_size );
    if( error != 0 )
    {
        perror( "Could not munmap mapped reads db" );
        assert( false );
        exit( -1 );
    }
    #endif

    rdb->mmapped_data = NULL;

    rdb->mmapped_data_size = 0;

    free( rdb->index );
    rdb->index = NULL;
    
    return;
}

int
cmp_mapped_reads_db_index_t(
        const struct mapped_reads_db_index_t* i1,
        const struct mapped_reads_db_index_t* i2
    )
{
    /* sort by read id */
    return i1->read_id - i2->read_id;
}

void
index_mapped_reads_db( struct mapped_reads_db* rdb )
{
    const int REALLOC_BLOCK_SIZE = 1000000;
    
    /* initialize dynamic array of mapped_reads_db_index_t */
    MPD_RD_ID_T num_allcd_reads = REALLOC_BLOCK_SIZE;
    rdb->index = malloc(sizeof(struct mapped_reads_db_index_t)*num_allcd_reads);

    /* Copy the reads data pointer (adding the offset from num_mapped_reads) */
    /* we use a char just to have a byte indexed memory block, meaning that we
       can use pointer ariuthmetic with sizeof */
    char* ptr = rdb->mmapped_data + sizeof(MPD_RD_ID_T);
    
    /* count mmapped reads to check they match the saved mapped reads count */
    MPD_RD_ID_T num_indexed_reads = 0;

    /* Loop through all of the reads */
    while( ((size_t)ptr - (size_t)rdb->mmapped_data)
            < rdb->mmapped_data_size )
    {
        /* start of a mapped_read_t */

        /* index the read_id_node's for this mapped_read_t */
        /* note - this assumes there is at least one read_id_node for each
         * mapped_read_t */
        while(true)
        {
            /* if necessary, resize the index */
            if( num_indexed_reads + 1 == num_allcd_reads )
            {
                num_allcd_reads += REALLOC_BLOCK_SIZE;
                rdb->index = realloc(rdb->index,
                        num_allcd_reads*sizeof(struct mapped_reads_db_index_t) );
            }
            assert( num_indexed_reads < num_allcd_reads );

            /* cast current location to read_id_node so it can be examined */
            read_id_node* node = (read_id_node*) ptr;

            /* add the read_id_node to the index */
            rdb->index[num_indexed_reads].read_id = node->read_id;
            rdb->index[num_indexed_reads].ptr = ptr;

            num_indexed_reads += 1;

            if( !(node->are_more) )
            {
                /* increment pointer to the start of the mapped read locations,
                 * and break */
                ptr += sizeof( read_id_node );
                break;
            }

            ptr += sizeof( read_id_node );
        }

        /* skip over the rest of the mapped_read_t in the mmapped memory */ 
        ptr = skip_mapped_read_locations_in_mapped_read_t( ptr );
    }

    /* reclaim any wasted memory */
    rdb->index = realloc( rdb->index,
            num_indexed_reads*sizeof(struct mapped_reads_db_index_t) );

    /* sort the index by read id - this restores synchronization with the
     * reads database */
    qsort( rdb->index,
           num_indexed_reads,
           sizeof(struct mapped_reads_db_index_t),
           (int(*)(const void*, const void*))cmp_mapped_reads_db_index_t
    );

    /* make sure that we have indexed every read */
    if( num_indexed_reads != rdb->num_mapped_reads )
    {
        fprintf( stderr, 
                 "FATAL           :  The number of indexed reads (%i) is not equal to the number of reads in the mapped read db ( %i). This may indicate that the mapped read db is corrupt.", 
                 num_indexed_reads, rdb->num_mapped_reads );

        assert(false);
        exit(-1);
    }
    
    return;
}

/* for debugging */
void
print_mapped_reads_db_index(
        struct mapped_reads_db* rdb
    )
{
    MPD_RD_ID_T i;
    for( i = 0; i < rdb->num_mapped_reads; i++ )
    {
        fprintf(stderr, "read_id: %u, ptr: %p\n",
                rdb->index[i].read_id, rdb->index[i].ptr );
    }
}

/*
 *  END Mapped DB Reads
 *
 **************************************************************************/


/*****************************************************************************
 *
 * Conditional Probs DB Code
 *
 ***************************************************************************/


void
init_cond_prbs_db_from_mpd_rdb( 
    struct cond_prbs_db_t** cond_prbs_db,
    struct mapped_reads_db* mpd_rdb
)
{
    *cond_prbs_db = malloc( sizeof( struct cond_prbs_db_t ) );
    
    /* reset the database position */
    rewind_mapped_reads_db( mpd_rdb );
    
    /* find the maximum readid */
    MPD_RD_ID_T max_rd_id = 0;
    mapped_read_t* mapped_rd;
    while( EOF != get_next_read_from_mapped_reads_db( mpd_rdb, &mapped_rd ) )
    {
        MPD_RD_ID_T read_id = get_read_id_from_mapped_read( mapped_rd );
        max_rd_id = MAX( max_rd_id, read_id );
    }
    (*cond_prbs_db)->max_rd_id = max_rd_id;
    
    /* allocate space for the prb start pointers */
    (*cond_prbs_db)->cond_read_prbs = calloc( max_rd_id+1, sizeof(ML_PRB_TYPE*) );
    
    /* allocate space for the prbs */
    rewind_mapped_reads_db( mpd_rdb );
    while( EOF != get_next_read_from_mapped_reads_db( mpd_rdb, &mapped_rd ) )
    {
        mapped_read_index* rd_index;
        init_mapped_read_index( &rd_index, mapped_rd );

        (*cond_prbs_db)->cond_read_prbs[rd_index->read_id] 
            = calloc( rd_index->num_mappings, sizeof( ML_PRB_TYPE  )  );

        free_mapped_read_index( rd_index );
    }

    return;
}

void
free_cond_prbs_db( struct cond_prbs_db_t* cond_prbs_db )
{
    MPD_RD_ID_T i;
    for( i = 0; i < cond_prbs_db->max_rd_id+1; i++ )
    {
        if( cond_prbs_db->cond_read_prbs[i] != NULL )
            free( cond_prbs_db->cond_read_prbs[i] );
    }
    
    free( cond_prbs_db->cond_read_prbs );
    free( cond_prbs_db );
}

/*
 *  END Conditional Prbs DB Code
 *
 **************************************************************************/

