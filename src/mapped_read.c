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

#include "log.h"

/*******************************************************************************
 *
 * Debugging mapped reads
 *
 */

void
print_mapped_read_location(
        mapped_read_location* loc
    )
{
    /* cast the mapped read location to mapped_read_location_prologue and
       get the prologue information */
    mapped_read_location_prologue* prologue
        = (mapped_read_location_prologue*) loc;

    fprintf( stderr,
        "chr:%i, strand:%i, are_more:%i, trimmed_length:%i, seq_error:%f\n",
        prologue->chr, prologue->strand, prologue->are_more,
        prologue->trimmed_length, prologue->seq_error );

    /* print out each of the sublocations */
    char* rd_ptr = loc;

    /* skip the prologue */
    rd_ptr += sizeof(mapped_read_location_prologue);

    /* read until there are no more sublocations */
    while(true)
    {
        mapped_read_sublocation* sub_loc
            = (mapped_read_sublocation*) rd_ptr;

        fprintf( stderr, 
"start:%i, length:%i, rev_comp:%i, is_full_contig:%i, next_subread_is_gapped:%i, next_subread_is_ungapped:%i\n",
            sub_loc->start_pos, sub_loc->length, sub_loc->rev_comp,
            sub_loc->is_full_contig, sub_loc->next_subread_is_gapped,
            sub_loc->next_subread_is_ungapped );

        if( last_sublocation_in_mapped_read_location(sub_loc) )
            break;

        rd_ptr += sizeof(mapped_read_sublocation);
    }
}

void
print_mapped_read(
        mapped_read_t* rd
    )
{
    /* Index the mapped read to make this easier (be careful, if there's a bug 
       in indexing mapped reads this will be messed up - then again, so will
       everything else) */
    mapped_read_index* rd_index = build_mapped_read_index(rd);

    fprintf(stderr, "==== MAPPED READ ID %i (%i mappings) ====\n",
        rd_index->read_id, rd_index->num_mappings);

    /* print out each of the mappings */
    int i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        print_mapped_read_location( rd_index->mappings[i] );
    }
}

/** MAPPED READ SUBLOCATIONS **/
mapped_read_sublocation*
get_start_of_sublocations_in_mapped_read_location(
        const mapped_read_location* loc )
{
    // skip the location prologue
    char* ptr = (char*) loc;
    ptr += sizeof(mapped_read_location_prologue);

    return (mapped_read_sublocation*) ptr;
}

int
get_start_for_sublocations_group(
        mapped_read_sublocation* group )
{
    /* simply return the start of the first sublocation in the group */
    return group->start_pos;
}

int
get_stop_for_sublocations_group(
        mapped_read_sublocation* group )
{
    /* loop through to the last sublocation, and return its
     * start + length */

    char* ptr = (char*) group;

    while( true )
    {
        if( last_sublocation_in_mapped_read_location(
                (mapped_read_sublocation*) ptr ) )
            break;

        ptr += sizeof( mapped_read_sublocation );
    }

    mapped_read_sublocation* final_sublocation
        = (mapped_read_sublocation*) ptr;

    return final_sublocation->start_pos + final_sublocation->length;
}

enum bool
last_sublocation_in_mapped_read_location( mapped_read_sublocation* sub )
{
    if( !sub->next_subread_is_ungapped &&
        !sub->next_subread_is_gapped )
    {
        /* then there is no "next" subread, and this is the last sublocation in
         * the mapped_read_location */
        return true;
    }

    return false;
}

enum bool
end_of_sublocations_group( mapped_read_sublocation* sub_loc )
{
    /* Sublocations are grouped by consecutive gapped sublocations.
     *
     * If a sublocation is either
     * 1) The last sublocation in a group of gapped sublocations
     * 2) The last sublocation in a mapped_read_location
     *
     * it is at the end of a group of sublocations.
     */
    if( sub_loc->next_subread_is_ungapped ||
            last_sublocation_in_mapped_read_location( sub_loc ))
        return true;

    return false;
}

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
        /* Assume there is at least one mapped_read_location in a mapped_read_t */

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

    /* reclaim any unused memory */
    index->mappings = realloc( index->mappings,
                               num_mapped_locations*
                               sizeof(mapped_read_location*) );

    /* save number of mapped locations */
    index->num_mappings = num_mapped_locations;

    return;
}

mapped_read_index*
build_mapped_read_index( mapped_read_t* rd )
{
    mapped_read_index* index = malloc( sizeof( mapped_read_index ) );
    assert( index != NULL );

    index->rd = rd;
    index->read_id = get_read_id_from_mapped_read( rd );
    index->num_mappings = 0;
    index->mappings = NULL;

    /* Index the locations in the mapped read */
    index_mapped_read( rd, index );

    return index;
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

/*******************************************************************************
*
*   Joining candidate mappings
*
*/

void
join_candidate_mappings_for_single_end( candidate_mappings* mappings, 
                                        candidate_mapping*** joined_mappings, 
                                        float** penalties,
                                        int* joined_mappings_len )
{
    /* Allocate memory */
    /* For the single end case, there are N potential mappings. Allocate N*2
     * pointers for the NULL delimiters */
    *joined_mappings = malloc( sizeof(candidate_mapping*) * mappings->length * 2 );
    *penalties = malloc( sizeof(float) * mappings->length );

    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        /* second to last pointer is the pointer to the current candidate
         * mapping */
        (*joined_mappings)[(i+1)*2-2] = mappings->mappings + i;
        /* last pointer is NULL pointer indicating the end of this set of
         * joined candidate mappings */
        (*joined_mappings)[(i+1)*2-1] = NULL;

        (*penalties)[i] = mappings->mappings[i].penalty;
    }

    *joined_mappings_len = mappings->length;

    return;
}

void
join_candidate_mappings_for_paired_end( candidate_mappings* mappings, 
                                        struct read* r,
                                        candidate_mapping*** joined_mappings, 
                                        float** penalties,
                                        int* joined_mappings_len )
{
    int pair_1_start = 0;
    int pair_2_start = -1;
    int num_pair_1 = 0;
    int num_pair_2 = 0;

    /* find the start of the second pair candidate mappings */
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        candidate_mapping* mapping = mappings->mappings + i;

        if( mapping->pos_in_template == 1 && pair_2_start == -1 )
        {
            pair_2_start = i;
            num_pair_1 = i;
        }
    }

    if( num_pair_1 <= 0 || pair_2_start == -1 )
    {
        /* Don't build a mapped read if each pair doesnt have a read */
        return;
    }
    /* sanity check */
    assert( pair_2_start > pair_1_start );
    num_pair_2 = (i-1) - pair_2_start;

    /* Build the initial memory allocation.
     * For each set of joined mappings, we need
     * to allocate 3 pointers - one for the first mapping, one for the second
     * mapping, and one for the NULL delimiter. */
    const int allcd_step_size = 100;
    int allcd_size = allcd_step_size;
    *joined_mappings = malloc(
            sizeof(candidate_mapping*) * (allcd_size*3));
    *penalties = malloc( sizeof(float) * allcd_size );

    int num_joined_mappings = 0;
    for( i = pair_1_start; i < pair_2_start; i++ )
    {
        int j;
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
            
            /* make sure this fragment length is allowed */
            int frag_start = MIN( pair_1_mapping->start_bp, 
                                  pair_1_mapping->start_bp );
            int frag_stop = MAX( 
                pair_1_mapping->start_bp + pair_1_mapping->mapped_length, 
                pair_2_mapping->start_bp + pair_2_mapping->mapped_length );
            int frag_len = frag_stop - frag_start;
            if( frag_len > r->prior.max_fragment_length )
                continue;
            
            /* third to last pointer is 1st candidate mapping */
            (*joined_mappings)[3*num_joined_mappings] = pair_1_mapping;

            /* second to last pointer is 2nd candidate mapping */
            (*joined_mappings)[3*num_joined_mappings + 1] = pair_2_mapping;

            /* last pointer is NULL pointer indicating the end of this set of
             * joined candidate mappings */
            (*joined_mappings)[3*num_joined_mappings + 2] = NULL;

            /* add the penalty for this (joined) candidate mapping. Simply
             * add the penalties from both candidate mappings (since these are
             * log probabilities, adding is equivalent to the product of the
             * marginal probabilities) */
            (*penalties)[num_joined_mappings]
                = pair_1_mapping->penalty
                + pair_2_mapping->penalty;

            /* After we've used this as an index, make it a count */
            num_joined_mappings += 1;
            
            /* Update the array indices */
            if( num_joined_mappings > MAX_NUM_CAND_MAPPINGS )
                goto cleanup;

            if( num_joined_mappings == allcd_size )
            {
                allcd_size += allcd_step_size;
                *joined_mappings = realloc( 
                    *joined_mappings, 
                    sizeof(candidate_mapping*)*(allcd_size*3)
                );
                *penalties = realloc( *penalties, sizeof(float)*allcd_size );
            }
        }
    }

cleanup:
    
    /* recalim unused memory */
    *joined_mappings = realloc( 
        *joined_mappings, 
        sizeof(candidate_mapping*)*(num_joined_mappings*3)
    );
    *penalties = realloc( *penalties, sizeof(float)*num_joined_mappings );


    *joined_mappings_len = num_joined_mappings;

    return;
}

void
join_candidate_mappings( candidate_mappings* mappings, 
                         struct read* r,
                         candidate_mapping*** joined_mappings, 
                         float** penalties,
                         int* joined_mappings_len
                       )
{
    /* sort the candidate_mappings so we can optimize */
    sort_candidate_mappings( mappings );
    
    /* decide which joining function to call */
    if( r->num_subtemplates == 1 )
    {
        join_candidate_mappings_for_single_end( mappings,
                                                joined_mappings,
                                                penalties,
                                                joined_mappings_len );
    } else if( r->num_subtemplates == 2 ) {
        join_candidate_mappings_for_paired_end( mappings,
                                                r,
                                                joined_mappings,
                                                penalties,
                                                joined_mappings_len );
    } else {
        statmap_log(LOG_FATAL, "Only 2 or less subtempaltes are currently supported in candidate read joining.");
    }

    return;
}

void
recheck_candidate_mapping(
        candidate_mapping* mapping,
        struct read_subtemplate* rst,
        struct genome_data* genome,
        struct error_model_t* error_model
    )
{
    float rechecked_penalty = 0;

    /* build penalty array from the underlying read subtemplate.
     *
     * Since we actually reverse complement the full *genome*
     * sequence in the recheck (to simplify the algorithm when rechecking
     * gapped reads), the read is being compared to the sequence that it mapped
     * to. Therefore, we only need to build the fwd penalty array. */
    struct penalty_array_t pa;
    build_penalty_array( &pa, rst, error_model );

    /* Get the entire region of the genome covered by this fragment */
    int full_fragment_length = get_length_from_cigar_string( mapping);
    char* genome_seq = find_seq_ptr( genome, mapping->chr, mapping->start_bp, 
            full_fragment_length );
    assert( genome_seq != NULL );

    /* Build the sequence to compare and a pick a penalty array
     * depending on strand */
    char* mapped_seq = calloc( full_fragment_length+1, sizeof(char) );
    assert( mapped_seq != NULL );

    if( mapping->rd_strnd == FWD )
    {
        memcpy( mapped_seq, genome_seq, full_fragment_length );
    } else {
        rev_complement_read( genome_seq, mapped_seq, full_fragment_length );
    }

    /* loop over the entries in this candidate mapping's cigar string */
    int ref_pos = 0;
    // skip any initial soft clipping (currently not represented in the cigar string)
    // TODO - use the cigar string?
    int read_pos = mapping->trimmed_length; 

    int i;
    for( i = 0; i < mapping->cigar_len; i++ )
    {
        int cigar_index;
        if( mapping->rd_strnd == FWD ) {
            cigar_index = i;
        } else {
            cigar_index = mapping->cigar_len - i - 1;
        }

        struct CIGAR_ENTRY entry = mapping->cigar[cigar_index];

        if( entry.op == 'M' )
        {
            rechecked_penalty += recheck_penalty(
                    mapped_seq + ref_pos,
                    rst->char_seq + read_pos,
                    pa.array + read_pos,
                    entry.len
                );

            /* update both the ref and read position */
            ref_pos += entry.len;
            read_pos += entry.len;
        } else if( entry.op == 'N' ) {
            /* just update the genome pointer to skip the intron */
            ref_pos += entry.len;
        } else if( entry.op == 'S' ) {
            /* just update the read pointer to skip the soft clipped bases */
            read_pos += entry.len;
        } else {
            statmap_log( LOG_FATAL, "Found unsupported CIGAR string op '%c'",  entry.op  );
            assert( false );
            exit(-1);
        }
    }

    /* The rechecked penalty is the total penalty of any match (M) segments,
     * plus the marginal probability for any untemplated G's (which are soft
     * clipped, but recorded in mapping->num_untemplated_gs) */
    mapping->penalty = rechecked_penalty + 
        (mapping->num_untemplated_gs*UNTEMPLATED_G_MARGINAL_LOG_PRB);

    /* cleanup memory */
    free( mapped_seq );
    free_penalty_array( &pa );

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
            = r->subtemplates + (*current_mapping)->pos_in_template;

        recheck_candidate_mapping( *current_mapping, rst, genome, error_model );

        /* Add the log probability (equivalent to product) */
        rechecked_group_penalty += (*current_mapping)->penalty;

        /* Move pointer to the next candidate mapping in the group */
        current_mapping++;
    }

    /* Add the penalty for the fragment length *if* this read is a full fragment
       NOTE: this may not be correct for paired RNA-Seq */
    if( r->prior.frag_type == FULL_GENOME_FRAGMENT ||
        r->prior.frag_type == FULL_TRANSCRIPTOME_FRAGMENT )
    {
        int fl = get_fl_for_joined_candidate_mappings( group_start );
        rechecked_group_penalty += get_fl_log_prb( fl_dist, fl );
    }

    *penalty = rechecked_group_penalty;

    return;
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
filter_joined_mapping_penalties( 
    float** joined_mapping_penalties, int joined_mappings_len )
{
    /* Rewrite the penalties array with just the penalties from the filtered
     * mappings */
    int num_filtered_penalties = 0;

    int i;
    for( i = 0; i < joined_mappings_len; i++ )
    {
        if( (*joined_mapping_penalties)[i] != 1 ) {
            num_filtered_penalties++;
        }
    }

    float* filtered_penalties = malloc( sizeof(float)*num_filtered_penalties );
    int fp_index = 0;
    for( i = 0; i < joined_mappings_len; i++ )
    {
        if( (*joined_mapping_penalties)[i] != 1 )
        {
            filtered_penalties[fp_index] = (*joined_mapping_penalties)[i];
            fp_index++;
        }
    }

    /* Free the original penalties array and replace it with a pointer to the
     * new, filtered penalties array */
    free( *joined_mapping_penalties );
    *joined_mapping_penalties = filtered_penalties;

    return;
}

void
filter_joined_candidate_mappings( candidate_mapping*** joined_mappings, 
                                  float** joined_mapping_penalties,
                                  int* joined_mappings_len,
                                  
                                  struct genome_data* genome,
                                  struct read* r,
                                  struct error_model_t* error_model,
                                  struct fragment_length_dist_t* fl_dist,

                                  struct mapping_params* mapping_params
    )
{
    /* Unpack the search parameters for the recheck */
    float min_match_penalty = mapping_params->recheck_min_match_penalty;
    float max_penalty_spread = mapping_params->recheck_max_penalty_spread;

    /* Initialize the max penalty to the smallest allowable penalty */
    float max_penalty = min_match_penalty;

    /* Pointer to the start of the current set of joined candidate_mappings */
    candidate_mapping** current_mapping = *joined_mappings;

    int i;
    /* recheck the locations, updating joined_mapping_penalties, 
       and find the max penalty */
    for( i = 0; i < *joined_mappings_len; i++ )
    {
        recheck_joined_candidate_mappings( 
            current_mapping, (*joined_mapping_penalties) + i,
            genome, r, error_model, fl_dist 
        );

        /* we may need to update the max penalty */
        if( (*joined_mapping_penalties)[i] > max_penalty ) {
            max_penalty = (*joined_mapping_penalties)[i];
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
            if( (*joined_mapping_penalties)[i] 
                < ( max_penalty - max_penalty_spread ) )
            {
                filter_current_group = true;
            }
        }

        /* Make sure that the penalty isn't too low */
        if( (*joined_mapping_penalties)[i] < min_match_penalty )
        {
            filter_current_group = true;
        }

        /* if either test failed, filter the current group */
        if( filter_current_group )
        {
            remove_candidate_mapping_group( current_mapping );
            /* Mark the invalid penalty with 1 - since the 
             * joined_mapping_penalties are log probabilities, any value > 0
             * is invalid. */
            (*joined_mapping_penalties)[i] = 1;
            filtered_mappings_len -= 1;
        }

        advance_pointer_to_start_of_next_joined_candidate_mappings(
                &current_mapping, i, *joined_mappings_len );
    }

    filter_joined_mapping_penalties( 
        joined_mapping_penalties, *joined_mappings_len );

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

    /* move the current_mapping pointer to the start of the first set of joined
     * candidate_mappings */
    while( *current_mapping == NULL )
        current_mapping++;

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
            /* add a sublocation for every M entry in a candidate mapping's
             * cigar string. These correspond to alignments around a gap in
             * a gapped assay. Since every candidate mapping has at least one
             * M entry, this code is correct for both gapped and ungapped
             * assays. */
            int j;
            for( j = 0; j < (*current_mapping)->cigar_len; j++ )
            {
                if( (*current_mapping)->cigar[j].op == 'M' )
                {
                    size += sizeof( mapped_read_sublocation );
                }
            }

            current_mapping++;
        }

        if( i < joined_mappings_len - 1 )
        {
            /* skip any NULL markers until we get to the start of the next set
             * of candidate_mappings */
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
populate_mapped_read_sublocations_from_candidate_mappings(
        candidate_mapping*** cm_ptr,
        char** rd_ptr )
{
    /* Get a pointer to the current candidate mapping */
    candidate_mapping* current_mapping = **cm_ptr;

    while( current_mapping != NULL )
    {
        /* Keep track of the start position of each mapped location (including
         * the cumulative offset from gaps) */
        int start_pos = current_mapping->start_bp;

        /* Add a sublocation for each M region in each candidate mappings in
         * the set of joined candidate mappings */
        int j;
        for( j = 0; j < current_mapping->cigar_len; j++ )
        {
            struct CIGAR_ENTRY cigar_entry = current_mapping->cigar[j];

            if( cigar_entry.op == 'M' )
            {
                mapped_read_sublocation* subloc = 
                    (mapped_read_sublocation*) *rd_ptr;

                subloc->start_pos = start_pos;
                subloc->length = cigar_entry.len;

                if( current_mapping->rd_strnd == FWD )
                {
                    subloc->rev_comp = 0;
                } else if( current_mapping->rd_strnd == BKWD ) {
                    subloc->rev_comp = 1;
                } else {
                    assert( current_mapping->rd_strnd == FWD ||
                            current_mapping->rd_strnd == BKWD );
                }

                /* TODO unused for now */
                subloc->is_full_contig = 0;

                /* Look ahead in the cigar string to set the
                 * next_subread_is_gapped is flag */
                if( j == current_mapping->cigar_len - 1 )
                {
                    /* if there is no next entry in the cigar string */
                    subloc->next_subread_is_gapped = 0;
                } else if ( current_mapping->cigar[j+1].op == 'N' ) {
                    /* if the next entry is a gap */
                    subloc->next_subread_is_gapped = 1;
                } else {
                    assert( false );
                }

                /* Look ahead in the set of joined candidate mappings to set
                 * the next_subread_is_ungapped flag */
                candidate_mapping* next_mapping = *(*cm_ptr + 1);
                if( next_mapping == NULL )
                {
                    subloc->next_subread_is_ungapped = 0;
                } else {
                    /* Only set this flag on the last subread in a set of
                     * gapped subreads */
                    if( j == current_mapping->cigar_len - 1 ) {
                        subloc->next_subread_is_ungapped = 1;
                    }
                }

                /* Advance pointer in mapped_read_t to the next sublocation to
                 * fill in */
                *rd_ptr += sizeof( mapped_read_sublocation );
            }

            /* Move the start position to the start of the next entry in the
             * cigar string */
            start_pos += cigar_entry.len;
        }

        /* Advance the candidate mapping pointers to the next candidate mapping
         * in the set of joined candidate mappings */
        (*cm_ptr)++;
        current_mapping = **cm_ptr;
    }

    return;
}

void
populate_mapped_read_locations_from_joined_candidate_mappings(
        mapped_read_t** rd,
        candidate_mapping** joined_mappings,
        int joined_mappings_len,
        float* joined_mapping_penalties )
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

        prologue->chr = (*current_mapping)->chr;

        /* The strand of a mapped_read_location (which may be built from
         * multiple candidate mappings with different strands) is determined by
         * the first candidate mapping (first read) in the joined set. */
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

        /* Set the mapped_read_location's trimmed_length from the first
         * candidate mapping
         * (corresponding to the first read in the fragment) */
        assert( (*current_mapping)->trimmed_length >= 0 );
        assert( (*current_mapping)->trimmed_length <= TRIMMED_LENGTH_MAX );
        prologue->trimmed_length = (*current_mapping)->trimmed_length;

        /* unused_bits are initialized to 0 by the calloc in
         * init_new_mapped_read_from_single_read_id */

        /* Convert the sum of the log joined_mapping_penalties to a probability
         * in standard [0,1] probability space */
        prologue->seq_error = pow( 10, joined_mapping_penalties[i] );
        assert( prologue->seq_error >= 0 && prologue->seq_error <= 1 );

        rd_ptr += sizeof( mapped_read_location_prologue );

        populate_mapped_read_sublocations_from_candidate_mappings(
                &current_mapping, &rd_ptr );

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
        int joined_mappings_len,
        float* joined_mapping_penalties )
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
           &rd, joined_mappings, joined_mappings_len, joined_mapping_penalties);
    
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
        statmap_log( LOG_FATAL, "Could not open mapped reads file %s",  fname  );
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
        statmap_log( LOG_FATAL, "Could not open unmappable reads file %s",  fname_buffer  );
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
        statmap_log( LOG_FATAL, "Could not open nonmapping reads file %s",  fname_buffer );
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
        statmap_log( LOG_FATAL, "Unrecognized mode '%c' for open mapped reads db.",  (*rdb)->mode  );
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

    if( (*rdb)->fl_dist != NULL ) {
        free_fl_dist( &((*rdb)->fl_dist) );
    }
    
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
        statmap_log( LOG_ERROR, "Mapped Reads DB is read-only - cannot add read." );
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
        statmap_log( LOG_FATAL, "Error writing to packed mapped reads db." );
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
        statmap_log( LOG_FATAL, "Can only rewind mapped reads db that is open for reading ( mode 'r' )." );
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
        statmap_log( LOG_FATAL, "Cannot get read from mapped reads db unless it is open for reading." );
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
    assert( mpd_rds_db != NULL );
    struct fragment_length_dist_t* fl_dist = mpd_rds_db->fl_dist;
    
    /* build an index for this mapped_read */
    mapped_read_index* rd_index = build_mapped_read_index( rd );

    if( 0 == rd_index->num_mappings )
        return;
    
    float *prbs = calloc( rd_index->num_mappings, sizeof(float) );
    
    /* prevent divide by zero */
    double prb_sum = ML_PRB_MIN;
    MPD_RD_ID_T i;
    for( i = 0; i < rd_index->num_mappings; i++ )
    {
        mapped_read_location* loc = rd_index->mappings[i];

        float cond_prob = get_seq_error_from_mapped_read_location( loc );

        if( mapped_read_location_is_paired( loc) )
            cond_prob *= get_fl_prb( fl_dist, get_fl_from_mapped_read_location( loc ) );

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
    if( fp_size == 0 ) {
        statmap_log( LOG_FATAL, "Cannot mmap empty mapped reads db" );
    }
    rewind( rdb->mapped_fp ); // reset fp

    /* find the size of the opened file */
    struct stat buf;
    fstat(fdin, &buf);
    rdb->mmapped_data_size = buf.st_size;
    
    #ifdef MALLOC_READS_DB
    fseek( rdb->mapped_fp, 0, SEEK_SET );

    statmap_log( LOG_NOTICE, "Allocating %zu bytes for the mapped reads db.", buf.st_size );
    rdb->mmapped_data = malloc( buf.st_size );
    if( NULL == rdb->mmapped_data ) {
        statmap_log( LOG_FATAL, "Failed to allocate %zu bytes for the mapped reads.", (size_t) buf.st_size );
    }
    
    fread( rdb->mmapped_data, buf.st_size, 1, rdb->mapped_fp );
           
    #else
    /* mmap the file */
    rdb->mmapped_data
        = mmap( NULL, rdb->mmapped_data_size,  
                PROT_READ|PROT_WRITE, 
                MAP_POPULATE|MAP_SHARED, fdin, (off_t) 0 );

    if( rdb->mmapped_data == (void*) -1 ) {
        statmap_log( LOG_FATAL, "Can not mmap the fdescriptor '%i'", fdin );
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
    if( error != 0 ) {
        statmap_log( LOG_FATAL, "Could not munmap mapped reads db (returned %i)", error );
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
    if( num_indexed_reads != rdb->num_mapped_reads ) {
        statmap_log( LOG_FATAL,
                "The number of indexed reads (%i) is not equal to the number of reads in the mapped read db ( %i). This may indicate that the mapped read db is corrupt.",
                num_indexed_reads, rdb->num_mapped_reads
            );
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
        mapped_read_index* rd_index = build_mapped_read_index( mapped_rd );

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

