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
 *  Mapped Read Locations
 *
 */

enum bool
first_sublocation_is_rev_complemented( mapped_read_location* loc )
{
    char* ptr = (char*) loc;

    /* Skip the prologue */
    ptr += sizeof( mapped_read_location_prologue );

    /* Get the first mapped_read_sublocation (assumes there's at least one) */
    mapped_read_sublocation* first_sublocation =
        (mapped_read_sublocation *) ptr;

    /* enum bool and single bit flags use the same encoding, so we can return
     * the value of the flag directly */
    return first_sublocation->rev_comp;
}

void
free_mapped_read_location( mapped_read_location* loc )
{
    if( loc == NULL ) return;

    free( loc );
    return;
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

void
join_candidate_mappings( candidate_mappings* mappings, 
                         candidate_mapping*** joined_mappings, 
                         float** penalty,
                         int* joined_mappings_len
    )
{
    /* joined_mappings is a list of pointers to candidate_mappings. Each
     * "joined mapping" is a list of pointers to candidate_mappings separated
     * by NULL poiners */

    /* TODO just handle single end reads for now */
    /* just loop through the mappings and add everything up */

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
        *penalty = realloc( *penalty,
                sizeof(float) * *joined_mappings_len );
        (*penalty)[*joined_mappings_len - 1] = mappings->mappings[i].penalty;
    }

    return;
}

void
filter_joined_candidate_mappings( candidate_mappings* mappings, 
                                  candidate_mapping*** joined_mappings, 
                                  int* joined_mappings_len
    )
{
    /* For now, just don't filter */
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

        /* advance beyond the NULL marker to the start of the next set of
         * joined candidate mappings */
        current_mapping++;
    }

    return size;
}

void
init_new_mapped_read_from_single_read_id( mapped_read_t** rd,
                                          MPD_RD_ID_T read_id,
                                          size_t mapped_read_size )
{
    /* Allocate memory for the mapped read */
    *rd = malloc( mapped_read_size );

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
         * are log probabilites, we can just multiply them, and do the
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

            /* look ahead to figure out if the next subread is gapped or not.
             * TODO: how to distinguish gapped vs. ungapped (probably metadata
             * on the cm) */
            if( *(current_mapping + 1) == NULL )
            {
                subloc->next_subread_is_gapped = 0;
                subloc->next_subread_is_ungapped = 0;
            } else {
                /* TODO - for now, only handle paired end case. */
                subloc->next_subread_is_gapped = 1;
                subloc->next_subread_is_ungapped = 0;
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

        /* advance beyond the NULL marker to the start of the next set of
         * joined candidate mappings */
        current_mapping++;
    }
}

void
build_mapped_read_from_joined_candidate_mappings(
        mapped_read_t** rd,
        MPD_RD_ID_T read_id,
        candidate_mapping** joined_mappings,
        int joined_mappings_len )
{
    /* TODO for now, assume there is only one read_id_node */
    size_t mapped_read_size = sizeof( read_id_node );
    mapped_read_size +=
        calculate_mapped_read_space_from_joined_candidate_mappings(
                joined_mappings, joined_mappings_len );

    init_new_mapped_read_from_single_read_id( rd, read_id, mapped_read_size );

    populate_mapped_read_locations_from_joined_candidate_mappings(
            rd, joined_mappings, joined_mappings_len );
}

void
build_mapped_read_from_candidate_mappings( 
        mapped_read_t** mpd_rd,
        candidate_mappings* mappings, 
        MPD_RD_ID_T read_id
    )
{
    /* 
     * 1) Join candidate mappings ( read subtemplates ) from the same 
     *    read template
     *    - join_candidate_mappings
     *      produces a list of pointers into 'mappings', where mappings
     *      that shouldn't be joined are seperated by NULL pointers
     *    - filter_joined_candidate_mappings
     *      calculate error estiamtes, and filter bad joins
     * 2) For each joined candidate mapping
     *    - Init a mapped read ( and allocate memory )
     *      - calculate_mapped_read_space_from_joined_candidate_mappings
     *      - init_new_mapped_read_from_single_read_id
     *          ( init the new mapped read, and add in the read id )
     *    - Populate the new mapped read
     *         FOR EACH set_of_joined_candidate_mappings in joined_mappings
     *             - build_mapped_read_location_prologue
     *             ( add this into mapped_read_t )
     *             FOR EACH candidate mapping in set_of_joined_candidate_mappings
     *                 - build_mapped_read_sublocation_from_candidate_mapping
     *                 ( add this into mapped_read_t )
     *               
     * 3) Fix the SAM code
     *    - almost trivial, because all of the information is inside of the 
     *      mapped eads structure
     *
     */

    /* DEBUG */
    fprintf( stderr, "READ_ID %u has %i candidate mappings\n", read_id, mappings->length );

    /* Build a list of joined candidate mappings */
    int joined_mappings_len = 0;
    candidate_mapping** joined_mappings = NULL;
    float* penalties = NULL;

    join_candidate_mappings( mappings,
                             &joined_mappings,
                             &penalties,
                             &joined_mappings_len );

    /* TODO stub */
    filter_joined_candidate_mappings( mappings,
                                      &joined_mappings,
                                      &joined_mappings_len );

    build_mapped_read_from_joined_candidate_mappings( mpd_rd,
                                                      read_id,
                                                      joined_mappings,
                                                      joined_mappings_len );

    free( joined_mappings );
    free( penalties );

    return;
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
    (*rdb)->fp = fopen( fname, mode );
    if( (*rdb)->fp == NULL )
    {
        perror("FATAL       :  Could not open mapped reads file");
        assert( false );
        exit(-1);
    }
    
    /* number of mapped reads in the DB */
    (*rdb)->num_mapped_reads = 0;
    
    /* Initialize the mode to 0, this will be set by
       the mode specific init function */
    (*rdb)->mode = 0;

    /* use a mutex to eliminate the spurious helgrind warnings that we got
       when using a spinlock */
    pthread_mutexattr_t mta;
    pthread_mutexattr_init(&mta);
    (*rdb)->mutex = malloc( sizeof(pthread_mutex_t) );
    pthread_mutex_init( (*rdb)->mutex, &mta );
    
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

    /* write placeholder for size of mapped reads db, this 
       will be updated when we close the mapped read db*/
    MPD_RD_ID_T placeholder = 0;
    fwrite( &placeholder, sizeof(MPD_RD_ID_T), 1, (*rdb)->fp );
    
    (*rdb)->mode = 'w';
}

void
build_fl_dist_from_file( struct mapped_reads_db* rdb, FILE* fl_fp )
{
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    return;
}

void
build_fl_dist_from_filename( struct mapped_reads_db* rdb, char* filename )
{
    FILE* fl_fp = fopen( filename, "r" );
    if( fl_fp == NULL )
    {
        fprintf( stderr, "Failed to open fl_dist from filename %s\n", filename);
        exit(-1);
    }
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    fclose( fl_fp );
}

void
close_reading_specific_portions_of_mapped_reads_db( struct mapped_reads_db* rdb)
{
    munmap_mapped_reads_db( rdb );

    if( rdb->fl_dist != NULL ) {
        free_fl_dist( &(rdb->fl_dist) );
        rdb->fl_dist = NULL;
    }

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
    fseek( rdb->fp, 0, SEEK_SET );
    fwrite( &(rdb->num_mapped_reads), sizeof(MPD_RD_ID_T), 1, rdb->fp );
    
    fclose( rdb->fp );

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
        perror( "Unrecognized mode for open mapped read db." );
        assert( false );
        exit( -1 );
    }
    
    pthread_mutex_destroy( (*rdb)->mutex );
    free( (*rdb)->mutex );
    
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

    pthread_mutex_lock( rdb->mutex );
    error = write_mapped_read_to_file( rd, rdb->fp );
    rdb->num_mapped_reads += 1;
    pthread_mutex_unlock( rdb->mutex );

    if( error < 0 )
    {
        fprintf(stderr, "FATAL       :  Error writing to packed mapped reads db.\n");
        assert( false );
        exit( -1 );
    }
    
    return;
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
    pthread_mutex_lock( rdb->mutex );
    /* if we have read every read */
    if( rdb->current_read == rdb->num_mapped_reads )
    {
        pthread_mutex_unlock( rdb->mutex );
        *rd = NULL;
        return EOF;
    }
    
    MPD_RD_ID_T current_read_id = rdb->current_read;
    rdb->current_read += 1;
    pthread_mutex_unlock( rdb->mutex );

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
    int fdin = fileno( rdb->fp );
    
    /* make sure the entire file has been written to disk */
    fflush( rdb->fp );

    /* check that the file is not empty before trying to mmap */
    fseek(rdb->fp, 0L, SEEK_END);   // seek to end
    long fp_size = ftell(rdb->fp);  // tell() to get size
    if( fp_size == 0 )
    {
        fprintf( stderr,
                 "FATAL       :  Cannot mmap empty mapped reads db.\n" );
        exit( 1 );
    }
    rewind( rdb->fp ); // reset fp

    /* find the size of the opened file */
    struct stat buf;
    fstat(fdin, &buf);
    rdb->mmapped_data_size = buf.st_size;
    
    #ifdef MALLOC_READS_DB
    fseek( rdb->fp, 0, SEEK_SET );

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
    
    fread( rdb->mmapped_data, buf.st_size, 1, rdb->fp );
           
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

