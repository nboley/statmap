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

/****** Temporary containers ******/

/*
 * Mapped read sublocations container is a utility data strucutre to add in
 * building mapped read locations
 */

void
init_mapped_read_sublocations_container(
        mapped_read_sublocations_container** c )
{
    (*c) = malloc( sizeof( mapped_read_sublocations_container ) );

    (*c)->container = NULL;
    (*c)->length = 0;
}

void
free_mapped_read_sublocations_container(
        mapped_read_sublocations_container* c )
{
    if( c == NULL ) return;

    if( c->container != NULL )
        free( c->container );

    free(c);

    return;
}

void
add_mapped_read_sublocation_to_container(
        mapped_read_sublocation* sub_loc,
        mapped_read_sublocations_container* c )
{
    c->length += 1;
    c->container = realloc( c->container,
            sizeof( mapped_read_sublocation ) * c->length );
    c->container[c->length-1] = *sub_loc;
}


/*************************************************************************
 *
 *  Mapped Read Locations
 *
 */

void
build_mapped_read_location( mapped_read_location** loc,
                            MRL_CHR_TYPE chr,
                            MRL_FLAG_TYPE flag,
                            ML_PRB_TYPE seq_error,
                            mapped_read_sublocations_container* sublocations )
{
    assert( sublocations->length > 0 );

    /* Start by allocating space for the prologue and the sublocations */
    *loc = malloc( sizeof( mapped_read_location_prologue ) +
                   ( sizeof( mapped_read_sublocation ) * sublocations->length ));

    mapped_read_location_prologue *prologue =
        (mapped_read_location_prologue *) *loc;

    assert( chr < CHR_NUM_MAX );
    prologue->chr = chr;
    prologue->flag = flag;
    prologue->seq_error = seq_error;
    /* we don't know if there are more mapped read locations, so this defaults
     * to 0. This is updated in build_mapped_read */
    prologue->are_more = 0;

    mapped_read_sublocation *current_sublocation =
        (mapped_read_sublocation *)
        ( (char*)(*loc) + sizeof(mapped_read_location_prologue) );

    int i;
    for( i = 0; i < sublocations->length; i++ )
    {
        *current_sublocation = sublocations->container[i];
        current_sublocation++; // pointer arithmetic
    }
}

void
free_mapped_read_location( mapped_read_location* loc )
{
    if( loc == NULL ) return;

    free( loc );
    return;
}

/***** Temporary container for mapped_read_location(s) *****/

void
init_mapped_read_locations_container(
        mapped_read_locations_container **c )
{
    (*c) = malloc( sizeof( mapped_read_locations_container ));

    (*c)->container = NULL;
    (*c)->length = 0;
}

void
free_mapped_read_locations_container(
        mapped_read_locations_container* c )
{
    if( c == NULL ) return;

    if( c->container != NULL )
    {
        /* free the mapped read locations */
        int i;
        for( i = 0; i < c->length; i++ )
        {
            free( c->container[i] ); // & (?)
        }

        free( c->container );
    }

    free(c);

    return;
}

void
add_mapped_read_location_to_container(
        mapped_read_location* loc,
        mapped_read_locations_container* c )
{
    c->length += 1;
    c->container = realloc( c->container,
            sizeof( mapped_read_location* ) * c->length );
    c->container[c->length-1] = loc;
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

size_t
get_size_of_mapped_read_locations_container(
        mapped_read_locations_container* locs )
{
    size_t total_size = 0;

    int i;
    /* loop over each mapped_read_location in the container */
    for( i = 0; i < locs->length; i++ )
    {
        mapped_read_location* loc = locs->container[i];
        size_t loc_size = get_size_of_mapped_read_location( loc );
        total_size += loc_size;
    }

    return total_size;
}

void
build_mapped_read( mapped_read_t** rd,
                   int* read_ids,
                   int num_read_ids,
                   mapped_read_locations_container* locs )
{
    assert( num_read_ids > 0 );
    assert( locs->length > 0 );

    /* Allocate memory for the read_id_nodes and mapped_read_locations */
    size_t read_id_nodes_bytes = sizeof( read_id_node ) * num_read_ids;
    size_t mapped_read_locations_bytes =
        get_size_of_mapped_read_locations_container( locs );
    *rd = malloc( read_id_nodes_bytes + mapped_read_locations_bytes );

    /* Get a pointer into the mapped_read_t to iterate over the bytes */
    char* ptr = (char*) *rd;

    /* Copy the read id nodes into the pseudo structure */
    int i;
    for( i = 0; i < num_read_ids; i++ )
    {
        read_id_node* curr_node = (read_id_node*) ptr;

        curr_node->read_id = read_ids[i];

        /* If this is the last read_id_node, set the are_more flag to 0 */
        if( i == num_read_ids - 1 )
        {
            curr_node->are_more = 0;
        }
        else /* Otherwise, set it to 1 */
        {
            curr_node->are_more = 1;
        }

        ptr += sizeof( read_id_node );
    }

    /* Copy the mapped_read_locations into the pseudo structure */
    for( i = 0; i < locs->length; i++ )
    {
        mapped_read_location* loc = locs->container[i];
        size_t loc_size = get_size_of_mapped_read_location( loc );
        memcpy( ptr, loc, loc_size );

        /* Set the are_more flag in each mapped_read_location */
        mapped_read_location_prologue* prologue =
            (mapped_read_location_prologue*) ptr;

        if( i == locs->length - 1 )
        {
            /* if this is the last mapped_read_location, set the are_more flag to 0 */
            prologue->are_more = 0;
        } else {
            /* otherwise, set it to 1 */
            prologue->are_more = 1;
        }

        ptr += loc_size;
    }

    return;
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
    index->mappings = malloc( num_allocated_locations*
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

/*************************************************************************
 *
 *  Convert Candidate Mappings to Mapped Reads
 *
 */

/* returns 1 for success, 0 for failure */
static inline int
convert_unpaired_candidate_mapping_into_mapped_read( 
    candidate_mapping* cm,
    mapped_read_location** loc    
)
{
    /* Ensure all of the flags are turned off */
    MRL_FLAG_TYPE flag = 0;
 
    /* Add the location */
    if( BKWD == cm->rd_strnd )
        flag |= FIRST_READ_WAS_REV_COMPLEMENTED;
    
    float seq_error = pow( 10, cm->penalty );
    assert( seq_error >= 0.0 && seq_error <= 1.0 );

    /* Build a single sublocation */
    mapped_read_sublocation subloc;

    subloc.start_pos = cm->start_bp;
    subloc.length = cm->rd_len;

    if( cm->rd_strnd == FWD ) {
        subloc.strand = 0;
    } else if( cm->rd_strnd == BKWD ) {
        subloc.strand = 1;
    } else {
        assert( cm->rd_strnd == FWD || cm->rd_strnd == BKWD );
    }

    subloc.next_subread_is_gapped = 0;
    subloc.next_subread_is_ungapped = 0;

    /* Build a container for the subloc */

    mapped_read_sublocations_container* sublocs = NULL;
    init_mapped_read_sublocations_container( &sublocs );
    add_mapped_read_sublocation_to_container( &subloc, sublocs );

    /* Build the mapped read location with the container of sublocations */
    build_mapped_read_location( loc, cm->chr, flag, seq_error, sublocs );

    free_mapped_read_sublocations_container( sublocs );
    
    return 1;
}


/* returns 1 for success, 0 for failure */
static inline int
join_two_candidate_mappings( 
    candidate_mapping* first_read,
    candidate_mapping* second_read,
    mapped_read_location** loc    
)
{
    /* Ensure all of the flags are turned off */
    MRL_FLAG_TYPE flag = IS_PAIRED;
    
    /* Set the appropriate flags */
    if( first_read->start_bp < second_read->start_bp )
        flag |= FIRST_PAIR_IS_FIRST_IN_GENOME;
    
    if( FWD == first_read->rd_strnd ) {
        // TODO this was the original code - is it correct?
        flag |= FIRST_READ_WAS_REV_COMPLEMENTED;
    }

    /* Get the chr */
    assert( first_read->chr == second_read->chr);
    int chr = first_read->chr;
    
    /* Get the start and stop */
    int start, stop;
    if( first_read->start_bp < second_read->start_bp )
    {
        start = first_read->start_bp;
        stop = second_read->start_bp + second_read->rd_len;
    } else {
        // TODO double check correctness
        start = second_read->start_bp;
        stop = first_read->start_bp + first_read->rd_len;
    }

    float seq_error = pow( 10, first_read->penalty + second_read->penalty );

    /* Build a single sublocation (TODO for now) */
    mapped_read_sublocation subloc;

    subloc.start_pos = start;
    subloc.length = stop - start;

    if( first_read->rd_strnd == FWD )
    {
        subloc.strand = 0;
    } else if( first_read->rd_strnd == BKWD ) {
        subloc.strand = 1;
    } else {
        assert( first_read->rd_strnd == FWD || first_read->rd_strnd == BKWD );
    }

    subloc.next_subread_is_gapped = 0;
    subloc.next_subread_is_ungapped = 0;

    /* Build a container for the subloc */

    mapped_read_sublocations_container* sublocs = NULL;
    init_mapped_read_sublocations_container( &sublocs );
    add_mapped_read_sublocation_to_container( &subloc, sublocs );

    /* Build the mapped read location with the container of sublocations */
    build_mapped_read_location( loc, chr, flag, seq_error, sublocs );

    free_mapped_read_sublocations_container( sublocs );
    
    /* ignore reads with zero probability ( possible with FL dist ) */
    if( seq_error > 2*ML_PRB_MIN ) {
        return 1;
    }

    return 0;
}

static inline double
mapped_read_from_candidate_mapping_arrays( 
    candidate_mapping* r1_array,
    int r1_array_len,
    candidate_mapping* r2_array,
    int r2_array_len,
    mapped_read_locations_container* locs )
{
    double prob_sum = 0;
    
    int i;
    for( i = 0; i < r1_array_len; i++ )
    {
        int j;
        for( j = 0; j < r2_array_len; j++ )
        {
            /* if the chrs mismatch, since these are sorted, we know that
             * there is no need to continue */
            if( r2_array[j].chr > r1_array[i].chr )
                break;
            
            if( r2_array[j].chr != r1_array[i].chr )
                continue;
            
            /* Determine which of the candidate mappings corresponds 
               with the first pair */
            candidate_mapping* first_read = NULL;
            candidate_mapping* second_read = NULL;
            if( PAIRED_END_1 == r1_array[i].rd_type ) {
                first_read = r1_array + i;
                second_read = r2_array + j;
            } else {
                assert( PAIRED_END_1 == r2_array[j].rd_type );
                first_read = r2_array + j;
                second_read = r1_array + i;
            }
            
            assert( first_read->chr == second_read->chr );
            
            mapped_read_location* loc = NULL;
            int rv = join_two_candidate_mappings( first_read, second_read, &loc );
            /* if the location is valid ( non-zero probability ) */
            if( rv == 1 )
            {
                prob_sum += get_seq_error_from_mapped_read_location( loc );
                add_mapped_read_location_to_container( loc, locs );
            }
        }
    }
    
    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_locations_from_paired_candidate_mappings( 
        mapped_read_locations_container* locs,
        candidate_mappings* mappings
    )
{
    double prob_sum = 0;
    
    assert( mappings->length > 0 );
    assert( mappings->mappings[0].rd_type > SINGLE_END ); 

    int p2_start = -1;

    /* 
     * Find the start of the pair 2 mapped reads
     */
    int i;
    for( i=0; i < mappings->length; i++ )
    {
        if( mappings->mappings[i].rd_type == PAIRED_END_2 )
        {
            p2_start = i;
            break;
        }
    }
    
    /* If there were no paired end 2 reads, then we say no reads mapped */
    if( -1 == p2_start ) {
        /* make sure that the number of reads is 0 */
        //assert( mpd_rd->num_mappings == 0);
        return 0;
    }

    prob_sum += mapped_read_from_candidate_mapping_arrays(
        mappings->mappings,
        p2_start,
        mappings->mappings + p2_start,
        mappings->length - p2_start,
        locs
    );

    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
double
build_mapped_read_locations_from_unpaired_candidate_mappings( 
        mapped_read_locations_container* locs,
        candidate_mappings* mappings )
{
    double prob_sum = 0;

    int i;    
    for( i = 0; i < mappings->length; i++ )
    {
        /* we expect a read to be either paired end or not - never both */
        assert( mappings->mappings[i].rd_type == SINGLE_END );

        /* if the mapping hasnt been determined to be valid, ignore it */
        if( (mappings->mappings)[i].recheck != VALID )
            continue;

        /* store local read location data */
        mapped_read_location* loc;
    
        int rv = 
            convert_unpaired_candidate_mapping_into_mapped_read( 
                mappings->mappings + i, &loc );
        
        /* if the conversion succeeded */
        if( 1 == rv )
        {
            prob_sum += get_seq_error_from_mapped_read_location( loc );
            add_mapped_read_location_to_container( loc, locs );
        }
    }

    return prob_sum;
}

void
build_mapped_read_from_candidate_mappings( 
        mapped_read_t** mpd_rd,
        candidate_mappings* mappings, 
        MPD_RD_ID_T read_id
    )
{
    /* 
     * Building mapped reads has several components:
     * First, the normal reads
     * Second, the paired end reads
     *    Paired end read can match iff:
     *    1) They are from the same read
     *    2) They are on opposite strands
     *    3) They are on the same chromosome
     *
     *    Luckily, every passed candidate is from the same
     *    read, so (1) is established automatically. Second, 
     *    the candidates are sorted ( in cmp_candidate_mappings )
     *    by read_type, strand, chr, bp_position. So, the plan is
     *
     * 1) Get and print all of the single ended reads ( rd_type == SINGLE_END )
     * 2) Find the Start of the PAIRED_END_1 and PAIRED_END_2 reads
     * 3) Merge join them 
     * 
     */
    
    if( mappings->length == 0 )
        return;
    
    /* store the sum of the marginal probabilities */
    double prob_sum;

    /* store the mapped locations */
    mapped_read_locations_container* mapped_locs = NULL;
    init_mapped_read_locations_container( &mapped_locs );

    /* this assumes all reads are either paired or not */
    if( mappings->mappings[0].rd_type == SINGLE_END )
    {
        prob_sum =
            build_mapped_read_locations_from_unpaired_candidate_mappings( 
                mapped_locs, mappings );
    } else {
        prob_sum =
            build_mapped_read_locations_from_paired_candidate_mappings( 
                mapped_locs, mappings );
    }    
    
    /* If there are no proper mappings ( this can happen if, 
       for instance, only one pair maps ) free the read and return */
    if( 0 == prob_sum )
    {
        free_mapped_read( *mpd_rd );
        *mpd_rd = NULL;
        return;
    }

    int read_ids[] = { read_id };
    build_mapped_read( mpd_rd, read_ids, 1, mapped_locs );
    
    free_mapped_read_locations_container( mapped_locs );

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
        if( get_flag_from_mapped_read_location( loc )&IS_PAIRED )
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
            /* cast current location to read_id_node so it can be examined */
            read_id_node* node = (read_id_node*) ptr;

            /* add the read_id_node to the index */

            /* if the array is full */
            if( num_indexed_reads + 1 == num_allcd_reads )
            {
                num_allcd_reads += REALLOC_BLOCK_SIZE;
                rdb->index = realloc(rdb->index,
                        num_allcd_reads*sizeof(struct mapped_reads_db_index_t) );
            }
            assert( num_indexed_reads < num_allcd_reads );

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

