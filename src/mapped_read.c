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

size_t
get_num_allocated_bytes_in_mapped_read_location( mapped_read_location* loc )
{
    mapped_read_location_prologue *prologue =
        (mapped_read_location_prologue *) loc;
    return prologue->num_allocated_bytes;
}

void
init_mapped_read_location( mapped_read_location** loc,
                           MRL_CHR_TYPE chr,
                           MRL_FLAG_TYPE flag,
                           ML_PRB_TYPE seq_error,
                           enum bool are_more )
{
    /* Start by allocating space for the prologue */
    *loc = malloc( sizeof( mapped_read_location_prologue ));

    /* Cast loc to a pointer to prologue so we can set it up */
    mapped_read_location_prologue *prologue =
        (mapped_read_location_prologue *) *loc;

    prologue->num_allocated_bytes = sizeof( mapped_read_location_prologue );

    assert( chr < CHR_NUM_MAX );
    prologue->chr = chr;

    prologue->flag = flag;
    prologue->seq_error = seq_error;

    /* enum bool: false = 0, true = 1, so we can assign to an unsigned
     * bitfield of length 1 */
    prologue->are_more = are_more; 
}

void
free_mapped_read_location( mapped_read_location* loc )
{
    if( loc == NULL ) return;

    free( loc );
    return;
}

// TODO store (start, stop) or (start, length) ?
// Originally stored start, stop
void
add_subtemplate_location_to_mapped_read_location( mapped_read_location* loc,
                                                  MRL_START_POS_TYPE start,
                                                  MRL_FL_TYPE length,
                                                  enum STRAND strand,
                                                  enum bool are_more )
{
    /* Assumes the prologue has already been allocated and assigned */
    mapped_read_location_prologue *prologue =
        (mapped_read_location_prologue *) loc;

    size_t original_size = prologue->num_allocated_bytes;
    size_t new_size = original_size + sizeof( mapped_read_subtemplate_location );

    loc = realloc( loc, new_size );
    if( loc == NULL )
    {
        fprintf( stderr, "Error allocating memory for mapped read location.\n");
        assert(false);
        exit(-1);
    }

    mapped_read_subtemplate_location* new_st_loc = 
        (mapped_read_subtemplate_location *) ( (char*)loc + original_size );

    new_st_loc->start_pos = start;
    new_st_loc->length = length;

    if( strand == FWD ) {
        new_st_loc->strand = 0;
    } else if ( strand == BKWD ) {
        new_st_loc->strand = 1;
    } else {
        assert( strand == FWD || strand == BKWD );
    }

    new_st_loc->are_more = are_more;
}


/*************************************************************************
 *
 *  Mapped Reads
 * 
 *  Reads that have been joined, but unlike mapped reads proper have not
 *  had the 'extra' read information attached.
 *
 */

void
set_num_allocated_bytes_in_mapped_read( size_t size,
                                        mapped_read_t* rd )
{
    *( (size_t*) rd ) = size;

    return;
}

size_t
get_num_allocated_bytes_in_mapped_read( mapped_read_t* rd )
{
    return *( (size_t*) rd );
}

void
realloc_mapped_read( mapped_read_t* rd, size_t size )
{
    rd = realloc( rd, size );
    if( rd == NULL && size > 0 )
    {
        fprintf( stderr, "FATAL       :  Error allocating memory for mapped read.\n");
        assert(false);
        exit(-1);
    }

    set_num_allocated_bytes_in_mapped_read( size, rd );
}

void
init_mapped_read( mapped_read_t** rd )
{
    *rd = NULL;
    /* allocate enough memory to store num_allocated_bytes */
    size_t alloc_size = sizeof( size_t );
    realloc_mapped_read( *rd, alloc_size );

    return;
}

void
free_mapped_read( mapped_read_t* rd )
{
    free( rd );
    return;
}

void
add_read_id_node_to_mapped_read( MPD_RD_ID_T read_id,
                                 enum bool are_more,
                                 mapped_read_t* rd )
{
    /* IMPORTANT: assumes that we are building the mapped_read_t, and no
     * mapped_read_locations have been added yet */

    /* reallocate memory to store the new read id node */
    size_t original_size = get_num_allocated_bytes_in_mapped_read( rd );
    size_t new_size = original_size + sizeof( read_id_node );
    realloc_mapped_read( rd, new_size );

    read_id_node *node = (read_id_node*) ((char*)rd + original_size);

    /* Since we are using 31 bits to store the read_id, but READ_ID_TYPE can't
     * be 31 bits in size, we check the value to make sure it will fit in the
     * packed structure */
    assert( read_id < MAX_READ_ID );

    /* Set the read id */
    node->read_id = read_id;

    // TODO for now, every mapped_read should only have 1 read_id_node
    assert( !are_more );

    /* Set the are_more flag */
    if( are_more )
    {
        node->are_more = 1;
    } else {
        node->are_more = 0;
    }
}

MPD_RD_ID_T
get_read_id_from_mapped_read( mapped_read_t* rd )
{
    /* TODO for now we assume each mapped_read_t has only one read id node.
     * Therefore, we just take the first read id node and return the stored
     * read id */
    char* rd_ptr = (char*) rd;

    /* skip num_allocated_bytes */
    rd_ptr += sizeof( size_t );

    /* cast to read_id_node pointer to dereference */
    read_id_node *node = (read_id_node*) rd_ptr;
    return node->read_id;
}

void
add_location_to_mapped_read( mapped_read_location* loc,
                             mapped_read_t* rd )
{
    /* Allocate new space for the location */
    size_t original_size = get_num_allocated_bytes_in_mapped_read( rd );
    size_t loc_size = get_num_allocated_bytes_in_mapped_read_location( loc );
    size_t new_size = original_size + loc_size;
    realloc_mapped_read( rd, new_size );

    /* Pointer to start of new mapped_read_location in mapped_read_t */
    mapped_read_location* new_loc =
        (mapped_read_location *) ( (char*) rd + original_size );

    /* Direct copy bytes from mapped_read_location into rd */
    memcpy( new_loc, loc, loc_size );
}

size_t 
write_mapped_read_to_file( mapped_read_t* read, FILE* of  )
{
    size_t num_written = 0;
    size_t num_allocated = get_num_allocated_bytes_in_mapped_read( read );

    num_written = fwrite( read, sizeof(char), num_allocated, of );
    if( num_written != num_allocated )
        return -num_written;

    return 0;
}

char*
skip_subtemplate_locations_in_mapped_read_location( char* rd )
{
    /* assumes the rd pointer is positioned at the start of
     * a mapped_read_subtemplate_location */
    char* rd_ptr = rd;

    /* skip the subtemplate locations */
    while(true)
    {
        mapped_read_subtemplate_location *st_loc =
            (mapped_read_subtemplate_location *) rd_ptr;

        if( !(st_loc->are_more) )
        {
            rd_ptr += sizeof( mapped_read_subtemplate_location );
            break;
        }

        rd_ptr += sizeof( mapped_read_subtemplate_location );
    }

    return rd_ptr;
}

char*
skip_mapped_read_location_in_mapped_read_t( char* rd )
{
    /* assumes the rd pointer is positioned at the start of
     * a mapped_read_location_prologue */

    char* rd_ptr = rd;

    /* skip the prologue and all subtemplate location entries to get to the
     * start of the next mapped_read_location */
    rd_ptr += sizeof( mapped_read_location_prologue );

    /* skip the subtemplate locations */
    rd_ptr = skip_subtemplate_locations_in_mapped_read_location( rd_ptr );

    return rd_ptr;
}

char*
skip_mapped_read_locations_in_mapped_read_t( char* rd )
{
    while(true)
    {
        mapped_read_location_prologue* loc =
            (mapped_read_location_prologue*) rd;

        if( !(loc->are_more) )
        {
            /* skip the current (final) location */
            rd = skip_mapped_read_location_in_mapped_read_t( rd );
            break;
        }

        rd = skip_mapped_read_location_in_mapped_read_t( rd );
    }

    return rd;
}

char*
get_start_of_mapped_read_locations_in_mapped_read_t( mapped_read_t* rd )
{
    /* Assumes rd is at the start of a mapped_read_t structure */

    /* Cast the mapped_read_t to a char* in order to iterate over the
     * bytepacked data structures using sizeof */
    char* rd_ptr = (char*) rd;

    /* skip num_allocated_bytes */
    rd_ptr += sizeof( size_t );

    /* skip read_id_node(s) */
    while(true)
    {
        /* cast current location to read_id_node so it can be examined */
        read_id_node* node = (read_id_node*) rd_ptr;

        if( !(node->are_more) )
        {
            /* increment pointer to the start of the mapped read locations,
             * and break */
            rd_ptr += sizeof( read_id_node );
            break;
        }

        rd_ptr += sizeof( read_id_node );
    }

    return rd_ptr;
}

/*
 * NOTE - this code assumes there is always at least one read_id_node and
 * one mapped_read_location in the pseudo structure
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

    /* Locate the start of the mapped read locations in this mapped read. Cast
     * the pointer to a char* in order to iterate over the bytepacked data
     * structure using sizeof */
    char* rd_ptr = get_start_of_mapped_read_locations_in_mapped_read_t( rd );

    int num_mapped_locations = 0;

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

        index->mappings[num_mapped_locations] = (mapped_read_location*) rd_ptr;
        num_mapped_locations += 1;

        /* cast read pointer to mapped_read_location_prologue */
        mapped_read_location_prologue* loc =
            (mapped_read_location_prologue *) rd_ptr;

        /* if this was the last mapped_read_location, we're done counting */
        if( !( loc->are_more ) )
        {
            break;
        }

        /* move rd_ptr to the start of the next mapped_read_location */
        rd_ptr = skip_mapped_read_location_in_mapped_read_t( rd_ptr );
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

    free( index->mappings );
    free( index );

    return;
}

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
    
    init_mapped_read_location( loc, cm->chr, flag, seq_error, false );
    add_subtemplate_location_to_mapped_read_location(
            *loc, cm->start_bp, cm->rd_len, cm->rd_strnd, false );
    
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
    
    /* Initialize the mapped read location */
    // TODO - are_more - we don't know - add code to add_mapped_read_location
    init_mapped_read_location( loc, chr, flag, seq_error, false );

    /* Add this location as the first read subtemplate location */
    // TODO are_more = false for now, pending candidate mappings rewrite
    add_subtemplate_location_to_mapped_read_location(
            *loc, start, stop-start, first_read->rd_strnd, false );
    
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
    mapped_read_t* mpd_rd
)
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
                add_location_to_mapped_read( loc, mpd_rd );
                prob_sum += get_seq_error_from_mapped_read_location( loc );
            }
        }
    }
    
    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_paired_candidate_mappings( 
        mapped_read_t* mpd_rd,
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
        mpd_rd
    );

    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_unpaired_candidate_mappings( 
        mapped_read_t* mpd_rd,
        candidate_mappings* mappings
    )
{
    double prob_sum = 0;

    /* store local read location data */
    mapped_read_location* loc;
    
    int i;    
    for( i = 0; i < mappings->length; i++ )
    {        
        /* we expect a read to be either paired end or not - never both */
        assert( mappings->mappings[i].rd_type == SINGLE_END );

        /* if the mapping hasnt been determined to be valid, ignore it */
        if( (mappings->mappings)[i].recheck != VALID )
            continue;

        int rv = 
            convert_unpaired_candidate_mapping_into_mapped_read( 
                mappings->mappings + i, &loc );
        
        /* if the conversion succeeded */
        if( 1 == rv )
        {
            prob_sum += get_seq_error_from_mapped_read_location( loc );
            add_location_to_mapped_read( loc, mpd_rd );
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
    
    /* Initialize the mapped read */
    init_mapped_read( mpd_rd );
    add_read_id_node_to_mapped_read( read_id, false, mpd_rd );

    if( mappings->length == 0 )
        return;
    
    /* store the sum of the marginal probabilities */
    double prob_sum;

    /* this assumes all reads are either paired or not */
    if( mappings->mappings[0].rd_type == SINGLE_END )
    {
        prob_sum = build_mapped_read_from_unpaired_candidate_mappings( 
            *mpd_rd, mappings );
    } else {
        prob_sum = build_mapped_read_from_paired_candidate_mappings( 
            *mpd_rd, mappings );
    }    
    
    /* If there are no proper mappings ( this can happen if, 
       for instance, only one pair maps ) free the read and return */
    if( 0 == prob_sum )
    {
        free_mapped_read( *mpd_rd );
        *mpd_rd = NULL;
        return;
    }
    
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
        free_mapped_read( rd );
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
        free_mapped_read( r );
    }
    
    free_mapped_read( r );
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
    char* read_start = rdb->mmapped_data + sizeof(MPD_RD_ID_T);
    
    /* count mmapped reads to check they match the saved mapped reads count */
    MPD_RD_ID_T num_indexed_reads = 0;

    /* Loop through all of the reads */
    while( ( (size_t)read_start - (size_t)rdb->mmapped_data ) 
           < rdb->mmapped_data_size )
    {
        /* start of a mapped_read_t */

        /* index the read_id_node's for this mapped_read_t */
        /* note - this assumes there is at least one read_id_node for each
         * mapped_read_t */
        while(true)
        {
            /* cast current location to read_id_node so it can be examined */
            read_id_node* node = (read_id_node*) read_start;

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
            rdb->index[num_indexed_reads].ptr = read_start;

            num_indexed_reads += 1;

            if( !(node->are_more) )
            {
                /* increment pointer to the start of the mapped read locations,
                 * and break */
                read_start += sizeof( read_id_node );
                break;
            }

            read_start += sizeof( read_id_node );
        }

        /* skip over the rest of the mapped_read_t in the mmapped memory */ 
        read_start = skip_mapped_read_locations_in_mapped_read_t( read_start );
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
        free_mapped_read( mapped_rd );
    }
    (*cond_prbs_db)->max_rd_id = max_rd_id;
    
    free_mapped_read( mapped_rd );
    
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
        free_mapped_read( mapped_rd );
    }

    free_mapped_read( mapped_rd ); // necessary?

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

