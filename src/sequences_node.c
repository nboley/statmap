/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include <float.h>

#include "sequences_node.h"
#include "log.h"
#include "pseudo_location.h"
#include "quality.h"

/**** Locations Node **********************************************************/
/* 
 * Use this instead of a sequences node if level == 0.
 *
 */

void
init_locations_node( locations_node** node )
{
    *node = malloc( sizeof(unsigned short) );
    **((unsigned short**) node) = 0;
    return;
}

size_t
size_of_locations_node( locations_node* node )
{
    return sizeof( unsigned short ) 
        + *((unsigned short*) node)*sizeof( INDEX_LOC_TYPE );
}

locations_node*
add_location_to_locations_node( locations_node* node,
                                INDEX_LOC_TYPE loc )
{
    /* allocate memory for the new location */
    int num_locations = *( (unsigned short*) node );
    node = realloc( node, sizeof(unsigned short)
                    + (num_locations+1)*sizeof(INDEX_LOC_TYPE) );
    assert( node != NULL );

    /* add the new location */
    INDEX_LOC_TYPE* locs = 
        (INDEX_LOC_TYPE*) (((unsigned short*) node) + 1 );
    locs[ num_locations ] = loc;
    assert( loc.is_maternal == 0 && loc.is_paternal == 0 );
    
    /* increment the number of locations */
    *( (unsigned short*) node ) += 1;

    return node;
}

void
get_locations_from_locations_node( const locations_node* const node, 
                                   mapped_locations* results,
                                   const float penalty,
                                   const enum STRAND strnd,
                                   struct genome_data* genome )
{
    int num_locations = *( (unsigned short*) node );
    INDEX_LOC_TYPE* locs = 
        (INDEX_LOC_TYPE*) (((unsigned short*) node) + 1 );

    int i;
    for( i = 0; i < num_locations; i++ )
    {
        assert( (locs + i)->is_maternal == 0 && (locs + i)->is_paternal == 0 );
        add_and_expand_location_from_index(
                results,
                locs + i,
                strnd,
                penalty,
                genome
            );
    }
    
    return;

}


/**** BITMAP functions ********************************************************/

 int 
check_bit( byte* bitmap, int index )
{
    return ((bitmap[index/CHAR_BIT])&( 1 << (CHAR_BIT-1-(index%CHAR_BIT))));
}

 void
set_bit( byte* bitmap, int index )
{
    assert( !((bitmap[index/CHAR_BIT]) & ( 1 << (CHAR_BIT-1-index%CHAR_BIT) )) );
    bitmap[index/CHAR_BIT] |= ( 1 << (CHAR_BIT-1-index%CHAR_BIT));
}

 void
clear_bit( byte* bitmap, int index )
{
    assert( ((bitmap[index/CHAR_BIT]) & ( 1 << (CHAR_BIT-1-index%CHAR_BIT) )) );
    bitmap[index/CHAR_BIT] &= (~( 1 << (CHAR_BIT-1-index%CHAR_BIT)));
}

 void 
print_bitmap( byte* bitmap, int size )
{
    int i;
    for( i = 0; i < size; i++ )
    {
        if( check_bit(bitmap, i) ) {
            printf("1");
        } else {
            printf("0");
        }
    }
    printf("\n");
}

/* determine the size of the bitmap in bytes */
 size_t
bitmap_size( int num_seqs )
{
    if( num_seqs % CHAR_BIT == 0 )
        return num_seqs/CHAR_BIT;

    return num_seqs/CHAR_BIT + 1;
}

/**** END bitmap functions ************************************************/

 NUM_SEQ_IN_SEQ_NODE_TYPE
get_num_sequence_types( const sequences_node* const seqs )
{
    return *((NUM_SEQ_IN_SEQ_NODE_TYPE*) seqs);
}

 void
set_num_sequence_types( sequences_node* seqs, 
                        NUM_SEQ_IN_SEQ_NODE_TYPE value )
{
    assert( value <= MAX_SEQ_NODE_ENTRIES );
    *((NUM_SEQ_IN_SEQ_NODE_TYPE*) seqs) = value;
    return;
}

 MEMORY_SEGMENT_SIZE
get_num_used_bytes( sequences_node* seqs )
{
    return *( (MEMORY_SEGMENT_SIZE*) ( seqs 
        + sizeof( NUM_SEQ_IN_SEQ_NODE_TYPE )
    ));
}

 void
set_num_used_bytes( sequences_node* seqs,
                    MEMORY_SEGMENT_SIZE value  )
{
    assert( value < MAX_SEQUENCES_NODE_SIZE );
    *( (MEMORY_SEGMENT_SIZE*) ( seqs 
        + sizeof( NUM_SEQ_IN_SEQ_NODE_TYPE )
    )) = value;

    return;
}

 MEMORY_SEGMENT_SIZE
get_num_allocated_bytes( sequences_node* seqs )
{
    return *( (MEMORY_SEGMENT_SIZE*) ( seqs 
        + sizeof( NUM_SEQ_IN_SEQ_NODE_TYPE )
        + sizeof( MEMORY_SEGMENT_SIZE )
    ));
}

 void
set_num_allocated_bytes( sequences_node* seqs,
                         MEMORY_SEGMENT_SIZE value  )
{
    assert( value < MAX_SEQUENCES_NODE_SIZE );
    *( (MEMORY_SEGMENT_SIZE*) ( seqs 
        + sizeof( NUM_SEQ_IN_SEQ_NODE_TYPE )
        + sizeof( MEMORY_SEGMENT_SIZE )
    )) = value;

    return;
}

 byte* 
get_bitmap_start( const sequences_node* const seqs )
{
    return (byte*) ( seqs 
        + sizeof( NUM_SEQ_IN_SEQ_NODE_TYPE )
        + 2*sizeof( MEMORY_SEGMENT_SIZE )
    );
}

 sequences_node* 
insert_bit_into_bitmap( sequences_node* seqs, int insert_position )
{
   /* 
    * shift the bits down in the bytes below the bit of interest. 
    * ie, if the bitmap is |00100110|10100000| and low is 4, then we have
    * |0010I0110|10100000| -> |0010I011|01010000|. So, the algorithm is 
    * 0) add a new byte to the array if necessary
    * 1) allocate memory to store a copy of all bytes below the byte to be
    *    changed of course this is just all bytes low/8 through max_bytes
    * 2) copy shifted values into the new array
    * 3) test for the 1 bit in the old and, if it's set, set the high
    *    bit in the copy
    * 4) DEAL WITH THE BYTE THAT IS BEING INSERTED
    * 5) check to see if the one bit is set - if so, set the high bit in
    *    the copied array
    * 6) build the bit mask for all of the bits <= the bit we want to set 
    * 7) set the shifted, bitmasked bits
    * 8) set the correct bit ( actually, we calloc so it is pre-set )
    */

    /* get the number of used bits in the bitmap */
    NUM_SEQ_IN_SEQ_NODE_TYPE old_size
        = get_num_sequence_types( seqs );

    /* store the size of the array in bytes */    
    size_t num_bytes = bitmap_size( old_size );

    /* get a pointer to the bitmap */
    byte* bitmap = get_bitmap_start( seqs );

    /* 0 - add a new byte to the bitmap if necessary */
    int size = old_size + 1;
    /* if the extra bit pushed us intoa  new byte */
    if( size % CHAR_BIT == 1 ) {
        /* realloc new memory if necessary */
        size_t num_allcd_bytes = get_num_allocated_bytes( seqs );
        seqs = realloc( seqs, num_allcd_bytes + 1 );
        assert( seqs != NULL );
        /* increment the number of allcd, and used bytes */
        set_num_allocated_bytes( seqs, num_allcd_bytes + 1 );

        /* if seqs was moved during the realloc, we need to deal with that */
        bitmap = get_bitmap_start( seqs );

        /* insert memory into the sequences node at the correct location */
        insert_memory( seqs, bitmap + num_bytes, 1, true );

        /* initialize the new byte to 0 */
        /* this is done automatically by insert memory */
        
        /* update the bytes use counters */
        set_num_used_bytes( seqs, get_num_used_bytes( seqs ) + 1 );
        /* the size of the bitmap, in bytes */
        num_bytes++;       
    }

    /* add a logic branch for the case that the bit is being added to the end */
    /* 
     * ie, if the insert position == size, then we really want to append, and
     * since we initialize to zero, this is done for us. 
     */
    if( insert_position == (size-1) ) 
    {
        return seqs;
    }
    
    /* 1 */
    /* store a copy of all the bits below the correct level */
    unsigned int target_byte_index = insert_position/CHAR_BIT;
    /* allocate memory for the copy */
    /* TODO - eliminate the full allocation */
    byte* bitmap_copy = calloc( num_bytes, sizeof(byte) );

    /* 2&3 - make the copy by looping through the bytes */
    unsigned int i;
    for( i = target_byte_index + 1; i < num_bytes; i++ )
    {
        /* 2 - copy over the shifted bits */
        bitmap_copy[i] = ( bitmap[i] >> 1 );
        /* 3 - if the one bit is set */
        if( 1 == ((bitmap[i-1])&1) ) {
            /* set the high bit */
            /* noting that the high bit *cant* be 1 due to the bit shift */
            /* TODO - use a static bitmask for this */
            bitmap_copy[i] += ( 1 << (CHAR_BIT - 1) );
        }
    }

    /* copy the shifted bits back into the original array */
    for( i = target_byte_index + 1; i < num_bytes; i++ )
    {
        bitmap[i] = bitmap_copy[i];
    }
    
    /* free the copy */
    free( bitmap_copy );

    /* 6 - build the bitmask for all of the bits >= the bit we want to set */
    unsigned int bitmask_size = CHAR_BIT - insert_position%CHAR_BIT - 1;
    byte bitmask = 0;
    for( i = 0; i < bitmask_size+1; i++ )
    {
        bitmask = ((bitmask << 1) + 1);
    }
    
    /* 7 - set the bitshifted bits */
    /* store the low bytes, and then shift them right */
    byte low_bits = ((bitmap[target_byte_index])&bitmask) >> 1;
    /* mask out the low bits in the real data */
    bitmap[target_byte_index] &= (~bitmask);
    /* combine the two */
    bitmap[target_byte_index] |= low_bits;
    
    /* return the potentially realloced seqs */
    return seqs;
}

 int
check_sequence_type_ptr( const sequences_node* const seqs, 
                         const int index )
{
   /* 
    * note that it is the callers responsibility to ensure that
    * index is sane, ie within the limits of what we would expect 
    * given the number of sequence types in seqs 
    *
    */
    
   /* first get the bitmap */
    byte* bitmap = get_bitmap_start( seqs );

    /* return true if it is set, and thus a ptr */
    return check_bit( bitmap, index );
}

 void
set_sequence_type_to_ptr( sequences_node* seqs, int index )
{
   /* 
    * note that it is the callers responsibility to ensure that
    * index is sane, ie within the limits of what we would expect 
    * given the number of sequence types in seqs 
    *
    */
    
   /* first get the bitmap */
    byte* bitmap = get_bitmap_start( seqs );

    /* set the correct bit */
    set_bit( bitmap, index );
    
    return;
}


 LETTER_TYPE* 
get_sequences_array_start( const sequences_node* 
                           const seqs, 
                           const LEVEL_TYPE seq_num_letters )
{
    /* 
     * IT IS THE CALLERS RESPONSIBILITY TO ensure that seq_length
     * is greater than 0. If it is not, the 
     */
    assert( seq_num_letters > 0 );
    
    /*
     * if this was the pseudo struct, this would just be
     * return seqs->sequences;
     */    
    LETTER_TYPE* loc = (LETTER_TYPE*) (
        /* the start, since it is a byte we can do size of */
        seqs
        /* 
         * deal with the memory allocated for the number of sequences in the
         * node type 
         */
        + sizeof(NUM_SEQ_IN_SEQ_NODE_TYPE)
        /* memory used */
        + sizeof(MEMORY_SEGMENT_SIZE)
        /* memory allocated */
        + sizeof(MEMORY_SEGMENT_SIZE)
        /* the sizeof the bitmap */
        + bitmap_size( get_num_sequence_types(seqs) )
     );

    return loc;
}

 locs_union* 
get_genome_locations_array_start( const sequences_node* const seqs, 
                                  const LEVEL_TYPE seq_num_letters )
{
    /*
     * if this was the pseudo struct, this would just be
     * return seqs->locations;
     */
    return (locs_union*) ( 
        /* start of the sequences array chunk */
        (byte*) get_sequences_array_start( seqs, seq_num_letters )
        /* the length of sequences array */
        + seq_num_letters*sizeof(LETTER_TYPE)*get_num_sequence_types(seqs)
    );
}

 INDEX_LOC_TYPE* 
get_overflow_genome_locations_array_start( 
    const sequences_node* const seqs, LEVEL_TYPE seq_num_letters
)
{
    return (INDEX_LOC_TYPE*) (
        /* start of the genome_locations array */
        (byte*) get_genome_locations_array_start( seqs, seq_num_letters )
        /* length of genome locations array */
        + get_num_sequence_types(seqs)*sizeof(locs_union)
    );
}

/* return the number of bytes between the start of seqs and ptr */
 MEMORY_SEGMENT_SIZE
bytes_before( sequences_node* seqs, void* ptr )
{
    return (MEMORY_SEGMENT_SIZE) ( ptr - (void*)seqs );
}

/* return the number of bytes between the start of ptr and the end of seqs */
 MEMORY_SEGMENT_SIZE
bytes_after( sequences_node* seqs, void* ptr )
{
    return (MEMORY_SEGMENT_SIZE) ( 
        get_num_used_bytes( seqs ) - 
        bytes_before( seqs, ptr )
    );
}

/* 
 * make emtpy space of size size at pointer. Ie, if the memory is 
 * 11111 ( in bytes ) then insert_empty_space( 2, 2  ) makes 1100111
 */
 void
insert_memory ( sequences_node* seqs, 
                void* start,
                MEMORY_SEGMENT_SIZE size, 
                enum bool initialize_to_zero
)
{
    /* 
     * we assume that the realloc has already been done for us, and that the 
     * extra allocated space is at the end of the memory block.
     */

    /* find the size of the memory that needs to be moved after this */
    size_t bytes_to_move = bytes_after( seqs, start  );

    /* move the bytes after the insert location */
    memmove( start+size, start, bytes_to_move  );

    /* TODO - make this a memset */
    /* initialize the new memory to zero */
    if( initialize_to_zero )
    {
        size_t i;
        for( i = 0; i < size; i++ )
        {
            ((byte*)start)[i] = 0;
        }
        return;
    }
}

void
init_sequences_node( sequences_node** seqs )
{
    /* 
     * We initialize the struct to have 0 genome locations, so we only 
     * need to allocate space to store the number of genome locations.
     * When we add a genome location, then we will grow the allocated size.
     */
    size_t alloc_size = 
        sizeof(NUM_SEQ_IN_SEQ_NODE_TYPE) 
        + 2*sizeof(MEMORY_SEGMENT_SIZE);

    assert( alloc_size < MAX_SEQUENCES_NODE_SIZE ); 

    /* note that the calloc initializes num_seq_types and the bitmap */
    *seqs = calloc( alloc_size, 1 );
    set_num_allocated_bytes( *seqs, alloc_size );
    set_num_used_bytes( *seqs, alloc_size );
}

void
free_seqs( sequences_node* seqs )
{
    free( seqs );
}

/*
 * Return the insert location for a new sequence fragment in units
 * of seq len. So, if the array consists of 2 sequences each with 
 * 2 letters, then the seqs_array is of length 4 LETTER_TYPE, and 
 * inserting into the middle position would correspond to an insert
 * location of 1 ( not 2 ).
 */
 insert_location
find_insert_location( sequences_node* seqs, 
                      LETTER_TYPE* new_seq, 
                      LEVEL_TYPE seq_len )
{
    /* find the start of the sequences array */
    LETTER_TYPE* seqs_array = get_sequences_array_start( seqs,  seq_len );

    /* determine if this sequence type already exists */
    /* 
     * we use a binary search to find if it matches, and the insert location
     * if it does not match 
     */  
    int low = 0;
    int high = get_num_sequence_types( seqs );
    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        int comparison = cmp_words( seqs_array + mid*seq_len, 
                                    new_seq, seq_len );
        if( comparison < 0 ) {
            low = mid + 1; 
        } else {
            /* can't be high = mid-1: here A[mid] >= value, */
            /* so high can't be < mid if A[mid] == value    */
            high = mid; 
        }
    }

    /* make sure the binary search is working */
    assert( low <= high );
    assert( low <= get_num_sequence_types(seqs) );
    assert( low >= 0 );

    insert_location rv;
    rv.location = low;    
    
    /* check to see if there is a match at this point */
    if( 0 < get_num_sequence_types( seqs )
        && low < get_num_sequence_types(seqs)
        && 0 == cmp_words( seqs_array+low*seq_len, new_seq, seq_len ) )
    { 
        rv.is_duplicate = true;
        
        if( !check_sequence_type_ptr( seqs, low )
            && ( get_genome_locations_array_start( seqs, seq_len )
                 + rv.location )->loc.chr == PSEUDO_LOC_CHR_INDEX )
        { 
            rv.is_pseudo = true; 
        }
        else { rv.is_pseudo = false; };
    } else {
        rv.is_duplicate = false;
        rv.is_pseudo = false;
    } 

    return rv;
}

/*
 * Add a sequence that doesnt currently exist into a sequences node. This means
 * that we dont need to deal with the extended GENOME_LOC array
 */
static  sequences_node*
add_new_sequence_to_sequences_node( sequences_node* seqs, 
                                    /* where in the seqs array this belongs */
                                    int insert_loc,
                                    /* the seq to add */
                                    LETTER_TYPE* new_seq,
                                    /* the length of seq in this node */
                                    int seq_len,
                                    /* the locations of this sequence */
                                    INDEX_LOC_TYPE loc  )
{
    /* allocate the space - we need room for 1 sequence and 1 genome_loc */
    size_t curr_size = get_num_used_bytes( seqs );
    size_t size_needed = curr_size 
        + seq_len*sizeof(LETTER_TYPE)
        + sizeof( locs_union );

    /* if we need to, allocate more memory */
    if( size_needed > get_num_allocated_bytes( seqs ) )
    {
        sequences_node* new_seqs 
            = realloc(seqs, size_needed);
        
        assert( new_seqs != NULL );
        set_num_allocated_bytes( new_seqs, size_needed );
        /* initialize the new bytes */
        memset( new_seqs + curr_size, 0, size_needed - curr_size );
        
        seqs = new_seqs;
    }

    /* make space in the genome locations array */
    /* get the start of the genome locations array */
    locs_union* loc_start = 
        get_genome_locations_array_start( seqs, seq_len );
    /* move the pointer to the new location start */
    loc_start += insert_loc;
    
    /* make space for the new genome location, and initialize it */
    insert_memory( seqs, loc_start, sizeof(locs_union), true );

    /* update the number of used bytes */
    curr_size += sizeof( locs_union );
    set_num_used_bytes( seqs, curr_size );
    
    /* this is done explicitly when the bitmap is added to */
    assert( loc.loc <= LOCATION_MAX && loc.loc >= 0 );
    loc_start->loc = loc;
    
    /* make space for the sequence, and initialize it */
    LETTER_TYPE* seqs_start =
        get_sequences_array_start( seqs, seq_len );
    seqs_start += (seq_len*insert_loc);

    insert_memory( seqs, seqs_start, seq_len*sizeof(LETTER_TYPE), true );

    /* update the number of used bytes */
    curr_size += seq_len*sizeof(LETTER_TYPE);
    set_num_used_bytes( seqs, curr_size );
    assert( get_num_used_bytes(seqs) == get_num_allocated_bytes(seqs) );
    
    memcpy( seqs_start, new_seq, sizeof(LETTER_TYPE)*seq_len );

    /* insert the bit into the bitmap */
    seqs = insert_bit_into_bitmap( seqs, insert_loc );
    
    /* increment the number of sequence types */
    NUM_SEQ_IN_SEQ_NODE_TYPE curr_num_seq_types = 
        get_num_sequence_types(seqs);
    set_num_sequence_types( seqs, curr_num_seq_types+1 );
    
    return seqs;
}

/*
 * Add a sequence that has already been added into a sequences node. This means
 * that we dont need to add a new sequence, but we do need to either, convert
 * the genome location into a pointer and create and extended genome location
 * list, or we need to just update the pointer, add an entry to the extended list
 * and update the previous entries.
 */
static  sequences_node*
add_duplicate_sequence_to_sequences_node(
    /* this is necessary when we need to add pseudo locations */
    struct genome_data* genome,
    struct pseudo_locations_t* ps_locs,
    sequences_node* seqs, 
    /* where in the seqs array this belongs, in seq_len units */
    int insert_loc,
    /* the length of seq in this node */
    LEVEL_TYPE seq_len,
    /* the locations of this sequence */
    INDEX_LOC_TYPE new_loc  )
{
    // DEBUG
    /*
    printf("adding duplicat: (add_duplicate_sequence_to_sequences_node)\n");
    printf("Chr: %i\t Loc: %i\n", new_loc.chr, new_loc.loc );
    */
    
    /* move the pointer to the new genome location start */
    locs_union* loc = 
        get_genome_locations_array_start( seqs, seq_len ) + insert_loc;
    
    /* if there is only 1 other copy of this ( at least so far ) */
    /* 
     * ie, the sequence at insert loc is *not* a pointer ( and we already
     * know that this has to be a duplicate given the function call 
     */
    if( !check_sequence_type_ptr(seqs, insert_loc) )
    {
        /* we add 2 genome locations 1 for the old, and 1 for new gen location */
        size_t num_used_bytes = get_num_used_bytes(seqs);
        size_t new_size = num_used_bytes + 2*sizeof(INDEX_LOC_TYPE);
        
        /* if necessary, realloc */
        if( new_size > get_num_allocated_bytes(seqs) )
        {
            /* allocate space for the extended genome locations array */
            seqs = realloc( seqs, new_size );
            assert( seqs != NULL );
            /* set the new allocated space size */
            set_num_allocated_bytes( seqs, new_size );
            /* initialize the new memory */
            memset( seqs + num_used_bytes, 0, new_size - num_used_bytes);
            /* update the pointer to the genomic locations array, if necessary */
            loc = get_genome_locations_array_start( seqs, seq_len ) + insert_loc;
        }

        /* 
         * the start of the new overflow locations - or, the end 
         * of the previous list. We dont care about the ordering of 
         * this list ( why would we ) so we simply append to the 
         * end. 
         *
         * Also, note that we havn't updated the used size yet so we
         * can use the get_used_bytes interface to find the end.
         */
        INDEX_LOC_TYPE* of_locs = 
            (INDEX_LOC_TYPE*) ( seqs + num_used_bytes );
        
        /* add the entries to the list */
        of_locs[0] = loc->loc;
        of_locs[1] = new_loc;

        /* now that the space is being used, set the new used size of seqs */
        set_num_used_bytes( seqs, new_size );

        assert( get_num_used_bytes(seqs) == get_num_allocated_bytes(seqs) );

        /** convert the genome_location_node to a pointer type */
        /* make the list a pointer */
        set_sequence_type_to_ptr( seqs, insert_loc );

        /* 
         * set the new locations array start to be the number of bytes after the
         * full array divided by the size of INDEX_LOC_TYPE. This is really
         * just the index of the first entry in the array.
         */
        size_t egla_len = bytes_after(
            seqs, get_overflow_genome_locations_array_start( seqs, seq_len )
        );

        loc->locs_array.locs_size = 2;

        assert( egla_len % sizeof( INDEX_LOC_TYPE ) == 0 ); 
        
        loc->locs_array.locs_start = 
            ( egla_len/sizeof(INDEX_LOC_TYPE) ) - 2;
            
    } 
    /* this has already been converted to a pointer */
    else {
        assert( loc->locs_array.locs_size >= 0 ); 
        assert( loc->locs_array.locs_start >= 0 ); 

        /* 
         * If this sequence has enough locations to warrant having
         * a pseudo location, then split it out
         */
        if( PSEUDO_LOC_MIN_SIZE - 1 == loc->locs_array.locs_size )
        {
            /*** add a pseudo location to the pseudo location DB */
            unsigned int psloc_index = 
                add_new_pseudo_location( ps_locs  );

            /* find the start */
            INDEX_LOC_TYPE* of_locs = 
                get_overflow_genome_locations_array_start( seqs, seq_len )
                + loc->locs_array.locs_start;
            
            int i;
            for( i = 0; i < loc->locs_array.locs_size; i++ )
            {
                add_new_loc_to_pseudo_location( 
                    ps_locs->locs + psloc_index, of_locs + i, genome );
                // printf("%i: Chr %i\tLoc: %u\n", i, of_locs[i].chr, of_locs[i].loc);
            }

            /* add the new location */
            add_new_loc_to_pseudo_location( 
                ps_locs->locs + psloc_index, &new_loc, genome );
            
            /* remove the array memory */
            /* 1) shift the other memory forward */
            size_t bytes_to_move = bytes_after( 
                seqs, of_locs + PSEUDO_LOC_MIN_SIZE - 1 ); 
            memmove( of_locs, of_locs + PSEUDO_LOC_MIN_SIZE - 1, bytes_to_move );
            /* 2) realloc to shrink used memory */
            size_t new_size = 
                get_num_used_bytes(seqs) 
                - (PSEUDO_LOC_MIN_SIZE-1)*sizeof(INDEX_LOC_TYPE);
            seqs = realloc( seqs, new_size );
            set_num_allocated_bytes( seqs, new_size );
            set_num_used_bytes( seqs, new_size );
            loc = get_genome_locations_array_start( seqs, seq_len ) + insert_loc;
            
            /* 
             * Because we just removed memory in the middle of this array, 
             * we need to move the start pointers forward for every position 
             * after this. This is a bit confusing - the array of gen locs
             * is ordered by sequence, not by the extended array. So, we need
             * to look at every genome location and if it is a pointer and
             * it starts after the array we just modified, increment the start 
             * location.
             */
            locs_union* locs_start = 
                get_genome_locations_array_start( seqs, seq_len );
            
            int new_start = loc->locs_array.locs_start;
            for( i = 0; i < get_num_sequence_types(seqs); i++ )
            {
                if( check_sequence_type_ptr( seqs, i )
                    && locs_start[i].locs_array.locs_start > new_start  )
                {
                    /* Make sure that the data still fits */
                    assert( locs_start[i].locs_array.locs_start 
                            < MAX_LOC_ARRAY_START );
                    locs_start[i].locs_array.locs_start 
                        -= ( PSEUDO_LOC_MIN_SIZE - 1 );
                }
            }
            
            /* unset the array bit */
            clear_bit( get_bitmap_start( seqs ), insert_loc );
            
            /* convert the old location to a proper pseudo locations. This means
               giving it the new pseudoloc_chr, the new pseudo position, and resetting
               the diploid mapping flags. */
            assert( psloc_index <= LOCATION_MAX );            
            memset( &(loc->loc), 0, sizeof(INDEX_LOC_TYPE) );
            loc->loc.chr = PSEUDO_LOC_CHR_INDEX;
            loc->loc.loc = psloc_index;
        } else {
            /* we add memory for the new genome location */
            size_t new_size = get_num_used_bytes(seqs) + sizeof(INDEX_LOC_TYPE);
            
            if( new_size > get_num_allocated_bytes(seqs) )
            {
                /* allocate space for the extended genome locations array */
                seqs = realloc( seqs, new_size );
                assert( seqs != NULL );
                /* set the new allocated space size */
                set_num_allocated_bytes( seqs, new_size );
                /* update the pointer to the genomic locations array, if necessary */
                loc = get_genome_locations_array_start( seqs, seq_len ) + insert_loc;
            }
            
            /* the start of the new overflow locations - note that we havnt updated
               the used size yet. */
            INDEX_LOC_TYPE* of_locs = 
                get_overflow_genome_locations_array_start( seqs, seq_len )
                + loc->locs_array.locs_start;
            
            /* move memory to make space for the new location */
            insert_memory( seqs, of_locs, sizeof(INDEX_LOC_TYPE), false );
            
            /* set the new value */
            of_locs[0] = new_loc;
            
            /* now that the space is being used, set the new used size of seqs */
            set_num_used_bytes( seqs, new_size );
            
            assert( get_num_used_bytes(seqs) == get_num_allocated_bytes(seqs) );
            
            /* Make sure that we havn't grown beyond what the array can hold */
            // assert( loc->locs_array.locs_size < MAX_LOC_ARRAY_SIZE );
            
            /* update the number of sequence */
            assert( loc->locs_array.locs_size < MAX_LOC_ARRAY_SIZE );
            loc->locs_array.locs_size += 1;        
            
            /* 
             * Because we just inserted memory in the middle of this array, 
             * we need to move the start pointers forward for every position 
             * after this. This is a bit confusing - the array of gen locs
             * is ordered by sequence, not by the extended array. So, we need
             * to look at every genome location and if it is a pointer and
             * it starts after the array we just modified, increment the start 
             * location.
             */
            locs_union* locs_start = 
                get_genome_locations_array_start( seqs, seq_len );
            
            /* TODO - move this declaration earlier */
            int new_start = loc->locs_array.locs_start;
            int i;
            for( i = 0; i < get_num_sequence_types(seqs); i++ )
            {
                if( check_sequence_type_ptr( seqs, i )
                    && locs_start[i].locs_array.locs_start > new_start  )
                {
                    /* Make sure that the data still fits */
                    assert( locs_start[i].locs_array.locs_start 
                            < MAX_LOC_ARRAY_START );
                    locs_start[i].locs_array.locs_start++;
                }
            }
        }
    }     
   return seqs;
}

sequences_node*
add_sequence_to_sequences_node(     
    /* this is necessary when we need to add pseudo locations */
    struct genome_data* genome,
    struct pseudo_locations_t* ps_locs,
    sequences_node* seqs, 
    LETTER_TYPE* new_seq,
    LEVEL_TYPE num_letters,
    INDEX_LOC_TYPE loc )
{
    /* find the location in the seqs node that the new sequence should go */
    insert_location il = 
        find_insert_location( seqs, new_seq, num_letters );
    
    /* 
     * if this is a pseudo location, then all we need to do is add
     * this real location to the pseudo locations DB 
     */
    if( il.is_pseudo == true )
    {
        int ps_loc_index = 
            get_genome_locations_array_start( seqs, num_letters )
            [il.location].loc.loc;
        add_new_loc_to_pseudo_location(
                ps_locs->locs+ps_loc_index, &loc, genome );
        return seqs;
    } 
    else if( il.is_duplicate == true )
    {
        /* debugging */
        /* Make sure that the sequence actually is a duplicate */
        /* That is, 
           1) Make sure that the sequence at il is actually a duplicate. Test
              via a simple mem compare 
        */
        assert( cmp_words( get_sequences_array_start( seqs,  num_letters )
                           + il.location*num_letters, 
                           new_seq, 
                           num_letters ) == 0 
        );
        
        return add_duplicate_sequence_to_sequences_node( 
            genome, ps_locs, seqs, il.location, num_letters, loc 
        );
    } else {
        return add_new_sequence_to_sequences_node( 
            seqs, il.location, new_seq, num_letters, loc 
        );
    }
    
    assert( false );
    return NULL;
}


float
find_sequences_in_sequences_node(
        const sequences_node* const seqs,
        /* the curr penalty for seq */
        float curr_penalty,
        /* the maximum allowable penalty */
        float min_match_penalty,

        /* the length of a full sequence */
        const int seq_length,
        /* the total num of letters in a seq */
        const int num_letters,
        /* the current level in the tree */
        const int node_level,
        /* the strand of the search seq */
        const enum STRAND strnd,

        mapped_locations* results,

        /* the penalty array */
        struct penalty_t* pa,

        struct genome_data* genome
    )
{

    /* store the updated maximum penalty, for if we find matches */
    float max_penalty = -FLT_MAX;

    /* we assume that the results are initialized */
    assert( results != NULL );
    assert( results->length <= results->allocated_length );

    if( results->allocated_length == USHRT_MAX ) {
        statmap_log( LOG_FATAL, "Overran max results length in find_sequences_in_sequences_node()" );
    }

    /* get the total number of sequences that we need to consider */
    NUM_SEQ_IN_SEQ_NODE_TYPE num_seqs = get_num_sequence_types( seqs );

    /* get the start of the sequences array */
    LETTER_TYPE* seq_array_start = get_sequences_array_start( 
        seqs, num_letters - node_level );

    int i;
    for( i = 0; i < num_seqs; i++ )
    {
        /* 
         * store the cumulative penalty - the total penalty up till now
         * for the i'th sequence 
         */
        float cum_penalty = multiple_letter_penalty(
            /* the start of the current seq in the array */
            seq_array_start + i*(num_letters-node_level), 
            node_level, 
            seq_length, 
            num_letters,
            min_match_penalty - curr_penalty,
            pa
        );
        
        /* // debugging code
           printf("%i: Ref %i%i\tNode %i%i\tPenalty: %e\n", 
           i, seqs->sequences[2*i], seqs->sequences[2*i + 1], 
           seq[0], seq[1], cum_penalty);
        */

        if( cum_penalty <= 0.5 )
        {
            /* 
             * if we havnt skipped to the next sequence, then the penalty should
             * still be above the minimum - so add it to the results
             */
            assert( cum_penalty >= min_match_penalty );
            
            /* update the cumulative penalty */
            cum_penalty += curr_penalty;

            /* if appropriate, update the max penalty */
            if( cum_penalty > max_penalty )
                max_penalty = cum_penalty;

            /* get the correct locaton */
            locs_union loc = 
                get_genome_locations_array_start( 
                    seqs, num_letters - node_level )[i];
            
            
            /* 
             * if the correct bit is set, then the genome location is a pointer.
             */        

            if( !check_sequence_type_ptr(seqs, i) )
            {
                add_and_expand_location_from_index(
                        results,
                        &(loc.loc),
                        strnd,
                        cum_penalty,
                        genome
                    );
            } else 
            /* Add all of the locs in the referenced array */
            {
                INDEX_LOC_TYPE* locs 
                    = get_overflow_genome_locations_array_start( 
                        seqs,  num_letters - node_level )
                    + loc.locs_array.locs_start;

                /*
                 * DEBUG the below memory error 
                 * TODO - remove this comment.
                 *
                for(j = 0; j < loc.locs_array.locs_size; j++ )
                {
                    printf( "Loc Ptr: %p\t Loc: %ui|%ui|%ui|%ui|%ui|%ui\n", 
                            locs + j,
                            ( (char*) (locs + j ))[0],
                            ( (char*) (locs + j ))[1],
                            ( (char*) (locs + j ))[2],
                            ( (char*) (locs + j ))[3],
                            ( (char*) (locs + j ))[4],
                            ( (char*) (locs + j ))[5]
                        );
                }
                */
                
                int j;                
                for(j = 0; j < loc.locs_array.locs_size; j++ )
                {   
                    /* 
                     * This is a very messy way of avoiding a valgrind
                     * error that was, potentially, a real problem. 
                     * Possibly because valgrind is stupid, and possibly
                     * because gcc is bad at alignment, if I try and get 
                     * i_loc by locs[j], then it attempts to read 2 full words,
                     * overrunning the boundary of the array. There shouldnt 
                     * be any problem with this but, to make sure, I declare a 
                     * temporary variable here, store it inside, and then 
                     * call the function. 
                     */
                    /*
                    INDEX_LOC_TYPE i_loc = *( 
                        (INDEX_LOC_TYPE*)
                        (((char*)locs) + j*sizeof(INDEX_LOC_TYPE))
                    );
                    */
                    
                    INDEX_LOC_TYPE tmp_loc = locs[j];
                    add_and_expand_location_from_index(
                            results,
                            &tmp_loc,
                            strnd,
                            cum_penalty,
                            genome
                        );
                }
            }
        }
        /* we are done with this sequence - move onto the next */
    }

    return max_penalty;
}

void
raw_print_sequences_node( sequences_node* seqs )
{
    size_t i;
    MEMORY_SEGMENT_SIZE num_used_bytes = get_num_used_bytes(seqs);
    for( i = 0; i < num_used_bytes; i++ )
    {
        printf("|%i", ((unsigned char*) seqs)[i]);
    };
    printf("|\n");

}

void
print_sequences_node( sequences_node* seqs, int seq_length )
{
    printf("Number Distinct Sequences: %i\n", get_num_sequence_types(seqs));
    printf("Number Used Bytes: %i\n", get_num_used_bytes(seqs));
    printf("Number Allocated Bytes: %i\n", get_num_allocated_bytes(seqs));
    byte* bitmap = get_bitmap_start( seqs );
    printf( "Bitmap: " );
    print_bitmap( bitmap, get_num_sequence_types(seqs) );
    
    if( seq_length > 0 )
    {
        /* print out the stored sequences and their locations */
        int i, j;

        /* store the beginning of the current sequence that we are on */
        LETTER_TYPE* curr_seq = get_sequences_array_start( seqs, seq_length );

        /* get the genome locations for this sequence */
        locs_union* locs = 
            get_genome_locations_array_start( seqs, seq_length );

        printf("Stored Sequences: \n");
        for( i = 0; i < get_num_sequence_types(seqs); i++ )
        {
            /* print out the current sequence */
            printf("%i\t", i+1);
            for( j = 0; j < seq_length; j++ )
            {
                printf("%i", curr_seq[j]);
            }
            printf("\n");

            /* print out the genome location */
            if( check_sequence_type_ptr( seqs, i ) )
            {
                printf("\t Is Pntr\tStart: %i\tSize: %i\n", 
                       locs[i].locs_array.locs_start,
                       locs[i].locs_array.locs_size
                );

                INDEX_LOC_TYPE* e_locs = 
                    get_overflow_genome_locations_array_start( seqs, seq_length )
                    + locs[i].locs_array.locs_start;

                for( j = 0; j < locs[i].locs_array.locs_size; j++ )
                {
                    printf("\t\t%i\n", e_locs[j].loc);
                }
            } else {
                if( locs[i].loc.chr == PSEUDO_LOC_CHR_INDEX ) 
                {
                    printf("\t Is PSEUDO Loc\tLoc: %i\n", locs[i].loc.loc);
                } else {
                    printf("\t Is Loc\tChr: %i\tLoc: %i\n", locs[i].loc.chr, locs[i].loc.loc);
                }
            }

            /* move tot he next sequence in the array */
            curr_seq += ( seq_length*sizeof(LETTER_TYPE) );
        }
    }
}


