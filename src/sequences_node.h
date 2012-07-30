/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef SEQUENCES_NODE
#define SEQUENCES_NODE

#include <stdio.h>

#include "config.h"
#include "quality.h"
#include "mapped_location.h"

/*** Locations Node Type ********************************************************/
typedef void locations_node;

void
init_locations_node( locations_node** node );

size_t
size_of_locations_node( locations_node* node );

locations_node*
add_location_to_locations_node( locations_node* node,
                                GENOME_LOC_TYPE loc );

inline void
get_locations_from_locations_node( const locations_node* const node, 
                                   mapped_locations* results,
                                   const float penalty,
                                   const enum STRAND strnd );

/*** END Locations Node *********************************************************/

/* the maximum number of unique sequences that will
   be added to a seq node before it splits */
/* FIXME - add 0 length sequence node */
// typedef unsigned char NUM_SEQ_IN_SEQ_NODE_TYPE;
// #define MAX_SEQ_NODE_ENTRIES 254
typedef unsigned char NUM_SEQ_IN_SEQ_NODE_TYPE;
#define MAX_SEQ_NODE_ENTRIES 200

/* BUG - make this unsigned short ( see below for problem )*/
/* we use this to store the total number of bytes allocated for a sequence node */
typedef unsigned int MEMORY_SEGMENT_SIZE;
/* BUG - make this USHRT_MAX
   the reason why it is not is that for reasonable low numbners of repeats
   ( on the order of a hundred ) then this can be overrun, and, since we
   dont have any detection mechanism, it is a memory overflwo error. I changed
   this to a bigger size to avoid this problem, but the only real solution is
   to split this out into a dynamic node ( and, eventually, location node ) 
   if the size gets too large. That is a bit complicated so, I've increased the 
   size but plan to come back. Theoretically this could be overrun as well, but
   I've added in a few asserts and it should be good up to tens of thousands of
  repeats, so ti shouldnt be a problem in practice. 
*/
/* this value could be adjusted to be extra efficient depending on cache sizes */
#define MAX_SEQUENCES_NODE_SIZE UINT_MAX

/*** Sequence Node Specific Data Types *******************************************/

typedef struct {
    enum bool is_duplicate;
    enum bool is_pseudo;
    int location;
} insert_location;


typedef union __attribute__((__packed__)) {
    /* 
     * store an actual location in the genome. 
     */
    GENOME_LOC_TYPE loc;
    /* store the locations of the genome locations array */
    struct __attribute__((__packed__)) {
        /* 
         * store the index of the start of the array for this location, 
         * in units of sizeof(GENOME_LOC_TYPE) bytes after the end of 
         * the standard genome loc array.
         */
        signed locs_start :29;
        /* store the size of the array, in units of the number of entries */
        signed locs_size :19;
    } locs_array ;
} locs_union;

//#define MAX_LOC_ARRAY_START 24
//#define MAX_LOC_ARRAY_SIZE 65535
#define MAX_LOC_ARRAY_SIZE 262143 // 19 bit - signed ( 2**18 - 1 )
#define MAX_LOC_ARRAY_START 268435455 // 29 bit - signed bit ( 2^28 - 1 )


/* 
 * this is a weird data type. It's important to keep all of the sequences 
 * in contiguous blocks of memory both for when we need to spill to disk and to
 * save space. However, there is no good way to allocate these sorts of structures
 * in C without resorting to pointers, which take up too much room.
 * Therefore, we define the seq node as a type char and then define functions to 
 * return correctly casted pointers to the individual positions in memory. Of 
 * course allocation needs to be done seperately as well, but if the pointers
 * are always requested through the assisting functions, then the access should
 * be equivalent to the naive pointer implementation. To illustrate the 
 * naive implementation, in the comment below I define the structure 

typedef struct {
    * 
    * store the number of sequences in the node. Access this attribute by calling
    * num_sequence_types( sequecnes_node )
    *
    NUM_SEQ_IN_SEQ_NODE_TYPE num_sequence_types;

    *
    * store the length of allocated memory in bytes. We use a unsigned short, but 
    * DO NOT DEPEND ON THIS - I MAY TRY AND SAVE THE MEMORY IN THE FURTURE. Its
    * only going to be used to speed up reallocs, and in the hope of eventually 
    * eliminating malloc. 
    * TODO get this value from the heap, which will be very platform specific...
    *
    MEMORY_SEGMENT_SIZE used_size; 
    
    * TODO - definitely eliminate this maybe 
    MEMORY_SEGMENT_SIZE allocated_size;
    
    * store a bitmap of whether or not a genome_location is a pointer to 
    * a list of genome_location_types ( true genome locations ) or if it
    * is actually a genome location. This could be stored as a byte in 
    * genome location, but we use a bitmap global to sequences node to try 
    * and save a little bit of memory. Use the bitmap interface functions to 
    * modify this bitmap. The size ( in bytes ) is given by bitmap_size, the 
    * size in bits is given by num_sequence_types ( an attribute of this 
    * 'struct'. 
    BYTE* bitmap;

    *
    * store the sequences. This is a LETTER_TYPE array who's size is dependent on
    * the level that we are at in the tree. Ie, if we are at the very bottom of
    * the tree, then all of the sequence information is stored in the brnach nodes
    * so this array is of length 0. 
    * Access this array through
    * get_sequences_array_start( sequences_node*, LEVEL_TYPE )
    * seq_length to get the length of each individual sequence. 
    * Of course, the number of sequences ( not letters ) in the array is just 
    * given by num_sequence_types.
    *
    LETTER_TYPE* sequences;

    *
    * store the initial genome locations array. This is described in detail in
    * the genome location struct. 
    * Access this array by calling
    * get_genome_locations_array_start
    * The length is, again, just given by the num_sequence_types
    *
    GENOME_LOC_TYPE* locations;

    *
    * store the array of expanded genome locations. This is for cases in which, 
    * either the genome location can't fit in 31 bits, or there are multiple 
    * genome locations for a given sequence type. The offsets of each of these 
    * arrays and their lengths are stored in the genome_location. Access this
    * by calling 
    * get_overflow_genome_locations_array_start(sequences_node*, LEVEL_TYPE) 
    * the length is unknown, but should be determined by the other locations. 
    *
    GENOME_LOC_TYPE* overflow_locations;

)

*/

/* the sequences node declaration, as the compiler sees it */
typedef void sequences_node;

struct mapped_locations;

/**** BITMAP functions ********************************************************/

inline int 
check_bit( byte* bitmap, int index );

inline void 
set_bit( byte* bitmap, int index );

inline void
clear_bit( byte* bitmap, int index );

inline void 
print_bitmap( byte* bitmap, int size );

inline size_t
bitmap_size( int num_seqs );

/**** attribute access methods for sequence node ******************************/
/* 
 * We store a sequences node as a byte string, mostly to save memory but also
 * to make on disk indexes a bit easier to build. However, this means that the
 * normal struct->attribute will not work. Instead, I've implemented a set of
 * methods that return access all of the correct values.
 *
 */

inline NUM_SEQ_IN_SEQ_NODE_TYPE
get_num_sequence_types( const sequences_node* const seqs );

inline void
set_num_sequence_types( sequences_node* seqs, 
                        NUM_SEQ_IN_SEQ_NODE_TYPE value );

inline MEMORY_SEGMENT_SIZE
get_num_used_bytes( sequences_node* seqs );

inline void
set_num_used_bytes( sequences_node* seqs,
                    MEMORY_SEGMENT_SIZE value  );

inline MEMORY_SEGMENT_SIZE
get_num_allocated_bytes( sequences_node* seqs );

inline void
set_num_allocated_bytes( sequences_node* seqs,
                         MEMORY_SEGMENT_SIZE value  );

inline byte* 
get_bitmap_start( const sequences_node* const seqs );

inline sequences_node* 
insert_bit_into_bitmap( sequences_node* seqs, int insert_position );

/* return true if the genome_location at the given index is a ptr ( vs loc ) */
inline int
check_sequence_type_ptr( const sequences_node* const seqs, 
                         int index );

inline void
set_sequence_type_to_ptr( sequences_node* seqs, int index );

inline LETTER_TYPE* 
get_sequences_array_start( const sequences_node* const seqs, 
                           LEVEL_TYPE seq_length );

inline locs_union* 
get_genome_locations_array_start( const sequences_node* const seqs, 
                                  LEVEL_TYPE seq_num_letters );

inline GENOME_LOC_TYPE* 
get_overflow_genome_locations_array_start( 
    const sequences_node* const seqs, 
    LEVEL_TYPE seq_num_letters
);

/* return the number of bytes between the start of seqs and ptr */
inline MEMORY_SEGMENT_SIZE
bytes_before( sequences_node* seqs, void* ptr );

/* return the number of bytes between the start of ptr and the end of seqs */
inline MEMORY_SEGMENT_SIZE
bytes_after( sequences_node* seqs, void* ptr );

/* 
 * make empty space of size size at pointer. Ie, if the memory is 
 * 11111 ( in bytes ) then insert_empty_space( 2, 2  ) makes 1100111
 */
inline void
insert_memory ( sequences_node* seqs, 
                void* start,
                MEMORY_SEGMENT_SIZE size, 
                enum bool initialize_to_zero
    );

/**** SEQUENCES functions *****************************************************/

void init_sequences_node( sequences_node** seqs );

void free_seqs( sequences_node* seqs );

void
raw_print_sequences_node( sequences_node* seqs );

void print_sequences_node( sequences_node* seqs, int seq_length );

float
find_sequences_in_sequences_node(
        const sequences_node* const seqs,
        /* the curr penalty for seq */
        float curr_penalty,
        /* the maximum allowable penalty */
        float min_match_penalty,

        /* the seq of interest */
        const LETTER_TYPE* const seq,
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
        struct penalty_t* pa
    );

/* FWD declaration */
struct pseudo_locations_t;
sequences_node*
add_sequence_to_sequences_node(     
    struct pseudo_locations_t* ps_locs,
    sequences_node* seqs, 
    LETTER_TYPE* new_seq,
    LEVEL_TYPE num_letters,
    GENOME_LOC_TYPE location );

/**** END SEQUENCES functions *************************************************/

#endif // SEQUENCES_NODE
