/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef INDEX_GENOME
#define INDEX_GENOME

#define DONT_PROFILE_MEMORY_USAGE

#include <stdio.h>

#include "config.h"
#include "statmap.h"
#include "genome.h"
#include "sequences_node.h"
#include "mapped_location.h"
#include "error_correction.h"

struct indexable_subtemplates*
build_indexable_subtemplates_from_read_subtemplate(
        struct read_subtemplate* rst,
        struct index_t* index,
        enum bool use_random_subtemplate_offset
    );

int
search_index_for_read_subtemplate(
        struct read_subtemplate* rst,
        struct mapping_params* mapping_params,
        
        mapped_locations*** search_results,

        struct genome_data* genome,
        enum bool use_random_subtemplate_offset);

/* a hash to store chromosome name */
char* chr_names[CHR_NUM_MAX+1];

/* NODE TYPES
*
* store the node type of the children nodes
* this is the way in which the children nodes are stored. The 3 possible
* values are
* d: dynamic branches
*    children is an array of structs that contain basepair type and a 
*    pointer to the child node. This is for cases in which there are not so 
*    enough children to justify using a static array, but there are enough
*    sequences below this level that we don't want to save them one by one.
*    we know we've reached the end when the child pointer is NULL.
* s: static branches
*    children is an array of pointers whose indexes are equal to the value
*    of the letter that they point two. We use this node type when the 
*    number of nodes in the dynamic branches gets big enough that the 
*    not having to store the letter saves memory.
* q: sequences
*    store the actual sequences. This is *always* the node type when we
*    get to the bottom of the tree, but it may appear higher up. Here, 
*    children is pointer to a sequences_node struct
*/

typedef struct
{
    /* store the types of the nodes that children points to */ 
    NODE_TYPE type;
    /* pointer to the children array - the form is determined by node_type */
    void* node_ref;
} static_node_child;

/* static node in the tree */
/* 
 * since static nodes lengths are determined by the letter size,
 * we dont need anything except a pointer to the beggining of the
 * array.
 */

typedef static_node_child static_node;

typedef struct {
    /* store the type of the node referenced by child */
    NODE_TYPE type;
    /* the child's letter */
    /*
     * Dynamic nodes work like the branch node of any search tree - you
     * compare the letter at position level_type to this letter and then
     * follow the pointer if necessary. child_type stores the letter at this
     * level for all children.
     */
    LETTER_TYPE letter;
    /* pointer to the child - the form is determined by node_type */
    void* node_ref;
} dynamic_node_child;

/* 
 * Dynamic node pseudo struct
 *
 * In order to save the 4 byte overhead for a dynamic node and to 
 * make it easier to marshall these to be saved on disk, we keep the 
 * dynamic node in a packed structure that is an unsigned short 
 * for the number of children followed by an array of dynamic_node_child
 * of length num_children. Use get_dnode_num_children( dynamic_node* dnode ) and
 * get_dnode_children( dynamic_node* dnode ) to get the two, rather than 
 * accessing the structure elements directly. 
 *
 * dynamic node in the tree *
 *

typedef struct {
    // number of children
    unsigned short children_length;
    dynamic_node_child* children;
} dynamic_node;

 *
 */

typedef char dynamic_node;

unsigned short
get_dnode_num_children( const dynamic_node* const dnode );

void
set_dnode_num_children( dynamic_node* const dnode, 
                        const unsigned short value );

dynamic_node_child*
get_dnode_children( const dynamic_node* const dnode );

size_t
size_of_dnode( const dynamic_node* const node );

size_t 
size_of_snode( );

/*
 * An entry in an execution queue. Each node is a potential match for a 
 * sequence ( or group of sequences. For instance, if a potential match
 * has a sequence node, ( node_type = 'q' ) at level 4 with a PENALTY of
 * -0.8, then we can call find match with this info and a sequence to output 
 * potential matches.
 */

#define STACK_GROWTH_FACTOR 100
#define STACK_INITIAL_FIRST_INDEX 50

typedef struct {
    void* node;
    NODE_TYPE node_type;
    int node_level;
    float penalty;
    float scaled_penalty;
    enum STRAND strnd;
} potential_match;

typedef struct {
    size_t first_index; /* inclusive */
    size_t last_index;  /* exclusive */
    /* so if first = last index, then the stack is empty */
    size_t allocated_size;
    potential_match* matches;
} potential_match_stack;

void
init_pmatch_stack( potential_match_stack** stack );

void
free_pmatch_stack( potential_match_stack* stack );

size_t
pmatch_stack_length( potential_match_stack* pmatch );

 potential_match_stack*
add_pmatch( potential_match_stack* stack, 
            const void* node, const NODE_TYPE node_type, const enum STRAND stnd,
            const int node_level, const int max_num_levels,
            const float penalty, const float min_match_penalty, 
            const float max_penalty_diff );

potential_match
pop_pmatch( potential_match_stack* stack );

void* 
tree_malloc( size_t size );

void* 
tree_realloc( void* ptr, size_t size );

void 
tree_free( void* ptr );

void
print_dynamic_node( dynamic_node* node );

 void init_static_node( static_node** node );

 void init_dynamic_node( dynamic_node** node );

extern void init_index( struct index_t** index,
                        int seq_length,
                        int num_diploid_chrs );

 dynamic_node*
add_child_to_dynamic_node( 
    dynamic_node* node, 
    LETTER_TYPE bp,
    int loc,
    int num_letters
);

 void add_child_to_static_node( 
    static_node* node,
    LETTER_TYPE bp
);

void free_node( void* node, char node_type );

void free_node_and_children( void* node, char node_type );

void free_tree( struct index_t* root );

 void 
build_static_node_from_dynamic_node( dynamic_node* dnode, 
                                     static_node** snode   );

void
build_dynamic_node_from_sequence_node(  sequences_node* qnode, 
                                        dynamic_node** dnode,
                                        const int num_levels,
                                        LEVEL_TYPE level,
                                        struct genome_data* genome );

 int 
find_child_index_in_dynamic_node ( 
    dynamic_node* node, 
    LETTER_TYPE bp
);

extern void
index_genome( struct genome_data* genome );

/*
find_child_index_in_static_node is done in add_sequence ( it's a simple hash )
*/

void 
add_sequence( struct genome_data* genome,
              LETTER_TYPE* seq, const int seq_length,
              INDEX_LOC_TYPE genome_loc );

extern void
add_junction_positions_from_chr_string( 
    struct index_t* index, int seq_length, char* chr_str, int chr_index );

extern void 
naive_add_genome_from_fasta_file( 
    char* filename, int seq_length, static_node* root, int add_junctions );

extern void
find_matches_from_root(
        struct index_t* index,

        struct index_search_params* search_params,
        mapped_locations* results,

        struct genome_data* genome,

        /* the length of the two reads ( below ) */
        const int read_len,
        
        struct penalty_t* fwd_penalties,
        struct penalty_t* rev_penalties,

        enum bool only_find_unique_sequences );

size_t
size_of_snode( );

size_t
calc_node_size( NODE_TYPE type, void* node  );

size_t 
sizeof_tree( struct index_t* index );


/* fwd declaration for the index type */
struct index_t;

void
free_ondisk_index( struct index_t* index );

void
load_ondisk_index( char* index_fname, struct index_t** index );

void
build_ondisk_index( struct index_t* index, char* ofname  );


/*
 * Check to see that, for a random sequence, the penalties returned for every
 * item in the tree is equal to the penalty returned by calcualting the penalty
 * on each sequence in order. This tests both that every correct sequence is
 * being returned, and the the tree's calculated penalties are equal to the 
 * full sequence penalties. If this test passes regularly, then you can test
 * the penalty function by calling multiple_letter_penalty 
 */
int 
check_full_tree_penalty( char* genome_str, 
                         int seq_length,
                         float min_match_penalty,
                         unsigned int num_tests,
                         float* lookuptable_position,
                         float* inverse_lookuptable_position,
                         float* lookuptable_bp
);

#endif // INDEX_GENOME


