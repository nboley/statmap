/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>

#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <sys/stat.h> /* permission bit defines */

#include "quality.h"
#include "index_genome.h"
#include "genome.h"
#include "pseudo_location.h"
#include "diploid_map_data.h"
#include "error_correction.h"
#include "mapped_location.h"

#include "log.h"

/* a global pointer offset */
/* This must be added to all pointers */
static size_t index_offset = 0;

#define DONT_MMAP_INDEX

/******************************************************************************/
/* Functions for finding seqeunces in thre tree                               */
/******************************************************************************/

float
subseq_penalty(
        struct read_subtemplate* rst,
        int subseq_offset,
        int subseq_length,
        struct penalty_array_t* penalties
    )
{
    assert( subseq_offset + subseq_length <= rst->length );

    float penalty = 0;
    /* loop over subsequence */
    int pos;
    for( pos = subseq_offset;
         pos < subseq_offset + subseq_length; 
         pos++ )
    {
        int bp = rst->char_seq[pos];

        if( bp == 'N' || bp == 'n' )
        {
            penalty += N_penalty;
        } else {
            /* take the match penalty for this bp - assuming we can find
             * a perfect match to a given subsequence in the genome, which
             * match would then have the best penalty (considering what we know
             * about error rates in positions, for the error scores, etc. from
             * the error model)? */
            penalty += penalties->array[pos].penalties[bp_code(bp)];
        }
   }

    return penalty;
}

int
find_optimal_subseq_offset( 
    struct read_subtemplate* rst,
    
    int subseq_len,

    /* define region of underlying read to search for optimal subsequences */
    int region_start,
    int region_length
) {
    /* Make sure the search region makes sense */
    assert( region_start >= 0 && region_start <= (rst->length - subseq_len) );
    assert( region_length <= rst->length );
    assert( (region_start + region_length - subseq_len) >= 0 );
    
    /*
       error_prb returns the inverse log probability of error:
       log10(1 - P(error)) for matches
    */
    int optimal_offset = region_start;
    float max_so_far = -FLT_MAX;
        
    int i;
    /* each possible start bp in the subsequence */
    for( i = region_start;
         i < (region_start + region_length - subseq_len);
         i++ )
    {
        float ss_pen = subseq_penalty(rst, i, subseq_len, rst->fwd_penalty_array);
        if( ss_pen > max_so_far ) {
            max_so_far = ss_pen;
            optimal_offset = i;
        }
    }

    /* Return the offset of the optimal index probe from the start of the read
     * subtemplate */
    return optimal_offset;
};

struct indexable_subtemplate*
build_indexable_subtemplate(
        struct read_subtemplate* rst,
        struct index_t* index,
        
        // area of the read subtemplate to take an indexable subtemplate from
        int range_start,
        int range_length,
        
        enum bool choose_random_offset
    )
{
    int subseq_length = index->seq_length;

    /* Our strategy for choosing indexable subtemplates depends on the type of
     * assay. For gapped assays (RNA_SEQ), we wish to maximize the distance
     * between probes in order to optimize intron finding. For ungapped assays,
     * we want to use the highest quality subsequence in the read for the index
     * search. */
    int subseq_offset = 0;

    if( _assay_type == RNA_SEQ ) // Gapped assay
    {
        /* FIXME soft clipping for gapped assays? */
        if( range_start == 0 ) {
            /* subseq_offset = softclip_len ? will that mess up the matching
             * code later? */
            subseq_offset = 0;
        } else {
            subseq_offset = range_start + range_length - subseq_length;
        }
    } else {
        if( choose_random_offset ) {
            double rn = ((double) rand() / ((double) RAND_MAX));
            subseq_offset = range_start + (int) (rn*(range_length - subseq_length));
        } else {
            subseq_offset = find_optimal_subseq_offset(
                rst,
                subseq_length,
                range_start,
                range_length
            );
        }
    }

    assert( subseq_offset >= 0 );
    struct indexable_subtemplate* ist = NULL;
    init_indexable_subtemplate( &ist, rst, subseq_length, subseq_offset);
    
    return ist;
}

struct indexable_subtemplates*
build_indexable_subtemplates_from_read_subtemplate(
        struct read_subtemplate* rst,
        struct index_t* index,
        enum bool use_random_subtemplate_offset
    )
{
    int subseq_length = index->seq_length;
    
    /* Make sure we the read is long enough for us to build index probes
     * (considering any softclipped bases from the start of the rst) */
    int indexable_length = rst->length - softclip_len;
    if( indexable_length < subseq_length )
    {
        statmap_log( LOG_WARNING,
                "Probe lengths must be at least %i basepairs short, to account for the specified --soft-clip-length (-S)",
                softclip_len
            );
        return NULL;
    }

   struct indexable_subtemplates* ists = NULL;
   init_indexable_subtemplates( &ists );

    /* for now, try to build the maximum number of indexable subtemplates up to
     * a maximum */
    int num_partitions;
    if( _assay_type == RNA_SEQ ) {
        /* for RNA-seq, always use two probes at either end of the read so we
         * can maximize the space for finding introns */
        num_partitions = 2;
    } else {
        num_partitions = MIN( indexable_length / subseq_length,
                              MAX_NUM_INDEX_PROBES );        
    }
    int partition_len = floor((float)indexable_length / num_partitions);
    
    int i;
    for( i = 0; i < num_partitions; i++ )
    {
        struct indexable_subtemplate* ist = NULL;

        int partition_start = softclip_len + i * partition_len;

        /* partition the read into equal sized sections and try to find the
         * best subsequence within each section to use as an index probe */
        ist = build_indexable_subtemplate( 
            rst, index,
            partition_start, partition_len,
            use_random_subtemplate_offset );

        if( ist == NULL ) {
            free_indexable_subtemplates( ists );
            return NULL;
        }

       // copy indexable subtemplate into set of indexable subtemplates
       add_indexable_subtemplate_to_indexable_subtemplates( ist, ists );
       // free working copy
       free_indexable_subtemplate( ist );
    }
    
    return ists;
}

int
search_index(
        struct genome_data* genome,
        struct indexable_subtemplate* ist,
        struct index_search_params* search_params,
        mapped_locations** results,
        enum bool only_find_unique_sequence
    )
{
    // reference to index
    struct index_t* index = genome->index;

    int subseq_length = index->seq_length;
    
    /* prepare the results container */
    init_mapped_locations( results, ist );
    
    /* Build bitpacked copies of the fwd and rev strand versions of this
     * indexable subtemplate */

    /* Store a copy of the read */
    /* This read has N's replaced with A's, and might be RC'd */
    char* sub_read = calloc(subseq_length + 1, sizeof(char));
    assert( sub_read != NULL );
    /* note that the NULL ending is pre-set from the calloc */
    memcpy( sub_read, ist->char_seq,
            sizeof(char)*(subseq_length) );
    replace_ns_inplace( sub_read, subseq_length );

    /** Deal with the read on the fwd strand */
    /* Store the translated sequences here */
    LETTER_TYPE *fwd_seq;
    fwd_seq = translate_seq( sub_read, subseq_length );
    /* If we couldnt translate it */
    if( fwd_seq == NULL )
    {
        // fprintf(stderr, "Could Not Translate: %s\n", st->char_seq);
        return -1;
    }
    assert( fwd_seq != NULL );
    
    /** Deal with the read on the opposite strand */
    LETTER_TYPE *bkwd_seq;
    char* tmp_read = calloc(subseq_length + 1, sizeof(char));
    rev_complement_read( sub_read, tmp_read, subseq_length );
    free( sub_read );
    
    bkwd_seq = translate_seq( tmp_read, subseq_length );
    assert( bkwd_seq != NULL );
    free( tmp_read );
    
    /* search the index */
    int rv = find_matches_from_root(
            index, 

            search_params,
            *results,

            genome,

            /* length of the reads */
            subseq_length,
            
            ist->fwd_penalties,
            ist->rev_penalties,

            only_find_unique_sequence
        );

    /* Cleanup memory */
    free( fwd_seq );
    free( bkwd_seq );
    
    return rv;
};


int
search_index_for_read_subtemplate(
        struct read_subtemplate* rst,
        struct mapping_params* mapping_params,
        
        mapped_locations*** search_results,

        struct genome_data* genome,
        enum bool use_random_subtemplate_offset
    )
{
    // build a set of indexable subtemplates from this read subtemplate
    rst->ists = build_indexable_subtemplates_from_read_subtemplate(
        rst, genome->index, use_random_subtemplate_offset );

    /* if we couldn't build indexable sub templates, ie the read was too short, 
       then don't try and map this read */
    if( rst->ists == NULL ) {
        return CANT_BUILD_READ_SUBTEMPLATES_ERROR;
    }
        
    /* Stores the results of the index search for each indexable subtemplate */
    int search_results_length = rst->ists->length;
    *search_results = calloc(sizeof(mapped_locations*),search_results_length+1);
        

    /* initialize search parameters for the index probes */
    struct index_search_params* index_search_params
        = init_index_search_params(rst->ists, mapping_params);

    int num_valid_index_searches = 0;
    int i;
    for( i = 0; i < rst->ists->length; i++ )
    {
        int rv = search_index(
                genome,
                rst->ists->container + i,
                index_search_params + i,
                &((*search_results)[i]),
                false
            );
        // if the index search returned an error, skip this index probe
        // if there are too many candidate mappings, skip this index probe
        if( rv != 0 
            || (*search_results)[i]->length > MAX_NUM_CAND_MAPPINGS)  
        {
            free((*search_results)[i]->locations);
            (*search_results)[i]->locations = NULL;
            (*search_results)[i]->length = 0;
            (*search_results)[i]->allocated_length = 0;
            continue;
        }        
        num_valid_index_searches++;
    }    

    if(num_valid_index_searches < MIN_NUM_INDEX_PROBES)
    {
        free( index_search_params );
        free_search_results( *search_results );
        *search_results = NULL;
        
        return NOT_ENOUGH_VALID_INDEX_PROBES_ERROR;
    }

    free( index_search_params );    
    return 0;
}


/******************************************************************************/
/* Memory function for  building the index      */
/******************************************************************************/

void* 
tree_malloc( size_t size )
{
    void* rv;
    rv = malloc( size );
    if( rv == NULL ) {
        printf("Out of memory\n");
        exit(-1);
    }

    return rv;
}   

void* 
tree_realloc( void* ptr, size_t size )
{
    void* rv;
    rv = realloc( ptr, size );
    if( rv == NULL ) {
        printf("Out of memory\n");
        exit(-1);
    }
    
    return rv;
}   

void 
tree_free( void* ptr ) 
{
    free( ptr );
}

/******************************************************************************/
/* Functions for queuing reads                                                 */
/******************************************************************************/

void
init_pmatch_stack( potential_match_stack** stack )
{
    *stack = malloc( sizeof(potential_match_stack) );
    (*stack)->first_index = STACK_INITIAL_FIRST_INDEX;
    (*stack)->last_index = STACK_INITIAL_FIRST_INDEX;
    (*stack)->allocated_size = STACK_GROWTH_FACTOR;
    (*stack)->matches = malloc( STACK_GROWTH_FACTOR*sizeof(potential_match) );
}

void
free_pmatch_stack( potential_match_stack* stack )
{
    free( stack->matches );
    free(stack);
}

void
fprint_pmatch( FILE* fp, potential_match* m )
{
    fprintf(fp, "%p\t %c\t %i\t %.3f\t %.3f\t %u\n", 
            m->node, m->node_type, m->node_level, 
            m->penalty, m->scaled_penalty, m->strnd );

}

void
fprint_pmatch_stack( FILE* fp , potential_match_stack* s )
{
    unsigned int i;
    fprintf( fp, "====START===================================================\n");
    for(i = s->first_index; i < s->last_index; i++ )
    {
        fprintf( fp, "%i: ", i );
        fprint_pmatch( fp, s->matches + i );
    }
    fprintf( fp, "====END=====================================================\n");
}

 float 
scale_penalty( float penalty, // potential_match_stack* stack, 
               int level,
               int max_num_levels,
               float max_penalty_diff,
               float min_penalty )
{
    /* 
     * if we are at the top level, give an optimistic estimate until
     * we are proven wrong because we can probably find *something*
     * ( or else this is a stupid experiment ). For lack of a better estimate,
     * assume that it has an expected penalty of negative the max top difference.
     * 
     */
    double min_starting_penalty 
        = ( max_penalty_diff < -0.1) ? max_penalty_diff: min_penalty;
    
    double rv = ( level < 8 ) ? min_starting_penalty * ( ( 8-level ) / 8 ) : 0;
    rv +=  penalty*(((double) max_num_levels)/(level+1));

    return rv;
}

size_t
pmatch_stack_length( potential_match_stack* pmatch )
{
    return (pmatch->last_index - pmatch->first_index);
};

 potential_match_stack*
add_pmatch( potential_match_stack* stack, 
            const void* node, const NODE_TYPE node_type, const enum STRAND strnd,
            const int node_level, const int max_num_levels,
            const float penalty, const float min_match_penalty, const float max_penalty_diff)
{
    /* if necessary, increase the available size */
    // Assert guaranteed by type definition
    // assert( pmatch_stack_length( stack ) >= 0 );
    if( stack->last_index == stack->allocated_size )
    {
        stack->allocated_size += STACK_GROWTH_FACTOR;
        stack->matches = realloc( stack->matches, 
                                  sizeof(potential_match)*
                                  (stack->allocated_size) 
        );
    }

    /* if necessary, shift the stack forward */
    if( stack->first_index == 0 )
    {
        /* if there isnt enough room for the shift, realloc */
        if( pmatch_stack_length(stack) + STACK_INITIAL_FIRST_INDEX 
                >= stack->allocated_size )
        {
            stack->allocated_size += STACK_GROWTH_FACTOR;
            stack->matches = realloc( stack->matches, 
                                      sizeof(potential_match)*
                                      (stack->allocated_size) 
            );
            assert( stack->matches != NULL );
        }

        /* move the stack forward */
        memmove( stack->matches + STACK_INITIAL_FIRST_INDEX, /* dst */
                 stack->matches /* src */, 
                 sizeof(potential_match)*pmatch_stack_length( stack  )  
        );
        
        stack->first_index += STACK_INITIAL_FIRST_INDEX;
        stack->last_index += STACK_INITIAL_FIRST_INDEX;
    }

    assert( pmatch_stack_length(stack) == 0
            || (stack->matches + stack->last_index - 1)->node_type == 'd' 
            || (stack->matches + stack->last_index - 1)->node_type == 's'
            || (stack->matches + stack->last_index - 1)->node_type == 'l'
            || (stack->matches + stack->last_index - 1)->node_type == 'q' );

    /* TODO - make this smaller about already known penalties */
    /* 
     * If the penalty is less than the penalty at the top of 
     * the stack when weighted by their relative locations, add
     * this to the top.
     */

    potential_match* new_location = NULL;
    float scaled_penalty;
    scaled_penalty 
        = scale_penalty( penalty, node_level, max_num_levels, 
                         min_match_penalty, max_penalty_diff  );

    /* if the stack is currently empty */
    if( pmatch_stack_length(stack) == 0 )
    {
        new_location = stack->matches + stack->last_index;
        stack->last_index += 1;
    } else {
        float top_scaled_penalty 
            = (stack->matches + stack->last_index - 1)->scaled_penalty;
                
        if( scaled_penalty > top_scaled_penalty
            && scaled_penalty > min_match_penalty     ) 
        {
            /* add this to the top */
            new_location = stack->matches + stack->last_index;
            stack->last_index += 1;
        } else {
            /* add to the bottom */
            new_location = stack->matches + stack->first_index - 1;
            stack->first_index -= 1;
        }
    }

    assert( new_location != NULL );
    /* add the match to the set location */
    new_location->node = (void*) node;
    new_location->node_type = node_type;
    new_location->node_level = node_level;
    new_location->strnd = strnd;
    new_location->penalty = penalty;
    new_location->scaled_penalty = scaled_penalty;

    return stack;
}

potential_match
pop_pmatch( potential_match_stack* stack )
{
    /* TODO - lock access for thread safety */
    assert( pmatch_stack_length( stack ) > 0 );
    stack->last_index -= 1;
    return *(stack->matches + stack->last_index);
}

/* FIXME - WHERE SHOULD THIS GO */

void
print_dynamic_node( dynamic_node* node )
{
    int i;
    printf("#\tLet\tCt\tPtr\n");
    
    dynamic_node_child* children =
        get_dnode_children( node );

    for( i = 0; i < get_dnode_num_children( node ); i++ )
    {
        printf("%i:\t%i\t%c\t%p\n", 
               i+1, 
               children[i].type,
               children[i].letter, 
               children[i].node_ref
        );
    }

    return;
}


/******************************************************************************/
/* Functions for the construction of the tree                                 */
/******************************************************************************/

/*
    unsigned char level;
    unsigned char node_type;
    * d: dynamic branches
    * s: static branches
    * q: sequences
    unsigned short children_length;
    void* children;
*/

 void init_static_node( static_node** node )
{
    /* note that the calloc initializes the children to NULL */
    /* FIXME check for NULL */
    
    *node = calloc( ALPHABET_LENGTH, sizeof(static_node_child) );

    #ifdef PROFILE_MEMORY_USAGE
    num_static_nodes++;
    #endif
}

 void init_dynamic_node( dynamic_node** node )
{
    *node = tree_malloc(sizeof(unsigned short));
    set_dnode_num_children( *node, 0 );
    
    #ifdef PROFILE_MEMORY_USAGE
    num_dynamic_nodes++;
    #endif
}

unsigned short
get_dnode_num_children( const dynamic_node* const dnode )
{
    return *((unsigned short*) dnode );
}

void
set_dnode_num_children( dynamic_node* const dnode, 
                        const unsigned short value )
{
    *((unsigned short*) dnode ) = value;
}

dynamic_node_child*
get_dnode_children( const dynamic_node* const dnode )
{
    return ( dynamic_node_child* ) ( dnode + sizeof( unsigned short) ) ;
}

size_t
size_of_dnode( const dynamic_node* const node )
{
    return sizeof( unsigned short ) 
        + get_dnode_num_children( node )*sizeof( dynamic_node_child );
}


extern void 
init_index( struct index_t** index, int seq_length, int num_diploid_chrs )
{   
    /* init tree */
    *index = malloc( sizeof( struct index_t ) );
    static_node* tree_root;
    init_static_node( &tree_root );
    (*index)->index = tree_root;
    (*index)->index_type = TREE;
    (*index)->seq_length = seq_length;

    /* init ps_locs */
    (*index)->ps_locs = NULL;
    init_pseudo_locations( &((*index)->ps_locs) );

    /* init diploid_maps */
    (*index)->diploid_maps = NULL;
    init_diploid_maps_t( &((*index)->diploid_maps), num_diploid_chrs );
}

 dynamic_node*
add_child_to_dynamic_node( 
    dynamic_node* node, 
    LETTER_TYPE bp,
    int loc,
    int num_letters
)
{
    /* allocate new memory for the child */
    /* The memory check is done in tree_realloc */
    node = tree_realloc(
        node, size_of_dnode(node) + sizeof(dynamic_node_child)
    );
    
    /* update the number of children */
    int num_children
        = get_dnode_num_children( node ) + 1;
    set_dnode_num_children( node, num_children );

    /* get a pointer to the children array */
    dynamic_node_child* children =
        get_dnode_children( node );
    
    
    #ifdef PROFILE_MEMORY_USAGE
    dynamic_children_bytes += 
        ( sizeof(char) + sizeof(LETTER_TYPE) + sizeof(void*));
    #endif

    /* 
     * if the loc is not at the end of the array, copy the new memory one child
     * forward in the list, to make a place for the new entry.
     */
    if( loc != ( num_children - 1 ) ) {
        /* the children array */
        memmove( children+loc+1, /* destination */
                 children+loc, /* source */
                 sizeof(dynamic_node_child)*(num_children-loc-1) /* size */
        );
    }

    /* If the new child is a locations node */
    if( num_letters - 1 == 0 )
    {
        assert( num_letters > 0 );
        /* set the node type of the child */
        children[loc].type = 'l';
        /* set the child type ( ie, in the 1 bp length letter, the bp type ) */
        children[loc].letter = bp;
        /* initialize the new child node */
        locations_node* new_node;
        init_locations_node( &new_node );
        children[loc].node_ref = new_node;
    } else {
        assert( num_letters > 0 );
        /* set the node type of the child */
        children[loc].type = 'q';
        /* set the child type ( ie, in the 1 bp length letter, the bp type ) */
        children[loc].letter = bp;
        /* initialize the new child node */
        sequences_node* new_node;
        init_sequences_node( &new_node );
        children[loc].node_ref = new_node;
    }
    
    return node;
}

 void add_child_to_static_node( 
    static_node* node,
    LETTER_TYPE bp
)
{
    /* set the node type of the child */
    node[bp].type = 'q';
    /* initialize the new child node */
    sequences_node* new_node;
    init_sequences_node( &new_node );
    node[bp].node_ref = new_node;
}


/**** CLEANUP functions *******************************************************/

void free_node( void* node, char node_type )
{
    switch( node_type )
    {
        case 's':
            tree_free( node );
            #ifdef PROFILE_MEMORY_USAGE
            num_static_nodes--;
            #endif
            break;
        case 'd':
            tree_free( node );
            #ifdef PROFILE_MEMORY_USAGE
            num_dynamic_nodes--;
            #endif
            break;
        case 'l':
            tree_free( node );
            #ifdef PROFILE_MEMORY_USAGE
            num_locations_nodes--;
            #endif
            break;
        case 'q':
            free_seqs( (sequences_node*) node );
            #ifdef PROFILE_MEMORY_USAGE
            num_sequence_nodes--;
            #endif
            break;
        default:
            statmap_log( LOG_FATAL, "Impossible branch '%c' reached in free_node.", node_type );
    }
}

void free_node_and_children( void* node, char node_type )
{
    node += index_offset;

    int i;
    switch( node_type )
    {
        case 's':
            /* free the children */
            for( i = 0; i < ALPHABET_LENGTH; i++ )
            {
                if( ((static_node*) node)[i].node_ref != NULL )
                {
                    free_node_and_children( 
                        ((static_node*) node)[i].node_ref,
                        ((static_node*) node)[i].type
                    );
                }
            }
            free_node( node, 's' );
            break;
        case 'd':
            /* free the children */
            for( i = 0; i < get_dnode_num_children( node ); i++ )
            {
                free_node_and_children( 
                    get_dnode_children( node )[i].node_ref,
                    get_dnode_children( node )[i].type
                );
            }
            free_node( node, 'd' );
            break;
        case 'q':
            free_seqs( (sequences_node*) node );
            break;
        case 'l':
            free_node( node, 'l' );
            break;
        default:
            statmap_log( LOG_FATAL, "Impossible branch '%c' reached in free_node_and_children.", node_type );
    }
}

void free_tree( struct index_t* index )
{
    assert( index->index_type == TREE );
    free_node_and_children( (static_node*) index->index, 's');

    //free_diploid_maps_t( index->diploid_maps );

    free( index );
}

/**** END CLEANUP functions ***************************************************/

 void 
build_static_node_from_dynamic_node( dynamic_node* dnode, 
                                     static_node** snode   )
{

    /* 
     * since it's a static node, we allocate the array 
     * ( noting calloc zeros the memory which sets the pointers to NULL ) 
     */
    *snode = calloc(ALPHABET_LENGTH, sizeof(void*));

    dynamic_node_child* children 
        = get_dnode_children( dnode );

    /* populate the children from the dynamic node */
    int i;
    for( i = 0; i < get_dnode_num_children( dnode ); i++ )
    {
        /* TODO - micro optimization ( is this optimized out? ) */
        (*snode)[children[i].letter].node_ref
            = children[i].node_ref;

        (*snode)[children[i].letter].type
            = children[i].type;
    }

    return;
}

void
build_dynamic_node_from_sequence_node(  sequences_node* qnode, 
                                        dynamic_node** dnode,
                                        const int num_levels,
                                        LEVEL_TYPE level,
                                        struct genome_data* genome )
{
    int i;
    
    init_dynamic_node(dnode);

    /* 
     * TODO - what the below says - for now we do it naively
     * populate the children from the sequence node
     *
     * NOT ACTUALLY IMPLEMENTED - THIS IS A COMMENT FOR A TODO
     * we could do this using the machinery that we've already developed, 
     * but it will be much faster if we can just peak inside of the sequence 
     * node, so that is what we do. Don't replicate this code for adding
     * sequences to a dynamic node - call add_child_to_dynamic_node with
     * find_child_index_in_dynamic_node
     * END NOT ACTUALLY IMPLEMENTED
     *
     */
    
    /* the number of letters in each tailing sequence in seqs */
    int num_letters = num_levels - level;

    /* make sure it makes sense to convert the node - TODO better error check */
    assert( get_num_sequence_types(qnode) > 1 );

    /* get the array of sequences form the (packed) sequences node */
    LETTER_TYPE* sequences = get_sequences_array_start( qnode, num_letters );

    /* loop through all of the children */
    for( i = 0; i < get_num_sequence_types(qnode); i++ )
    {
        /* 
         * Note that the sequences store all of the letters in the sequence 
         * node in sequential order. Since we only care about the first letter
         * in each sequence, we start i at 0 and take the i*num_letters letter
         */
        LETTER_TYPE bp = sequences[i*num_letters];

        /* 
         * since the sequence node is sorted, we dont have to find the child,
         * we just need to check to see if it's at the end. If it is, then 
         * add the sequences to the child. If not, then add a new node and 
         * then add the sequences.
         */        
        int clen = get_dnode_num_children( *dnode );
        dynamic_node_child* children = get_dnode_children( *dnode );
        /* 
         * if there are no children, or this child has a first letter that is 
         * different than the first letter of the last child, then add a child
         * to the dynamic node in the last position
         */
        if( clen == 0 ||
            bp != children[clen-1].letter )
        {
            /* BUG!!!! fix this fn */
            *dnode = add_child_to_dynamic_node( *dnode, bp, clen, num_letters );
            /* update children in case the realloc moved the pointer */
            children = get_dnode_children( *dnode );
        } else {
            /* make clen the correct child index */
            clen = clen - 1;
        }

        /**** Add the sequence to the leaf *********************************/
        /* TODO - cleanup this code block */
        /* If we are at a locations node, then the seq length is 0  */
        if( num_letters - 1  == 0 )
        {
            /* 
             * Add the sequence to the leaf - we just assume that it is a genome 
             * location for now. After it's added. we check the bit to make sure 
             * and, if it's a pointer, we fix it.
             */
            locations_node** child_seqs  = 
                (locations_node**) &((children + clen)->node_ref);

            /* 
             * Now, check if there are multiple genome locations associated
             * with this sequence.
             */
            locs_union loc = 
                get_genome_locations_array_start( qnode, num_letters )[i];


            if( !check_sequence_type_ptr(qnode, i) )
            {
                *child_seqs = add_location_to_locations_node(   
                    *child_seqs, 
                    loc.loc
                );
            } 
            else 
            /* Add all of the locs in the referenced array */
            {
                INDEX_LOC_TYPE* gen_locs 
                    = get_overflow_genome_locations_array_start( 
                        qnode,  num_letters )
                    + loc.locs_array.locs_start;

                int j;
                for(j = 0; j < loc.locs_array.locs_size; j++ )
                {
                    *child_seqs = add_location_to_locations_node(   
                        *child_seqs, 
                        gen_locs[j]
                    );
                }
            }
        }
        /* Otherwise, this is a sequence node */
        else {
            /* 
             * Add the sequence to the leaf - we just assume that it is a genome 
             * location for now. After it's added. we check the bit to make sure 
             * and, if it's a pointer, we fix it.
             */
            sequences_node** child_seqs  = 
                (sequences_node**) &((children + clen)->node_ref);

            /* 
             * Now, check if there are multiple genome locations associated
             * with this sequence.
             */
            locs_union loc = 
                get_genome_locations_array_start( qnode, num_letters )[i];

            /* 
             * if the correct bit is set, then the genome location is a pointer.
             */        
            if( !check_sequence_type_ptr(qnode, i) )
            {
                *child_seqs = add_sequence_to_sequences_node(
                    /* we pass NULL for the pseudo_locs, because we know that
                       it is only used to add pseudo locations and, since we 
                       are building a dynamic node from a sequences node, 
                       anything that is already a pseudo loc wont change */
                    genome,
                    NULL,
                    *child_seqs, 
                    sequences + i*num_letters + 1, 
                    num_letters-1, 
                    loc.loc
                );
            } 
            else 
            /* Add all of the locs in the referenced array */
            {
                INDEX_LOC_TYPE* gen_locs 
                    = get_overflow_genome_locations_array_start( qnode,  num_letters )
                    + loc.locs_array.locs_start;

                int j;
                for(j = 0; j < loc.locs_array.locs_size; j++ )
                {
                    INDEX_LOC_TYPE loc = gen_locs[j];

                    *child_seqs = add_sequence_to_sequences_node(   
                        genome,
                        NULL,
                        *child_seqs, 
                        sequences + i*num_letters + 1, 
                        num_letters-1, 
                        loc
                    );
                }
            }
        }
    }

    return;
}

 int 
find_child_index_in_dynamic_node ( 
    dynamic_node* node, 
    LETTER_TYPE bp
) 
{
    /* 
     * find the correct insert location - we use a binary sort to keep 
     * the items in order - we dont use the stdlib because of the overhead.
     */
    int low = 0;
    int high = get_dnode_num_children( node );
    dynamic_node_child* children = get_dnode_children( node );

    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        int comparison = cmp_letters( children[mid].letter, bp );
        if( comparison < 0 ) {
            low = mid + 1; 
        } else {
            /* can't be high = mid-1: here A[mid] >= value, */
            /* so high can't be < mid if A[mid] == value    */
            high = mid; 
        }
    }

    return low;
}

/*
find_child_index_in_static_node is done in add_sequence ( it's a simple hash )
*/


void 
add_sequence( struct genome_data* genome,
              LETTER_TYPE* seq, const int seq_length,
              INDEX_LOC_TYPE genome_loc ) 
{
    /* direct references to index and pseudo_locations */
    struct index_t* index = genome->index;
    struct pseudo_locations_t* ps_locs = index->ps_locs;

    assert( index->index_type == TREE );
    static_node* root = index->index;

    const int num_levels = calc_num_letters( seq_length );
    LEVEL_TYPE level = 0;
    /* 
     * store a pointer to the address that references curr_node. Since there is no
     * node pointing to root, it's initialized to NULL. This shouldn't matter because
     * it's purpose is to keep a record for node conversions, and root is always static
     * and thus never converted.
     */
    void** node_ref = NULL;
    NODE_TYPE* node_type_ref = NULL; 

    void* curr_node = (void*) root;
    char curr_node_type = 's';

    /* traverse child nodes until we hit a leaf node */
    while( level < num_levels )
    {
        LETTER_TYPE bp = seq[level];
        if( curr_node_type == 's' )
        {
            /* check to make sure this child exists - if it doesn't, create it*/
            if( ((static_node*) curr_node)[bp].node_ref == NULL ) 
            {
              add_child_to_static_node( (static_node*) curr_node, bp );
            }
            
            /* FIXME - but this is ok for now */
            node_type_ref = &(((static_node*) curr_node)[bp].type);
            curr_node_type = *node_type_ref;

            node_ref = &(((static_node*) curr_node)[bp].node_ref);
            curr_node = *node_ref;
            
            // debug
            // printf("LeveL: %i\tStatic Node: %p %p\t%c\n", 
            //       level, node_type_ref, node_ref, *node_type_ref);
            
        } else if( curr_node_type == 'd' )
        {
            int child_index = find_child_index_in_dynamic_node(
                (dynamic_node*) *node_ref, bp);
            /* 
             * make sure that the child is of the correct type - child_index 
             * might just be the insert location. If it is, it's the current 
             * node. If not, we need to build the new node.
             */
            
            int num_children = get_dnode_num_children( *node_ref );
            dynamic_node_child* children = get_dnode_children( *node_ref );

            /* if the node is incorrect, then insert the correct child */
            if ( child_index == num_children
                 || children[child_index].letter != bp )
            {
              *node_ref = add_child_to_dynamic_node( 
                  *node_ref, bp, child_index, num_levels - level
              );
            }

            assert( child_index < get_dnode_num_children(*node_ref) );
            assert( child_index >= 0 );
            
            node_type_ref = 
                &( get_dnode_children(*node_ref)[child_index].type);
            curr_node_type = *node_type_ref;

            node_ref = &( get_dnode_children(*node_ref)[child_index].node_ref );
            curr_node = *node_ref;

            // debug
            // printf("LeveL: %i\tDynamic Node: %p %p\t%c\n", 
            //       level, node_type_ref, node_ref, *node_type_ref);
            

        /* 
         * If the node type is sequence or locations, then we have reached 
         * the bottom of the tree 
         */
        } else { 
            assert( curr_node_type == 'q' || curr_node_type == 'l' );
            break;
        } 
        
        /* move to the next level */
        level++;

        assert( curr_node_type == 's' 
                || curr_node_type == 'q' 
                || curr_node_type == 'd'
                || curr_node_type == 'l' );

        assert( *node_type_ref = curr_node_type );    
        assert( *node_ref = curr_node );    
    }

    /* If we are at a level node */
    if( curr_node_type == 'l' )
    {
        // debug
        // printf("LeveL: %i\tLocations Node: %p %p\t%c\n", 
        //       level, node_type_ref, node_ref, *node_type_ref);
        
        /* 
         * Location nodes dont store any sequence information 
         * thus, if this is a location node we better be at the
         * bottom level of the tree
         */
        assert( level == num_levels );
        
        *node_ref = add_location_to_locations_node( 
            (locations_node*) curr_node, genome_loc );
        
        return;
    }

    /* make sure that we are at a sequence node */

    assert( curr_node_type == 'q' );
    assert( *node_type_ref == 'q' );
    /* add the sequence to current node */
    *node_ref = add_sequence_to_sequences_node(   
        genome,
        ps_locs,
        (sequences_node*) *node_ref, 
        seq + level,
        num_levels - level,
        genome_loc                  
    );

    /* 
     * check if the current sequence node has too many entries. If it does, 
     * then convert it to a dynamic node. 
     */

    /* If we are at the final level, then this should be a locations node */
    assert( level < num_levels );
    
    /* If the sequence node has too many entries, then convert it to
       a dynamic node. */

    /* 
     * If the first letter of every child of the sequence node is the 
     * same, then the child of the resulting dynamic node will have too
     * many children. As such, we need to keep splitting until the dynamic
     * node has more than one child.
     *
     * After the first split, if num_sequence_types == MAX_SEQ_NODE_ENTRIES
     * then we know that the first letters were all the same. Therefore, 
     * at the end of this while loop we set curr_node to be the first child
     * of the resulting dynamic node. 
     *
     * We need to recheck level because, if we have iterated to the bottom,
     * the first child will be a locations node.
     *
     */

    int num_sequence_types = 
        get_num_sequence_types( (sequences_node*) *node_ref );

    assert( num_sequence_types <= MAX_SEQ_NODE_ENTRIES );
    
    while( level < num_levels 
           && num_sequence_types == MAX_SEQ_NODE_ENTRIES )
    {
        // debug
        // printf("Level: %i\tSequences Node: %p %p\t%c\n", 
        //       level, node_type_ref, node_ref, *node_type_ref);

        
        /* Make sure the node is a sequence node */
        assert( 'q' == *node_type_ref );
        
        /* Make a dynamic node from *node_ref */
        dynamic_node* new_node = NULL; 
        build_dynamic_node_from_sequence_node( 
            (sequences_node*) *node_ref, 
            &new_node, 
            num_levels, 
            level,
            genome
        );

        /* Free the old sequences node */
        free_seqs( (sequences_node*) *node_ref );
        
        /* Set the pointers to point to the newly created node */
        *node_ref = new_node;
        *node_type_ref = 'd';
        
        /* get the first child */
        node_ref = &get_dnode_children( new_node )[0].node_ref;
        /* get the pointer to the first child's node type */
        node_type_ref = &get_dnode_children( new_node )[0].type;
        /* move to the next level */
        level += 1;
        
        if( *node_type_ref == 'q' )
        {
            num_sequence_types = 
                get_num_sequence_types( (sequences_node*) *node_ref );
        }

    assert( *node_type_ref == 'l' 
            || num_sequence_types < MAX_SEQ_NODE_ENTRIES );

    } 
    
    return;
}


int 
find_matches( void* node, NODE_TYPE node_type, int node_level, 
              const int seq_length,
              float curr_penalty, 

              struct index_search_params* search_params,
              mapped_locations* results,

              struct genome_data* genome,

              struct penalty_t* fwd_penalties,
              struct penalty_t* rev_penalties,

              /*
                 we pass this flag all the way down to optimize the error data
                 bootstrap by terminating early on multimappers
               */
              enum bool only_find_unique_sequences
    )
{
    // Initialize the return value to 0, for no error
    int rv = 0;
    
    const int num_letters = calc_num_letters( seq_length );

    /* Unpack the search parameters */
    float min_match_penalty = search_params->min_match_penalty;
    float max_penalty_spread = search_params->max_penalty_spread;

    /* initialize the stack */
    potential_match_stack* stack;
    init_pmatch_stack( &stack );
    stack = add_pmatch( stack, node, node_type, FWD,
                        node_level, num_letters, 
                        curr_penalty, 
                        min_match_penalty, max_penalty_spread );

    stack = add_pmatch( stack, node, node_type, BKWD,
                        node_level, num_letters, 
                        curr_penalty, 
                        min_match_penalty, max_penalty_spread );

    // get current time
    struct timespec start;
    struct timespec stop;
    int err;
    err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
    assert( err == 0 );
    
    int cntr;
    for(cntr = 0; pmatch_stack_length( stack ) > 0; cntr++ )
    {
        if( cntr%1000 == 0 )
        {
            double elapsed;
            err = clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stop);
            elapsed = (stop.tv_sec - start.tv_sec);
            elapsed += (stop.tv_nsec - start.tv_nsec) / 1000000000.0;
            if( elapsed > MAX_SEARCH_TIME )
            {
                results->length = 0;
                rv = INDEX_SEARCH_TOOK_TOO_LONG_ERROR;
                goto cleanup;            
            }
        }
        potential_match match = pop_pmatch( stack );

        /* add the index offset so that this is a phsyical location */
        /* 
         * because every pointer that we consider passes through this 
         * location for index searches, this is the only time that we
         * need consider it.
         *
         */
        void* node = index_offset + match.node;

        NODE_TYPE node_type = match.node_type; 
        int node_level = match.node_level;
        float curr_penalty = match.penalty;
        enum STRAND strnd = match.strnd;
        
        /* select sequence and penalty_array depending on fwd/bkwd strand */
        struct penalty_t* penalties;
        if( strnd == FWD ) {
            penalties = fwd_penalties;
        } else {
            penalties = rev_penalties;
        }

        /* TODO */
        /* 
         * If the sequence_length is less than the current node level, then
         * that means the index seed length is longer than the seq length. 
         * If this is the case, add every child below the current node as 
         * a match.
         */
           
        /* check to make sure that this branch is still valid */
        /* avoid rounding errors with the small number - it is your fault
           if you put in a negative number for the max penalty spread */
        if( max_penalty_spread > -0.01 &&
            curr_penalty < min_match_penalty )
        {
            continue;
        }
        
        if( node_type == 's')
        {
            /* deal with static node */
            /* TODO - optimize this case */
            unsigned int letter;
            /* loop through each potential child */
            for( letter = 0; letter < ALPHABET_LENGTH; letter++ )
            {
                /* if the child is null, keep going */
                if( ((static_node*) node)[letter].node_ref == NULL )
                    continue;

                /* this should be optimized out */
                float penalty_addition = compute_penalty( letter, node_level,
                        seq_length, min_match_penalty - curr_penalty,
                        penalties );

                /* if this letter exceeds the max, continue */
                if( penalty_addition > 0 ) {
                    /* skip forward to the next possible letter */
                    // int num_to_skip = ( 3 << ( LETTER_LENGTH - ((int) penalty_addition - 1 ) )
                    continue;
                }
                /* otherwise, find the penalty on the child function */
                else {
                  /* // debugging code
                    printf("Level: %i\t Node Type: %c\t Letter: %i\tPenalty: %e\n", 
                           node_level, node_type, letter, curr_penalty + penalty_addition);
                  */

                    /* add this potential match to the stack */
                    stack = add_pmatch( 
                        stack, 
                        ((static_node*) node)[letter].node_ref,
                        ((static_node*) node)[letter].type,
                        strnd,
                        node_level+1, num_letters, 
                        curr_penalty + penalty_addition, 
                        min_match_penalty, max_penalty_spread
                    );
                }
            }
        }
        /* deal with dynamic nodes */
        else if( node_type == 'd')
        {
            /* hopefully this will be optimized out */
            int num_children = get_dnode_num_children( node );
            dynamic_node_child* children = get_dnode_children( node );

            int i;
            /* loop through each potential child */
            for( i = 0; i < num_children; i++ )
            {
                /* this should be optimized out */
                LETTER_TYPE letter = children[i].letter;

                float penalty_addition = compute_penalty( letter, node_level,
                        seq_length, min_match_penalty - curr_penalty,
                        penalties);

                /* 
                 * If compute_penalty returns a value >= 1
                 * ( we say > 0.5 to deal with rounding errors )
                 * then that means that the letter at the level returned,
                 * minus 1 is the letter that exceeded the penalty. As
                 * such, since the letters are in sorted order, we dont 
                 * need to go back to penalty function. Rather, we can keep
                 * skipping letters until the letter at the returned level
                 * or any letter above that level changes. The code
                 * below implements this optimiztion.
                 */
                if( penalty_addition > 0.5 ) {
                    /* FIXME - consider the performace implications of this */
                    int break_index = (int) penalty_addition + 0.5;
                    if( break_index == 1 ) {
                        while( i + 1 < num_children &&
                               ((children[i].letter)&3) == 
                               ((children[i+1].letter)&3)  )
                        {
                            i++;
                        }
                    } else if( break_index == 2 ) {
                        while( i + 1 < num_children &&
                               ((children[i].letter)&15) == 
                               ((children[i+1].letter)&15)  )
                        {
                            i++;
                        }
                    } else if( break_index == 3 ) {
                        while( i + 1 < num_children &&
                               ((children[i].letter)&63) == 
                               ((children[i+1].letter)&63)  )
                        {
                            i++;
                        }
                    } 

                    continue;
                }
                /* otherwise, find the penalty on the child function */
                else {
                  /* debugging code
                    printf("Level: %i\t Node Type: %c\t Letter: %i\tPenalty: %e\n", 
                           node_level, node_type, letter, curr_penalty + penalty_addition);
                  */
                    stack = add_pmatch( stack, 
                                children[i].node_ref,
                                children[i].type,
                                strnd,
                                node_level+1, num_letters, 
                                curr_penalty + penalty_addition,
                                min_match_penalty, max_penalty_spread
                    );
                }
            }
        } 
        /* If this is a locations node */
        else if( node_type == 'l')
        {
            get_locations_from_locations_node(
                node, 
                results,
                curr_penalty,
                strnd,
                genome
            );

            /* 
             * if we are using the max penalty spread criteria, and 
             * the newest penalty 0 the amx spread is greater than 
             * the current minimum penalty, update it.
             * 
             */
            if( max_penalty_spread > -0.1 &&
                curr_penalty - max_penalty_spread > min_match_penalty )
            {
                min_match_penalty = curr_penalty - max_penalty_spread;
            }

        }
        /* deal with the sequence nodes */
        else {
            assert( node_type == 'q' );
            
            /* keep track of how many rsults we currently have */
            int old_results_len = results->length;
            
            /* debugging
            printf("Level: %i\t Node Type: %c\t\t\tPenalty: %e\n", 
                   node_level, node_type, curr_penalty );
            */
            float max_added_penalty = 
            find_sequences_in_sequences_node( 
                node, 
                curr_penalty, min_match_penalty,
                seq_length, num_letters, node_level, strnd,
                results,
                penalties,
                genome
            );            
            
            /* 
               If we are only looking for unique sequences,
               and we have already found some locations, then 
               finding more locations implies that we have 
               found some differing sequence ( because we know that
               each node stores identical sequence ). This is not
               to say that there are not non-identical sequence 
               within results already, but adding any more will 
               certainly be different. Thus, we set the results length
               to 0 and continue. 
             */
            if( only_find_unique_sequences 
                && old_results_len > 0 
                && results->length > old_results_len )
            {
                /* we can set this to zero and not worry about
                   a memory leak because the allocated length
                   will still be non-zero, so we can clean this
                   up properly */
                results->length = 0;
                goto cleanup;
            }
            
            // DEBUG
            //fprintf( stderr, "max_added_penalty: %f\n", max_added_penalty);

            /* 
             * if we are using the max penalty spread criteria, and 
             * the newest penalty 0 the amx spread is greater than 
             * the current minimum penalty, update it.
             * 
             */
            if( max_penalty_spread > -0.1 &&
                max_added_penalty - max_penalty_spread > min_match_penalty )
            {
                min_match_penalty = max_added_penalty - max_penalty_spread;
            }

            /* debugging
            print_mapped_locations( results );
            */
        }
    }

cleanup:
    free_pmatch_stack( stack );
    return rv;
}


extern int
find_matches_from_root(
        struct index_t* index,

        struct index_search_params* search_params,
        mapped_locations* results,

        struct genome_data* genome,

        /* the length of the two reads ( below ) */
        const int read_len,
        
        struct penalty_t* fwd_penalties,
        struct penalty_t* rev_penalties,

        enum bool only_find_unique_sequences
)
{
    assert( index->index_type == TREE );
    assert( index->seq_length == read_len );
    static_node* root = index->index;

    return find_matches( (void*) root, 's', 0, 
                         index->seq_length,
                         
                         0,
                         
                         search_params,
                         results,

                         genome,
                         
                         fwd_penalties,
                         rev_penalties,

                         only_find_unique_sequences ); 
}

size_t
size_of_snode( )
{
    return sizeof(static_node_child)*(1<<(2*LETTER_LEN));

}

size_t
calc_node_and_children_size( NODE_TYPE type, void* node  )
{
    size_t size = 0;
    int i;

    if( node == NULL )
        return 0;
    
    node += index_offset;
    
    switch ( type ) 
    {
    case 'q':
        /* sequence nodes dont have any children */
        return get_num_used_bytes( node );
    case 'l':
        /* locations nodes dont have any children */
        return size_of_locations_node( node );
    case 's':
        for( i = 0; i < 1<<(2*LETTER_LEN); i++ )
        {
            size += calc_node_and_children_size( 
                ((static_node*)node)[i].type, 
                ((static_node*)node)[i].node_ref 
            );
        }
        return size + size_of_snode();
    case 'd':
        for( i = 0; i < get_dnode_num_children( node ); i++ )
        {
            size += calc_node_and_children_size( 
                get_dnode_children(node)[i].type, 
                get_dnode_children(node)[i].node_ref 
            );
        }
        return size + size_of_dnode( node );
    }

    statmap_log( LOG_FATAL, "Unrecognized Node Type: '%c'",  type );
    assert( false );
    exit( -1 );
}

size_t
calc_node_size( NODE_TYPE type, void* node  )
{
    if( node == NULL )
        return 0;
    
    node += index_offset;
    
    switch ( type ) 
    {
    case 'q':
        /* sequence nodes dont have any children */
        return get_num_used_bytes( node );
    case 'l':
        /* locations nodes dont have any children */
        return size_of_locations_node( node );
    case 's':
        return size_of_snode();
    case 'd':
        return size_of_dnode( node );
    }

    statmap_log( LOG_FATAL, "Unrecognized Node Type: '%c'",  type );
    assert( false );
    exit( -1 );
}


size_t 
sizeof_tree( struct index_t* index )
{
    assert( index->index_type == TREE );
    return calc_node_and_children_size( 's', index->index );
}

/*****************************************************************************
 * Code to write this to disk 
 *
 * the rough algorithm is:
 *
 * 1) Build the tree in memory
 * 2) Calculate the size
 * 3) Open the mmap file
 * 4) Grow the file size
 * 5) MMAP it
 * HOW TO WRITE A NODE TO A STACK - 
   Add a new node.
   For each child ( if applicable ) 
      1) Add the child to the end of the stack
      2) Change the parent's pointer to NULL
      3) Add the pointer to the parents pointer to the child
 * 6) Write the tree to the file
 *    1) Build a to-do stack
 *    2) Write the root to the stack
 *    3) While the stack is not empty
 *        1) Write it
 *        2) Update the parent's pointer
 *        3) Append the children to the stack, with ptrs
 * 7) fsync
 * 8) free the mmapped tree
 * 9) free the in memory tree
 *
 */

/* On Disk Index Stack Items */
typedef struct {
    void** node_ref;
    NODE_TYPE node_type;
    LEVEL_TYPE level;
} ODI_stack_item;

/* On Disk Index Stack Items */
struct ODI_stack {
    size_t allocated_size;
    size_t size;
    ODI_stack_item* stack;
};

#define ODI_stack_GROW_FACTOR 100000

void
init_ODI_stack( struct ODI_stack** stack )
{
    *stack = malloc( sizeof( struct ODI_stack ) );
    (*stack)->allocated_size = ODI_stack_GROW_FACTOR;
    (*stack)->size = 0;
    (*stack)->stack = 
        malloc( sizeof( ODI_stack_item )*ODI_stack_GROW_FACTOR );
    
    return;
}

void
free_ODI_stack( struct ODI_stack* stack )
{
    free( stack->stack );
    free( stack );
}

void
add_ODI_stack_item( struct ODI_stack* stack, 
                    void** node_ref, 
                    NODE_TYPE node_type, 
                    LEVEL_TYPE level )
{
    /* 
     * It's possible for a root node to pass a NULL ptr ( indicating a 
     * non-existent child ). If this happens, do nothing.
     */
    if( *node_ref == NULL ) return;

    stack->size++;
    
    if( stack->size == stack->allocated_size )
    {
        stack->allocated_size += ODI_stack_GROW_FACTOR;
        stack->stack = realloc( 
            stack->stack, sizeof(ODI_stack_item)*(stack->allocated_size) );
        
        if( stack->stack == NULL )
        {
            statmap_log( LOG_FATAL, "Error allocating memory for the ODI Stack." );
            exit( -1 );
        }
    }
    
    size_t index = stack->size - 1;
    stack->stack[index].node_ref = node_ref;
    stack->stack[index].node_type = node_type;
    stack->stack[index].level = level;
    
    return;
}

int
cmp_ODI_stack_item( const ODI_stack_item* a,
                    const ODI_stack_item* b )
{
    return b->level - a->level;
} 

void
sort_ODI_stack( struct ODI_stack* stack )
{
    qsort ( 
        stack->stack,
        stack->size, 
        sizeof(ODI_stack_item), 
        (int(*)( const void*, const void* ))cmp_ODI_stack_item
    );
    
    return;
}

/* the size of the index header */
/* MAGIC_NUMBER + SEQ_LEN + INDEX_SIZE */
#define HEADER_SIZE ((size_t)( 1 + 1 + sizeof(size_t) ))

void
free_ondisk_index( struct index_t* index ) {
    if( index->index != NULL ) {
        free( index->index );
    } else {
        #ifndef MMAP_INDEX
        free( (void*)(index_offset - HEADER_SIZE) );
        #endif
    }
    
    if( NULL != index->ps_locs ) {
        free_pseudo_locations( index->ps_locs );
    }

    free( index );

    return;
}

void
load_ondisk_index( char* index_fname, struct index_t** index )
{
    int rv;
    
    /* 
       first, 
       open the file containing the index to ensure that
       the magic number is correct and to get the size.
    */
    
    statmap_log( LOG_NOTICE, "Loading the index file '%s'",  index_fname  );

    FILE* index_fp = fopen( index_fname, "r" );
    if( index_fp == NULL )
    {
        statmap_log( LOG_FATAL, "Cannot open the index '%s' for reading",  index_fname  );
        assert( false );
        exit( 1 );
    }
    
    char magic_number;
    unsigned char indexed_seq_len;
    size_t index_size;

    rv = fread( &magic_number, 1, 1, index_fp );
    assert( rv == 1 );
    
    if( magic_number != 0 )
    {
        statmap_log( LOG_FATAL, "We appear to have loaded an invalid index" );
        exit( 1 );
    }
    
    /* Find the indexed sequence length */
    rv = fread( &indexed_seq_len, 1, 1, index_fp );
    assert( 1 == rv );
    assert( indexed_seq_len > 0 );

    /* allocate space for the index */
    *index = malloc( sizeof( struct index_t ) );
    
    /* set the index type */
    (*index)->index_type = TREE;

    /* Set the indexed sequence length */
    (*index)->seq_length = indexed_seq_len;

    rv = fread( &index_size, sizeof(size_t), 1, index_fp );
    assert( rv == 1 );
    
    fclose( index_fp );

    /* store the index */
    void* OD_index;
    
    #ifdef MMAP_INDEX
    int fd;
    if ((fd = open(index_fname, O_RDONLY )) < 0)
        statmap_log( LOG_FATAL, "can't create %s for reading",  index_fname );
    
    if ((OD_index = mmap (0, index_size + HEADER_SIZE, 
                             // index_size + sizeof(size_t) + sizeof(char), 
                          PROT_READ,
                          MAP_SHARED|MAP_POPULATE, fd, 0)) == (caddr_t) -1)
        statmap_log( LOG_FATAL, "mmap error for index file" );
    #else
    FILE* fp = NULL;
    fp = fopen( index_fname, "r"  );
    if( fp == NULL ) {
       statmap_log( LOG_FATAL, "Could not open index file at %s", index_fname );
    }
    
    OD_index = calloc( index_size + HEADER_SIZE, 1  );
    size_t res = fread( OD_index, 1, index_size + HEADER_SIZE, fp );
    statmap_log( LOG_DEBUG, "Read %zu bytes ( out of %zu + %zu )", res, HEADER_SIZE, index_size );
    fclose( fp );
    #endif
    
    index_offset = ((size_t) OD_index) + HEADER_SIZE; 
             // sizeof(size_t) + sizeof(char) + sizeof( unsigned char );
    
    /* the root of the tree is always at 0, before being offset 
       by index_offset */
    (*index)->index = 0;
    
    /* load the pseudo locations */
    char pseudo_loc_ofname[500];
    sprintf(pseudo_loc_ofname, "%s.pslocs", index_fname  );
    FILE* ps_fp = fopen( pseudo_loc_ofname, "r" );
    assert( NULL != ps_fp );
    load_pseudo_locations( ps_fp, &((*index)->ps_locs) );
    fclose( ps_fp );

    /* load the diploid mappings */
    char dmap_ofname[500];
    sprintf( dmap_ofname, "%s.dmap", index_fname );
    FILE* dmap_fp = fopen( dmap_ofname, "r" );
    assert( NULL != dmap_fp );
    load_diploid_maps_from_file( &((*index)->diploid_maps), dmap_fp );
    fclose( dmap_fp );
    
    return;
}

void
build_ondisk_index( struct index_t* index, char* ofname  )
{
    size_t total_freed_memory = 0;

    /* first, calculate the size of the tree in bytes */
    size_t index_size = sizeof_tree( index );

    unsigned char indexed_seq_len = index->seq_length;
    assert( index->seq_length < 255 );

    /* open/create the output file */
    int fdout;
    if ((fdout = open(ofname, O_RDWR | O_CREAT | O_TRUNC, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0)
    {
        statmap_log( LOG_FATAL, "can't create index file for writing" );
        exit( 1 );
    }
    
    /* go to the location corresponding to the last byte */
    if (lseek (fdout, HEADER_SIZE + index_size - 1, SEEK_SET) == -1)
    {
        statmap_log( LOG_FATAL, "lseek error" );
        exit( 1 );
    }
    
    /* write a dummy byte at the last location */
    if (write (fdout, "", 1) != 1)
    {
        statmap_log( LOG_FATAL, "write error" );
        exit( 1 );
    }
    
    /* mmap the output file */
    void* OD_index;
    if ((OD_index = mmap (0, index_size + HEADER_SIZE, 
                          PROT_READ | PROT_WRITE,
                          MAP_SHARED, fdout, 0)) == (caddr_t) -1)
    {
        statmap_log( LOG_FATAL, "mmap error for index output" );
        exit( 1 );
    }

    /* write the header information to the file */
    /* first, write a magic number. This allows us to tell if the file
       we just opened is a binary index or not */
    ((char*) OD_index)[0] = 0;
    OD_index += 1;
    
    /* next, write the indexed sequence length */
    ((unsigned char*) OD_index)[0] = indexed_seq_len;
    OD_index += 1;

    /* next, write the size of the index in bytes */
    /* 
     * this is a size_t on whatever architecture the index was created on. 
     * this is kind of a problem for portability, but there are much 
     * bigger problems if one wants to move binary indexes between 
     * architectures.
     *
     */
    ((size_t*) OD_index)[0] = index_size;
    OD_index += sizeof( size_t );

    /* Store the current position in mmapped file */
    void* curr_pos = OD_index;

    /* initialize a stack to store nodes that still need to be written */
    struct ODI_stack* stack;
    init_ODI_stack( &stack );
    
    /* initialize the current node to be the root node */
    void* root_node = index->index;
    add_ODI_stack_item( stack, &root_node, 's', 0 );
    
    /* while there are still nodes to add */
    unsigned int stack_index = 0;
    while( stack_index < stack->size )
    {
        /* index to loop through children */
        int i;

        /* write the node to disk */
        ODI_stack_item* curr_node = stack->stack + stack_index;
        size_t node_size = calc_node_size(
            curr_node->node_type, *(curr_node->node_ref));
        memcpy( curr_pos, *(curr_node->node_ref), node_size);
        /* free the memory */
        // memset( *(curr_node->node_ref), 0, node_size );
        free_node( *(curr_node->node_ref), curr_node->node_type);
        total_freed_memory += node_size;
        // fprintf( stderr, "==================================FREED %zu BYTES\n", total_freed_memory );    

        /* update the parent pointer */
        *(curr_node->node_ref) = curr_pos - (size_t)OD_index;
        
        /* DEBUG */
        /*
        fprintf( stdout, "NODE %i---- Level: %i\tType: %c\tSize: %zu\tPtr: %p\n",
                 stack_index, curr_node->level, curr_node->node_type, 
                 node_size, curr_pos );
        */
        
        int level = curr_node->level;
        
        /* add all of the children to the stack */
        switch( curr_node->node_type )
        {
        case 's':
            for( i = 0; i < 1<<(2*LETTER_LEN); i++ )
            {
                add_ODI_stack_item( 
                    stack, 
                    &(((static_node*)(curr_pos))[i].node_ref),
                    ((static_node*)curr_pos)[i].type,
                    level + 1
                );
            }
            break;
        case 'd':
            for( i = 0; 
                 i < get_dnode_num_children((dynamic_node*)curr_pos); 
                 i++ )
            {
                add_ODI_stack_item( 
                    stack, 
                    &(get_dnode_children((dynamic_node*)curr_pos)[i].node_ref),
                    get_dnode_children((dynamic_node*)curr_pos)[i].type,
                    level + 1
                );
            }
            break;
        default:
            break;
        }
        
        /* move to the next stack item */
        stack_index += 1;
        /* move the pointer to free space */
        curr_pos = ((char*) curr_pos) + node_size;
        // printf( "%p\t%i\t%i\n", curr_pos, stack_index, stack->size );
    }
    
    free_ODI_stack( stack );

    munmap( OD_index - HEADER_SIZE, 
            index_size + HEADER_SIZE );

    close( fdout );

    /* write the pseudo locations to file */
    char pseudo_loc_ofname[500];
    sprintf(pseudo_loc_ofname, "%s.pslocs", ofname  );
    FILE* ps_locs_of = fopen(pseudo_loc_ofname, "w");
    if( NULL == ps_locs_of ) {
        statmap_log( LOG_FATAL, "Error opening '%s' for writing.", PSEUDO_LOCATIONS_FNAME );
    }

    size_t size_written = 
        write_pseudo_locations_to_file( index->ps_locs, ps_locs_of );
    statmap_log( LOG_NOTICE, "wrote %zu bytes to pseudo locs file.", size_written );
    fclose(ps_locs_of);
    
    /* write diploid mappings to file */
    char diploid_mapping_ofname[500];
    sprintf( diploid_mapping_ofname, "%s.dmap", ofname );
    FILE* dmap_of = fopen( diploid_mapping_ofname, "w" );
    if( NULL == dmap_of ) {
        statmap_log( LOG_FATAL, "Error opening '%s' for writing", diploid_mapping_ofname );
    }

    write_diploid_maps_to_file( index->diploid_maps, dmap_of );
    fclose( dmap_of );

    return;
}


/****************************************************************************/

#ifdef PROFILE_MEMORY_USAGE
extern void 
print_memory_usage_stats()
{
    fprintf(stderr, "\nMemory Usage Stats:\n");

    fprintf(stderr, "\nStatic Nodes Size:               %i kB\n", 
            (sizeof(static_node)*num_static_nodes)/1024);
    fprintf(stderr, "Size of Static Node:             %i Bytes\n", 
            sizeof(static_node)); 
    fprintf(stderr, "static_children_bytes:           %i\n", 
            static_children_bytes);

    fprintf(stderr, "\nDynamic Nodes Size:              %i kB\n",  
            (sizeof(dynamic_node)*num_dynamic_nodes)/1024);
    fprintf(stderr, "Size of Dynamic Node:            %i Bytes\n", 
            sizeof(dynamic_node)); 
    fprintf(stderr, "Dynamic Children Bytes:          %i kB\n",  
            dynamic_children_bytes/1024);
    fprintf(stderr, "\nTotal Dynamic Nodes Size:        %i kB\n",  
            (sizeof(dynamic_node)*num_dynamic_nodes)/1024
            + dynamic_children_bytes/1024);

    fprintf(stderr, "\nSequence Nodes Size:             %i kB\n", 
            (sizeof(sequences_node)*num_sequence_nodes)/1024);
    fprintf(stderr, "Size of Sequence Node:           %i Bytes\n", 
            sizeof(sequences_node)); 
    fprintf(stderr, "Sequences Bytes:                 %i kB\n", 
            sequence_node_sequence_bytes/1024);
    fprintf(stderr, "Genome Location Bytes:           %i kB\n", 
            genome_location_bytes/1024);
    fprintf(stderr, "Dynamic Genome Location Bytes:   %i kB\n", 
            dynamic_genome_location_bytes/1024);
    fprintf(stderr, "\nTotal Sequences Nodes Size:      %i kB\n",  
            (sizeof(sequences_node)*num_sequence_nodes)/1024
            + sequence_node_sequence_bytes/1024);

    fprintf(stderr, "\n");

    return;
}
#endif
