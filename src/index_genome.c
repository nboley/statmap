/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <ctype.h>

#include "quality.h"
#include "index_genome.h"
#include "snp.h"
#include "pseudo_location.h"

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

inline float 
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

inline size_t
pmatch_stack_length( potential_match_stack* pmatch )
{
    return (pmatch->last_index - pmatch->first_index);
};

inline potential_match_stack*
add_pmatch( potential_match_stack* stack, 
            const void* node, const NODE_TYPE node_type, const enum STRAND strnd,
            const int node_level, const int max_num_levels,
            const float penalty, const float min_match_penalty, const float max_penalty_diff)
{
    /* TODO - lock access to the stack ( or move to berkeley db? ) */
    /* BUG!!!! change growth factor and catch memory corruption */

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

inline void init_static_node( static_node** node )
{
    /* note that the calloc initializes the children to NULL */
    /* FIXME check for NULL */
    
    *node = calloc( ALPHABET_LENGTH, sizeof(static_node_child) );

    #ifdef PROFILE_MEMORY_USAGE
    num_static_nodes++;
    #endif
}

inline void init_dynamic_node( dynamic_node** node )
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



extern void init_tree( index_t** index, int seq_length )
{   
    *index = malloc( sizeof( index_t ) );
    // assert( (*index)->index_type == TREE );
    static_node* tree_root;
    init_static_node( &tree_root );
    (*index)->index = tree_root;
    (*index)->index_type = TREE;
    (*index)->seq_length = seq_length;
}

inline dynamic_node*
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

inline void add_child_to_static_node( 
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

extern void
index_genome( struct genome_data* genome, int indexed_seq_len )
{
    /* initialize the tree structure */
    init_tree( &(genome->index), indexed_seq_len );

    /* TODO - use snps directly */
    struct snp_db_t* snp_db = genome->snp_db;
    /* 
     * If there is no initialized snp DB, then we set the
     * number of snps to zero so that we dont try and 
     * add any downstream. Maybe, in the future, I will 
     * just always the db but set the number of snps to 0.
     *
     */
    int num_snps;
    if( snp_db == NULL )
    {
        num_snps = 0;
    } else {
        num_snps = snp_db->num_snps;
    }
    
    /* initialize the constant loc elements */
    GENOME_LOC_TYPE loc;
    /* Not a junction read */
    loc.read_type = 0;
    
    int seq_len = genome->index->seq_length;
    char* tmp_seq = malloc(seq_len*sizeof(char));

    int chr_index;
    unsigned int bp_index;
    
    int snp_lb = 0;
    int snp_ub = 0;

    /* 
     * Iterate through each indexable sequence in the genome. If the 
     * sequence covers snps, then add them.
     */
    for( chr_index = 1; chr_index < genome->num_chrs; chr_index++ )
    {
        if( genome->chr_lens[chr_index] > LOCATION_MAX )
        {
            fprintf( stderr, 
                     "FATAL ERROR       :  Max indexable chr size is '%i'\n", 
                     LOCATION_MAX );
            exit( -1 );
        }

        fprintf(stderr, "NOTICE      :  Indexing '%s'\n", (genome->chr_names)[chr_index] );

        /* Set the chr index in the soon to be added location */
        loc.chr = chr_index;

        for( bp_index = 0; bp_index < genome->chr_lens[chr_index] - seq_len; bp_index += 1 )
        {            
            /* Set the basepair index */
            loc.loc = bp_index;
            
            /* 
             * Update the snp_lb. While the snps are strictly below bp_index, the
             * lower bound, then the snps can not fall into the sequence and there
             * is no need to consider these
             */
            while( snp_lb < num_snps
                   && snp_db->snps[snp_lb].loc.chr == chr_index 
                   && snp_db->snps[snp_lb].loc.loc < bp_index  )
            {
                snp_lb += 1;
            }
                        
            /* 
             * Update the snp_ub. While the snps are strictly below bp_index, the
             * lower bound, then the snps can not fall into the sequence and there
             * is no need to consider these
             */
            snp_ub = MAX( snp_lb, snp_ub );
            while( snp_ub < num_snps
                   && snp_db->snps[snp_ub].loc.chr == chr_index 
                   && snp_db->snps[snp_ub].loc.loc < bp_index + seq_len )
            {
                snp_ub += 1;
            }
                   
            /* 
             * If lb == ub, then we know this sequence doesnt cover any snps. 
             */
            if( snp_lb == snp_ub )
            {
                memcpy( tmp_seq, genome->chrs[chr_index] + bp_index, sizeof(char)*seq_len );
                
                /* Add the normal sequence */
                LETTER_TYPE *translation;
                translate_seq( tmp_seq, seq_len, &(translation));
                
                /* if we cant add this sequence ( probably an N ) continue */
                if( translation == NULL ) {
                    continue;
                }
                
                loc.covers_snp = 0;
                loc.snp_coverage = 0;
                
                /* Add the sequence into the tree */
                add_sequence(genome->index, genome->ps_locs, translation, seq_len, loc);
                
                free( translation );                                
            } else {
                /* If there are too many snps in this sequence, print a warning */
                if( snp_ub - snp_lb > MAX_NUM_SNPS )
                {
                    fprintf(stderr, "ERROR       :  Can not add the sequence at 'chr %i: bp %i' - it contains %i SNPs ( %i max )\n", 
                            chr_index, bp_index, snp_ub - snp_lb, MAX_NUM_SNPS );
                    continue;
                }
                
                /* 
                 * We need to iterate through every combination of snps in the following
                 * list. Ie, if there are 2 snps, then we need to add all of the 
                 * subsequences with s1 on, s2 on, sn1 on, s2 off, s1 off, s2 off, 
                 * and both off.  This might be hard in general, but there
                 * is an easy way of dealing with this. Since there are 2^NS -1 total
                 * combinations, we just count from 1 to 2^NS. Then, if the ith bit
                 * is set, we know that the ith snp should be on. 
                 */
                
                /* Make sure there is room in the bitmap to store all of the snips */
                assert( sizeof(unsigned int)*8 > MAX_NUM_SNPS  );
                /* bm stands for bitmap */
                unsigned int bm;
                assert( snp_ub > snp_lb );
                for( bm = 0; bm < (((unsigned int)1) << (snp_ub - snp_lb)); bm++ )
                {
                    /* Make a copy of the sequence that we will be mutating */
                    /* 
                     * TODO - Make this more efficient. We shouldnt need to recopy the sequence
                     * every single time. Although, since snps should typically be pretty 
                     * sparse, this may not matter in practice.
                     */
                    memcpy( tmp_seq, genome->chrs[chr_index] + bp_index, sizeof(char)*seq_len );
                    
                    /* 
                     * Loop through each possible snp. If the bit is set, then 
                     * set the bp location in the seq to the alternate.
                     */
                    int snp_index; 
                    for( snp_index = 0; snp_index < snp_ub - snp_lb; snp_index++ )
                    {
                        /* If the correct bit is set */
                        if( (bm&(1<<snp_index)) > 0 )
                            tmp_seq[ snp_db->snps[snp_lb+snp_index].loc.loc - bp_index  ] 
                                = snp_db->snps[snp_lb+snp_index].alt_bp;                        
                    }
                    
                    LETTER_TYPE *translation;
                    translate_seq( tmp_seq, seq_len, &(translation));
                    
                    /* if we cant add this sequence ( probably an N ) continue */
                    if( translation == NULL ) {
                        continue;
                    }
                    
                    loc.covers_snp = 1;
                    loc.snp_coverage = bm;
                    
                    /* Add the sequence into the tree */
                    add_sequence(genome->index, genome->ps_locs, 
                                 translation, seq_len, loc);
                    
                    free( translation );                
                }
            }
        }
    }
    free( tmp_seq );

    /* sort all of the pseudo locations */
    sort_pseudo_locations( genome->ps_locs );


    return;
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
            perror("Impossible branch reached in free_node.\n");
            abort();
    }
}

void free_node_and_children( void* node, char node_type )
{
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
            fprintf(stderr, "FATAL       :  Impossible Branch: %c\n", node_type);
            perror("Impossbile branch reached in free_node_and_children.\n");
            abort();
    }
}

void free_tree( index_t* root )
{
    assert( root->index_type == TREE );
    free_node_and_children( (static_node*) root->index, 's' );
    free( root );
}

/**** END CLEANUP functions ***************************************************/

inline void 
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

inline void 
build_dynamic_node_from_sequence_node(  sequences_node* qnode, 
                                        dynamic_node** dnode,
                                        const int num_levels,
                                        LEVEL_TYPE level   )
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
                GENOME_LOC_TYPE* gen_locs 
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
                GENOME_LOC_TYPE* gen_locs 
                    = get_overflow_genome_locations_array_start( qnode,  num_letters )
                    + loc.locs_array.locs_start;

                int j;
                for(j = 0; j < loc.locs_array.locs_size; j++ )
                {
                    GENOME_LOC_TYPE loc = gen_locs[j];

                    *child_seqs = add_sequence_to_sequences_node(   
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

inline int 
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


inline void 
add_sequence( index_t* index, struct pseudo_locations_t* ps_locs,
              LETTER_TYPE* seq, const int seq_length, 
              GENOME_LOC_TYPE genome_loc ) 
{
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
            level 
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


inline void 
find_matches( void* node, NODE_TYPE node_type, int node_level, 
              const int seq_length,
              float curr_penalty, 
              float min_match_penalty,
              /* 
               * the maximum spread between the highest penalty 
               * match and any other valid match 
               */
              float max_penalty_spread,
              mapped_locations* results,

              /* fwd stranded data  */
              LETTER_TYPE* seq_1, 
              const float* const position_mutation_prs_1,
              const float* const inverse_position_mutation_prs_1,

              /* rev stranded data */
              LETTER_TYPE* seq_2, 
              const float* const position_mutation_prs_2,
              const float* const inverse_position_mutation_prs_2,

              const float* const lookuptable_bp
    )

{
    const int num_letters = calc_num_letters( seq_length );

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


    while( pmatch_stack_length( stack ) > 0 )
    {
        potential_match match = pop_pmatch( stack );
        
        void* node = match.node;
        NODE_TYPE node_type = match.node_type; 
        int node_level = match.node_level;
        float curr_penalty = match.penalty;
        enum STRAND strnd = match.strnd;
        
        LETTER_TYPE* seq;
        const float* position_mutation_prs;
        const float* inverse_position_mutation_prs;
 
        if( strnd == FWD )
        {
            seq = seq_1;
            position_mutation_prs = position_mutation_prs_1;
            inverse_position_mutation_prs = inverse_position_mutation_prs_1;
        } else {
            seq = seq_2;
            position_mutation_prs = position_mutation_prs_2;
            inverse_position_mutation_prs = inverse_position_mutation_prs_2;
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
                float penalty_addition = penalty_func(
                    letter, seq[node_level], 
                    node_level, seq_length,
                    min_match_penalty - curr_penalty,
                    position_mutation_prs,
                    inverse_position_mutation_prs,
                    lookuptable_bp
                );

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

                float penalty_addition = penalty_func(
                    letter, seq[node_level], 
                    node_level, seq_length,
                    min_match_penalty - curr_penalty,
                    position_mutation_prs,
                    inverse_position_mutation_prs,
                    lookuptable_bp
                );

                /* 
                 * If penalty_func returns a value >= 1
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
                strnd
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

            /* debugging
            printf("Level: %i\t Node Type: %c\t\t\tPenalty: %e\n", 
                   node_level, node_type, curr_penalty );
            */
            float max_added_penalty = 
            find_sequences_in_sequences_node( 
                node, 
                curr_penalty, min_match_penalty,
                seq, seq_length, num_letters, node_level, 
                strnd, results,
                position_mutation_prs,
                inverse_position_mutation_prs,
                lookuptable_bp
            );

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

    free_pmatch_stack( stack );
    return;
}


extern void
find_matches_from_root( index_t* index,
             
                        float min_match_penalty,
                        float max_penalty_spread,
                        mapped_locations* results,

                        /* the length of the two reads ( below ) */
                        const int read_len,

                        /* the fwd stranded data */
                        LETTER_TYPE* seq_1,
                        float* lookuptable_position_1,
                        float* inverse_lookuptable_position_1,

                        /* the bkwd stranded data */
                        LETTER_TYPE* seq_2,
                        float* lookuptable_position_2,
                        float* inverse_lookuptable_position_2,

                        float* lookuptable_bp
)
{
    assert( index->index_type == TREE );
    /* FIXME - allow for different read lens and index lens */
    assert( index->seq_length == read_len );
    static_node* root = index->index;

    return find_matches( (void*) root, 's', 0, 
                         index->seq_length,
                         
                         0, min_match_penalty, 
                         max_penalty_spread,
                         results,

                         seq_1,
                         lookuptable_position_1,
                         inverse_lookuptable_position_1,

                         seq_2,
                         lookuptable_position_2,
                         inverse_lookuptable_position_2,

                         lookuptable_bp
    ); 
}

void
search_index( index_t* index, 
              
              float min_match_penalty,
              float max_penalty_spread,
              mapped_locations* results,

              struct rawread* r,
              float* bp_mut_rates
    )
{
    /**** Prepare the read for the index search */
    /* Build the quality lookup tables */
    float* lookuptable_position = malloc(sizeof(float)*r->length);
    float* inverse_lookuptable_position = malloc(sizeof(float)*r->length);
    float* reverse_lookuptable_position = malloc(sizeof(float)*r->length);
    float* reverse_inverse_lookuptable_position = malloc(sizeof(float)*r->length);
    build_lookup_table_from_rawread(
        r, 
        lookuptable_position, inverse_lookuptable_position,
        reverse_lookuptable_position, reverse_inverse_lookuptable_position
    );
    
    /* Store a copy of the read */
    /* This read has N's replaced with A's, and might be RC'd */
    char* tmp_read = calloc(r->length + 1, sizeof(char));
    assert( tmp_read != NULL );
    memcpy( tmp_read, r->char_seq, sizeof(char)*(r->length) );
    /* note that the NULL ending is pre-set from the calloc */
    replace_ns_inplace( tmp_read, r->length );

    /** Deal with the read on the fwd strand */
    /* Store the translated sequences here */
    LETTER_TYPE *fwd_seq;
    fwd_seq = translate_seq( r->char_seq, r->length, &fwd_seq );
    /* If we couldnt translate it */
    if( fwd_seq == NULL )
    {
        // fprintf(stderr, "Could Not Translate: %s\n", r->char_seq);
        goto cleanup_table;
    }
    assert( fwd_seq != NULL );
    
    /** Deal with the read on the opposite strand */
    LETTER_TYPE *bkwd_seq;
    rev_complement_read( r->char_seq, tmp_read, r->length );
    bkwd_seq = translate_seq( tmp_read, r->length, &bkwd_seq );
    replace_ns_inplace( tmp_read, r->length );
    assert( bkwd_seq != NULL );
    
    /* map the full read */
    find_matches_from_root( index, 
                            
                            min_match_penalty,
                            max_penalty_spread,
                            results,

                            /* length of the reads */
                            r->length,
                            
                            /* the fwd stranded sequence */
                            fwd_seq, 
                            lookuptable_position,
                            inverse_lookuptable_position,
                            
                            /* the bkwd stranded sequence */
                            bkwd_seq, 
                            reverse_lookuptable_position,
                            reverse_inverse_lookuptable_position,
                            
                            bp_mut_rates
        );
    
    /* Free the allocated memory */
    free( fwd_seq );
    free( bkwd_seq );
    free( tmp_read );

cleanup_table:

    free( lookuptable_position );
    free( inverse_lookuptable_position );
    free( reverse_lookuptable_position );
    free( reverse_inverse_lookuptable_position );

    return;
};


size_t
size_of_snode( )
{
    return sizeof(static_node_child)*(2<<(2*LETTER_LEN));

}

size_t
calc_node_size( NODE_TYPE type, void* node  )
{
    size_t size = 0;
    int i;

    if( node == NULL )
        return 0;

    switch ( type ) 
    {
    case 'q':
        /* sequence nodes dont have any children */
        return get_num_used_bytes( node );
        break;
    case 'l':
        /* locations nodes dont have any children */
        return size_of_locations_node( node );
        break;
    case 's':
        for( i = 0; i < 1<<(2*LETTER_LEN); i++ )
        {
            size += calc_node_size( 
                ((static_node*)node)[i].type, 
                ((static_node*)node)[i].node_ref 
            );
        }
        return size + size_of_snode();
        break;
    case 'd':
        for( i = 0; i < get_dnode_num_children( node ); i++ )
        {
            size += calc_node_size( 
                get_dnode_children(node)[i].type, 
                get_dnode_children(node)[i].node_ref 
            );
        }
        return size + size_of_dnode( node );
        break;
    }

    fprintf(stderr, "FATAL       :  Unrecognized Node Type: '%c'\n", type);
    assert( false );
    exit( -1 );
}

size_t 
sizeof_tree( index_t* index )
{
    assert( index->index_type == TREE );
    return calc_node_size( 's', index->index );
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
typedef struct {
    size_t allocated_size;
    size_t size;
    ODI_stack_item* stack;
} ODI_stack;

#define ODI_stack_GROW_FACTOR 100000

void
init_ODI_stack( ODI_stack** stack )
{
    *stack = malloc( sizeof( ODI_stack ) );
    (*stack)->allocated_size = ODI_stack_GROW_FACTOR;
    (*stack)->size = 0;
    (*stack)->stack = 
        malloc( sizeof( ODI_stack_item )*ODI_stack_GROW_FACTOR );
    
    return;
}

void
add_ODI_stack_item( ODI_stack* stack, 
                    void** node_ref, 
                    NODE_TYPE node_type, 
                    LEVEL_TYPE level )
{
    stack->size++;
    
    if( stack->size == stack->allocated_size )
    {
        stack->allocated_size += ODI_stack_GROW_FACTOR;
        stack->stack = realloc( 
            stack->stack, sizeof(ODI_stack_item)*(stack->allocated_size) );
        if( stack->stack == NULL )
        {
            fprintf( stderr, "FATAL       :  Error allocating memory for the ODI Stack.\n" );
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
cmp_ODI_stack_item( const void* a, const void* b )
{
    int rv = 0;
    rv = ((ODI_stack_item*)b)->level - ((ODI_stack_item*)a)->level;
    return rv;
} 

void
sort_ODI_stack( ODI_stack* stack )
{
    qsort ( 
        stack->stack,
        stack->size, 
        sizeof(ODI_stack_item), 
        cmp_ODI_stack_item
    );
    
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




