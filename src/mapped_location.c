/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <assert.h>
#include <string.h> // memcpy()

#include "mapped_location.h"
#include "diploid_map_data.h"

#include "log.h"

/*************************************************************************
 *
 *  Mapped Location 
 *
 *  Mapped locations are the data structures returned by an index 
 *  lookup. They just store the seq penalty, and the data type that
 *  is stoed inside of the index, GEN_LOC_TYPE.
 *
 *
 */

int
cmp_mapped_locations_by_location( const mapped_location* loc1, 
                                  const mapped_location* loc2 )
{
    if( loc1->strnd != loc2->strnd )
        return loc1->strnd - loc2->strnd;
    
    if( loc1->chr != loc2->chr )
        return loc1->chr - loc2->chr;
    
    return loc1->loc - loc2->loc;
}

int
cmp_mapped_locations_by_penalty( const mapped_location* loc1, 
                                 const mapped_location* loc2 )
{
    if( loc1->penalty > loc2->penalty )
        return -1;

    if( loc1->penalty < loc2->penalty )
        return 1;
    
    return 0;

}

void 
sort_mapped_locations_by_location( mapped_locations* results )
{
  qsort( results->locations, 
         results->length, 
         sizeof(mapped_location),
         (int(*)(const void*, const void*))cmp_mapped_locations_by_location
  );
}


void 
sort_mapped_locations_by_penalty( mapped_locations* results )
{
  qsort( results->locations, 
         results->length, 
         sizeof(mapped_location),
         (int(*)(const void*, const void*))cmp_mapped_locations_by_penalty
  );
}

void
init_mapped_locations(
        mapped_locations** results,
        struct indexable_subtemplate* probe
    )
{
    *results = malloc( sizeof(mapped_locations) );
    (*results)->locations = 
        malloc(RESULTS_GROWTH_FACTOR*sizeof(mapped_location));
    
    (*results)->length = 0;
    (*results)->allocated_length = RESULTS_GROWTH_FACTOR;

    (*results)->probe = probe;
    
    return;
}


void
free_mapped_locations( mapped_locations* results )
{
    if( results == NULL )
        return;
    
    free( results->locations  );
    //free_indexable_subtemplates( results->probe );
    free( results  );
    return;
}

void
add_new_mapped_location( mapped_locations* results, 
                         unsigned int chr,
                         unsigned int loc,
                         enum STRAND strnd,
                         float penalty )
{
    // use 0.1 to avoid rounding error false asserts 
    assert( penalty < 0.1 );

    /* 
     * test to see if there is enough allocated memory in results
     * if there isn't then realloc
     */
    if( results->length == results->allocated_length )
    {
        results->allocated_length += RESULTS_GROWTH_FACTOR;
        results->locations = realloc(
            results->locations,
            results->allocated_length*sizeof(mapped_location)
        );
        
        if( results->locations == NULL )
        {
            exit(1);
        }
    }

    /* This should be optimized out */
    mapped_location* new_loc = results->locations + results->length;

    /* set the location */
    new_loc->chr = chr;
    new_loc->loc = loc;
    // assert( new_loc->loc >= 0 ); // always true
    
    /* set the read strand */
    new_loc->strnd = strnd;

    /* set the penalty */
    new_loc->penalty = penalty;

    /* add the new results to the end of the results set */
    results->length++;

    return;
}

/* Wrapper to copy an existing mapped_location using add_mapped_location */
void
add_mapped_location(
    mapped_location* loc,
    mapped_locations* locs
)
{
    add_new_mapped_location( locs,
                             loc->chr,
                             loc->loc,
                             loc->strnd,
                             loc->penalty );
}

void
copy_mapped_location( mapped_location* dest, mapped_location* src ) 
{
    *dest = *src;
    return;
}

void
expand_diploid_index_location(
        INDEX_LOC_TYPE* iloc,
        mapped_locations* results,
        enum STRAND strnd,
        float penalty,
        struct genome_data* genome
    )
{
    /* This should only be called on a shared diploid location */
    assert( iloc->is_paternal && iloc->is_maternal );
    /* This should only be called on real locations */
    // TODO ? really? or should we just skip them for now?
    assert( iloc->chr != PSEUDO_LOC_CHR_INDEX );

    /* add the paternal location.
     * the shared diploid location uses the paternal chr and loc, so we don't
     * need to do anything extra here */
    add_new_mapped_location( results,
                             iloc->chr,
                             iloc->loc,
                             strnd,
                             penalty );

    /* build the maternal location */
    /* lookup the maternal location from the paternal location information
     * stored on the shared diploid location */
    int paternal_chr = iloc->chr;
    int paternal_loc = iloc->loc;

    int maternal_chr = -1;
    int maternal_loc = -1;
    build_maternal_loc_from_paternal_loc(
            &maternal_chr, &maternal_loc,
            paternal_chr, paternal_loc,
            genome
        );
    assert( maternal_chr != -1 );
    assert( maternal_loc != -1 );

    add_new_mapped_location( results,
                             maternal_chr,
                             maternal_loc,
                             strnd,
                             penalty );
    return;
}

void
add_and_expand_location_from_index(
        mapped_locations* results,
        INDEX_LOC_TYPE* iloc,
        enum STRAND strnd,
        float penalty,
        struct genome_data* genome
    )
{
    if( iloc->is_paternal && iloc->is_maternal )
    {
        expand_diploid_index_location( iloc, results, strnd, penalty, genome );
    } else {
        add_new_mapped_location( results,
                                 iloc->chr,
                                 iloc->loc,
                                 strnd,
                                 penalty );
    }
}

void
print_mapped_locations( mapped_locations* results )
{
    int i;
    
    for( i = 0; i < results->length; i++)
    {
        printf("\t%i:%i\t%d\t%.6f\n", 
               results->locations[i].chr,
               results->locations[i].loc,
               results->locations[i].strnd,
               results->locations[i].penalty
        );
    }
    
    return;
}

void
sort_search_results(mapped_locations** search_results)        
{
    int i;
    for(i = 0; search_results[i] != NULL; i++ )
    {
        sort_mapped_locations_by_location( search_results[i] );
    }

    return;
}

void
free_search_results( mapped_locations** search_results )
{
    /* free each of the mapped_locations stored in the array */
    int i;
    for( i = 0; search_results[i] != NULL; i++ )
    {
        free_mapped_locations( search_results[i] );
    }

    /* free the array of pointers */
    free( search_results );

    return;
}

/*
 *  END Mapped Location
 *
 **************************************************************************/


/***** ml_match *****/

void
init_ml_match( struct ml_match** match, int match_len )
{
    *match = malloc( sizeof( struct ml_match ));

    /* Note - the length of the arrays will always be equal to the number of
     * indexable subtemplates, since we must be able to match across all of the
     * indexable subtemplates for a valid match. */
    (*match)->len = match_len;
    (*match)->matched = 0;

    (*match)->locations = calloc( match_len, sizeof( mapped_location ));
    (*match)->subseq_lengths = calloc( match_len, sizeof( int ));
    (*match)->subseq_offsets = calloc( match_len, sizeof( int ));

    (*match)->cum_penalty = 0;

    return;
}

struct ml_match*
copy_ml_match( struct ml_match* match )
{
    /* Return a copy of match */
    struct ml_match* match_copy;
    init_ml_match( &match_copy, match->len );
    match_copy->matched = match->matched;

    /* copy the arrays */
    memcpy( match_copy->locations, match->locations,
            sizeof( mapped_location )*match->len );
    memcpy( match_copy->subseq_lengths, match->subseq_lengths,
            sizeof( int )*match->len );
    memcpy( match_copy->subseq_offsets, match->subseq_offsets,
            sizeof( int )*match->len );

    /* copy the cumulative penalty */
    match_copy->cum_penalty = match->cum_penalty;

    return match_copy;
}

void
free_ml_match( struct ml_match* match )
{
    if( match == NULL )
        return;

    free( match->locations );
    free( match->subseq_lengths );
    free( match->subseq_offsets );

    free( match );

    return;
}

void
add_location_to_ml_match(
        mapped_location* location,
        struct ml_match* match, 
        int subseq_length,
        int subseq_offset )
{
    assert( match->matched < match->len );

    match->locations[match->matched] = *location;
    match->subseq_lengths[match->matched] = subseq_length;
    match->subseq_offsets[match->matched] = subseq_offset;

    match->cum_penalty += location->penalty;

    match->matched++;

    return;
}

/***** ml_matches *****/

void
init_ml_matches( struct ml_matches** matches )
{
    *matches = malloc( sizeof( struct ml_matches ));

    (*matches)->matches =
        malloc( ML_MATCHES_GROWTH_FACTOR*sizeof(struct ml_match *) );
    (*matches)->length = 0;
    (*matches)->allocated_length = ML_MATCHES_GROWTH_FACTOR;

    return;
}

void
free_ml_matches( struct ml_matches* matches )
{
    /* Free the individual ml_matches */
    int i;
    for( i = 0; i < matches->length; i++ )
    {
        free_ml_match( matches->matches[i] );
    }

    free( matches->matches );
    free( matches );
    return;
}

void
copy_ml_match_into_matches(
        struct ml_match* match,
        struct ml_matches* matches )
{
    /* 
     * test to see if there is enough allocated memory in results
     * if there isn't then realloc
     */
    if( matches->length == matches->allocated_length )
    {
        matches->allocated_length *= ML_MATCHES_GROWTH_FACTOR;
        matches->matches = realloc(
                matches->matches,
                matches->allocated_length*sizeof(struct ml_match*)
            );

        if( matches->matches == NULL )
        {
            statmap_log( LOG_FATAL, "Failed realloc in add_ml_match" );
            assert(false);
            exit(1);
        }
    }

    /* copy the ml_match */
    (matches->matches)[matches->length] = copy_ml_match( match );

    /* Update the number of matches */
    matches->length++;

    return;
}

/***** ml_match_stack *****/

void
init_ml_match_stack(
        struct ml_match_stack** stack )
{
    *stack = malloc( sizeof( struct ml_match_stack ));
    (*stack)->top = -1; // empty
    return;
}

void
free_ml_match_stack(
        struct ml_match_stack* stack )
{
    free( stack );
    return;
}

enum bool
ml_match_stack_is_empty(
        struct ml_match_stack* stack )
{
    if( stack->top == -1 )
        return true;

    return false;
}

void
ml_match_stack_push(
        struct ml_match_stack* stack,
        struct ml_match* match )
{
    /* Make sure we won't exceed the maximum size of the stack */
    assert( (stack->top + 1) < MAX_ML_MATCH_STACK_LEN );

    stack->top += 1;
    stack->stack[stack->top] = match;

    return;
}

struct ml_match*
ml_match_stack_pop(
        struct ml_match_stack* stack )
{
    /* Make sure we don't try to pop off an empty stack */
    assert( !ml_match_stack_is_empty( stack ));

    stack->top -= 1;
    return stack->stack[stack->top + 1];
}
