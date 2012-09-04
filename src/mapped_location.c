/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <assert.h>

#include "mapped_location.h"
#include "diploid_map_data.h"

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
            fprintf(stderr, "Failed realloc in add_mapped_locations\n");
            exit(1);
        }
    }

    /* This should be optimized out */
    mapped_location* new_loc = results->locations + results->length;

    /* set the location */
    new_loc->chr = chr;
    new_loc->loc = loc;
    assert( new_loc->loc >= 0 );
    
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
        assert( genome != NULL );
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
    // printf("Num:\tPenalty\tLoc\n");
    
    for( i = 0; i < results->length; i++)
    {
        printf("\t%i:%i\t%.6f\n", 
               results->locations[i].chr,
               results->locations[i].loc,
               results->locations[i].penalty
        );
    }
    
    return;
}

/*
 *  END Mapped Location
 *
 **************************************************************************/

void
init_matched_mapped_locations(
        struct matched_mapped_locations **m,
        mapped_location* base,
        struct indexable_subtemplate* base_probe,
        int num_match_containers )
{
    (*m) = malloc( sizeof( struct matched_mapped_locations ));

    (*m)->base = base;
    (*m)->base_probe = base_probe;

    (*m)->num_match_containers = num_match_containers;
    /* We calloc here so every pointer is initialized to NULL. We will leave the
     * pointer for the base's indexable_subtemplate NULL */
    (*m)->match_containers = calloc( num_match_containers,
            sizeof( mapped_locations* ));

    return;
}

void
free_matched_mapped_locations(
        struct matched_mapped_locations *m )
{
    /* free the mapped_locations */
    int i;
    for( i = 0; i < m->num_match_containers; i++ )
    {
        free_mapped_locations( m->match_containers[i] );
    }

    /* free the array of pointers to mapped locations */
    free( m->match_containers );

    /* free the whole structure */
    free( m );

    return;
}

void
add_matches_to_matched_mapped_locations(
        mapped_locations* matching_subset,
        int originating_ist_index,
        struct matched_mapped_locations* all_matches )
{
    all_matches->match_containers[originating_ist_index] = matching_subset;
    return;
}
