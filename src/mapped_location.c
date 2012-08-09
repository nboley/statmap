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
    
    if( loc1->location.chr != loc2->location.chr )
        return loc1->location.chr - loc2->location.chr;
    
    return loc1->location.loc - loc2->location.loc;
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
                         INDEX_LOC_TYPE location, 
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
    mapped_location* loc = results->locations + results->length;

    /* set the location */
    loc->location = location;
    assert( loc->location.loc >= 0 );
    
    /* set the read strand */
    loc->strnd = strnd;

    /* set the penalty */
    loc->penalty = penalty;
    
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
                             loc->location,
                             loc->strnd,
                             loc->penalty
        );
}

void
expand_diploid_mapped_location(
        mapped_location* loc,
        mapped_locations* locs,
        struct genome_data* genome
    )
{
    assert( loc->location.is_paternal && loc->location.is_maternal );
    /* This should only be called on real locations */
    assert( loc->location.chr != PSEUDO_LOC_CHR_INDEX );

    /* build the paternal location */
    mapped_location paternal;
    copy_mapped_location( &paternal, loc );
    
    /* Set the diploid flags */
    paternal.location.is_paternal = 1;
    paternal.location.is_maternal = 0;

    /* Since the shared diploid location uses the chr and loc from the paternal
     * copy, we don't need to change anything else */
    add_mapped_location( &paternal, locs );

    /* build the maternal location */
    mapped_location maternal;
    copy_mapped_location( &maternal, loc );
    maternal.location.is_paternal = 0;
    maternal.location.is_maternal = 1;
    /* lookup the maternal location from the paternal location information
     * used on the shared diploid location */

    int paternal_chr = loc->location.chr;
    int paternal_loc = loc->location.loc;

    int maternal_chr = -1;
    int maternal_loc = -1;
    build_maternal_loc_from_paternal_loc(
            &maternal_chr, &maternal_loc,
            paternal_chr, paternal_loc,
            genome
        );

    /* finished building the maternal mapped_location */
    maternal.location.chr = maternal_chr;
    maternal.location.loc = maternal_loc;
    add_mapped_location( &maternal, locs );

    return;
}

void
add_and_expand_mapped_location(
    mapped_location* loc,
    mapped_locations* locs,
    struct genome_data* genome
)
{
    /* if this is a shared diploid location, expand it into the paternal and
     * maternal copies */
    if( loc->location.is_paternal && loc->location.is_maternal )
    {
        assert( genome != NULL );
        expand_diploid_mapped_location( loc, locs, genome );
    } else {
        add_mapped_location( loc, locs );
    }
}

void
copy_mapped_location( mapped_location* dest, mapped_location* src ) 
{
    *dest = *src;
    return;
}

void
print_mapped_locations( mapped_locations* results )
{
    int i;
    // printf("Num:\tPenalty\tLoc\n");
    
    for( i = 0; i < results->length; i++)
    {
        printf("\t%i:%i\t%.6f\n", 
               results->locations[i].location.chr,
               results->locations[i].location.loc,
               results->locations[i].penalty
        );
    }
    
    return;
}

/*
 *  END Mapped Location
 *
 **************************************************************************/


/*** Mapped locations container ***/
void
init_mapped_locations_container(
        mapped_locations_container** mlc
    )
{
    *mlc = malloc( sizeof( mapped_locations_container ) );

    (*mlc)->container = NULL;
    (*mlc)->length = 0;
}

void
free_mapped_locations_container(
        mapped_locations_container* mlc
    )
{
    if( mlc == NULL ) return;

    // free the mapped locations stored in this container
    int i;
    for( i = 0; i < mlc->length; i++ )
    {
        free_mapped_locations( mlc->container[i] );
    }

    free( mlc->container );

    free( mlc );
}

void
add_mapped_locations_to_mapped_locations_container(
        mapped_locations* ml,
        mapped_locations_container* mlc
    )
{
    mlc->length += 1;
    mlc->container = realloc( mlc->container,
            sizeof( mapped_locations* ) * mlc->length );

    mlc->container[ mlc->length-1 ] = ml;
}
