/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <assert.h>

#include "mapped_location.h"

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
cmp_mapped_locations_by_location( void* loc1, void* loc2 )
{
    /* first test the chromosome identifier */
    if( ((mapped_location*) loc1)->location.chr
        > ((mapped_location*) loc2)->location.chr )
    {
        return 1;
    }
  
    if( ((mapped_location*) loc1)->location.chr 
        < ((mapped_location*) loc2)->location.chr )
    {
        return -1;
    }

    /* since the chromosomes must be identical... */
    if( ((mapped_location*) loc1)->location.loc
        > ((mapped_location*) loc2)->location.loc )
    {
        return 1;
    }
  
    if( ((mapped_location*) loc1)->location.loc 
        < ((mapped_location*) loc2)->location.loc )
    {
        return -1;
    }

    
    return 0;
}

int
cmp_mapped_locations_by_penalty( void* loc1, void* loc2 )
{
    if( ((mapped_location*) loc1)->penalty
        > ((mapped_location*) loc2)->penalty )
    {
        return -1;
    } 
 
    if( ((mapped_location*) loc1)->penalty 
        < ((mapped_location*) loc2)->penalty )
    {
        return 1;
    }
    
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
init_mapped_locations( mapped_locations** results )
{
    *results = malloc( sizeof(mapped_locations) );
    (*results)->locations = 
        malloc(RESULTS_GROWTH_FACTOR*sizeof(mapped_location));
    
    (*results)->length = 0;
    (*results)->allocated_length = RESULTS_GROWTH_FACTOR;
    
    (*results)->subseq_len = -1;
    (*results)->subseq_offset = -1;

    (*results)->skip = false;
    
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
add_mapped_location( mapped_locations* results, 
                     GENOME_LOC_TYPE location, 
                     enum STRAND strnd,
                     int trim_offset,
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

    /* set the trim offset */
    loc->trim_offset = trim_offset;

    /* set the penalty */
    loc->penalty = penalty;
    
    /* add the new results to the end of the results set */
    results->length++;

    return;
}

/* Wrapper to copy an existing mapped_location using add_mapped_location */
void
copy_mapped_location(
        mapped_location* loc,
        mapped_locations* locs
    )
{
    add_mapped_location( locs,
            loc->location,
            loc->strnd,
            loc->trim_offset,
            loc->penalty
        );
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
