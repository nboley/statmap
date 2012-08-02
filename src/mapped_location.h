/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef MAPPED_LOCATION
#define MAPPED_LOCATION

#include <stdio.h>

#include "config.h"
#include "genome.h"
#include "rawread.h"
#include "candidate_mapping.h"

/* block size for newly allocated results memory */
#define RESULTS_GROWTH_FACTOR 100



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

typedef struct {
    GENOME_LOC_TYPE location;
    enum STRAND strnd;
    int trim_offset;
    float penalty;
} mapped_location;

typedef struct {
    /* An array of mapped locations */
    mapped_location* locations;
    int length;
    int allocated_length;

    int subseq_len;
    int subseq_offset;

    enum bool skip;
} mapped_locations;

int
cmp_mapped_locations_by_location( void* loc1, void* loc2 );

int
cmp_mapped_locations_by_penalty( void* loc1, void* loc2 );

void 
sort_mapped_locations_by_location( mapped_locations* results );

void 
sort_mapped_locations_by_penalty( mapped_locations* results );

/* Deal with mapped locations arrays */

void
init_mapped_locations( mapped_locations** results );

void
free_mapped_locations( mapped_locations* results );

void
add_mapped_location( mapped_locations* results, 
                     GENOME_LOC_TYPE location, 
                     enum STRAND strnd,
                     int trim_offset,
                     float penalty );

void
copy_mapped_location(
        mapped_location* loc,
        mapped_locations* locs
    );

void
print_mapped_locations( mapped_locations* results );

/*
 *  END Mapped Location
 *
 **************************************************************************/

/*** mapped_locations_container ***/

typedef struct {
    /* container of pointers to mapped_locations */
    mapped_locations** container;
    int length;
} mapped_locations_container;

void
init_mapped_locations_container(
        mapped_locations_container** mlc
    );

void
free_mapped_locations_container(
        mapped_locations_container* mlc
    );

void
add_mapped_locations_to_mapped_locations_container(
        mapped_locations* ml,
        mapped_locations_container* mlc
    );

#endif /* define MAPPED_LOCATION */
