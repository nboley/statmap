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
    unsigned int chr; // TODO use MRL_CHR_TYPE, etc. ?
    unsigned int loc; // maybe put some unified types in config.h?
    enum STRAND strnd;
    float penalty;
} mapped_location;

typedef struct {
    /* An array of mapped locations */
    mapped_location* locations;
    int length;
    int allocated_length;

    struct indexable_subtemplate* probe;
} mapped_locations;

int
cmp_mapped_locations_by_location( const mapped_location* loc1, 
                                  const mapped_location* loc2 );

int
cmp_mapped_locations_by_penalty( const mapped_location* loc1, 
                                 const mapped_location* loc2 );

void 
sort_mapped_locations_by_location( mapped_locations* results );

void 
sort_mapped_locations_by_penalty( mapped_locations* results );

/* Deal with mapped locations arrays */

void
init_mapped_locations(
        mapped_locations** results,
        struct indexable_subtemplate* origin
    );

void
free_mapped_locations( mapped_locations* results );

void
add_new_mapped_location(
    mapped_locations* results, 
    unsigned int chr,
    unsigned int loc,
    enum STRAND strnd,
    float penalty );

void
add_mapped_location(
    mapped_location* loc,
    mapped_locations* locs );

void
add_and_expand_location_from_index(
        mapped_locations* results,
        INDEX_LOC_TYPE* iloc,
        enum STRAND strnd,
        float penalty,
        struct genome_data* genome
    );

void
copy_mapped_location( mapped_location* dest, mapped_location* src );

void
print_mapped_locations( mapped_locations* results );

/*
 *  END Mapped Location
 *
 **************************************************************************/

/***** ml_match *****/

/* Stores a matching set of mapped_locations */
struct ml_match {
    mapped_location** locations;
    int* subseq_lengths;
    int* subseq_offsets;
    int len;

    /* The cumulative length of gaps in the reference genome for this set of
     * matched mapped_locations */
    int cum_ref_gap;
};

void
init_ml_match( struct ml_match** match, int match_len );

void
free_ml_match( struct ml_match* match );

void
add_location_to_ml_match(
        mapped_location* location,
        struct ml_match* match, 
        int subseq_length,
        int subseq_offset,
        int location_index,
        int cum_ref_gap );

enum bool
ml_match_is_valid( struct ml_match* match );

/***** ml_matches ******/

#define ML_MATCHES_GROWTH_FACTOR 10
struct ml_matches {
    int allocated_length;
    int length;
    struct ml_match* matches;
};

void
init_ml_matches( struct ml_matches** matches );

void
free_ml_matches( struct ml_matches* matches );

void
add_ml_match(
        struct ml_matches* matches,
        struct ml_match* match );

/***** ml_match_stack *****/

#define MAX_ML_MATCH_STACK_LEN 200
struct ml_match_stack {
    struct ml_match stack[MAX_ML_MATCH_STACK_LEN];
    int top;
};

void
init_ml_match_stack(
        struct ml_match_stack** stack );

void
free_ml_match_stack(
        struct ml_match_stack* stack );

enum bool
ml_match_stack_is_empty(
        struct ml_match_stack* stack );

void
ml_match_stack_push(
        struct ml_match_stack* stack,
        struct ml_match* match );

struct ml_match
ml_match_stack_pop(
        struct ml_match_stack* stack );

#endif /* define MAPPED_LOCATION */
