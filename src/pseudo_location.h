#include "config.h"

/* realloc with error checking */
inline void* 
realloc_CE( void* ptr, size_t size );

/*
 * A single pseudo location. This stores all of the real 
 * locations associated with a pseudo location.
 */

struct pseudo_location_t {
    int num;
    GENOME_LOC_TYPE* locs;
};

/* We assume that space has already been allocated for the struct */
inline void
init_pseudo_location( struct pseudo_location_t* loc );

inline void 
free_pseudo_location( struct pseudo_location_t* loc );

inline void
add_loc_to_pseudo_location( 
    struct pseudo_location_t* ps_loc, const GENOME_LOC_TYPE* const loc);

/*
 * Pseudo locations container
 *
 */

struct pseudo_locations_t {
    int num;
    struct pseudo_location_t* locs;
};

void
init_pseudo_locations( 
    struct pseudo_locations_t** locs );
    
void
free_pseudo_locations( struct pseudo_locations_t* locs );

void
sort_pseudo_locations( struct pseudo_locations_t* locs );

/* return the index of the new entry */
unsigned int 
add_new_pseudo_location( struct pseudo_locations_t* locs );
