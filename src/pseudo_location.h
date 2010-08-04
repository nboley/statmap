#include "config.h"

/*
 * A single pseudo location. This stores all of the real 
 * locations associated with a pseudo location.
 */

struct pseudo_location_t {
    int num;
    GENOME_LOC_TYPE* locs;
};

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

void
fprint_pseudo_locations( FILE* of, struct pseudo_locations_t* ps_locs );

void
load_pseudo_locations( FILE* fp, struct pseudo_locations_t** ps_locs );

size_t
write_pseudo_locations_to_file( struct pseudo_locations_t* ps_locs, FILE* of );

void
load_pseudo_locations_from_mmapped_data( 
    struct pseudo_locations_t** ps_locs, char* data );

/* return the index of the new entry */
unsigned int 
add_new_pseudo_location( struct pseudo_locations_t* locs );
