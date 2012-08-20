#ifndef SAM_HEADER
#define SAM_HEADER

#include "config.h"
#include "read.h"

struct mapped_read_t;
struct genome_data;
struct rawread;
struct rawread_db_t;
struct mapped_reads_db;
struct genome_data;
struct cond_prbs_db_t;
enum bool; 

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* mpd_rd,
    struct cond_prbs_db_t* cond_prbs_db,    
    struct genome_data* genome,
    struct read* r,
    enum bool expand_pseudo_locations
);

void
write_nonmapping_reads_to_fastq( 
        struct rawread_db_t* rdb,
        struct mapped_reads_db* mappings_db
    );

void
write_mapped_reads_to_sam( 
        struct rawread_db_t* rdb,
        struct mapped_reads_db* mappings_db,
        struct cond_prbs_db_t* cond_prbs_db,
        struct genome_data* genome,
        enum bool reset_cond_read_prbs,
        /* whether or not to print out pseudo locations
           as real locations, or to print out each real loc
           that makes ups the pseudo location */
        enum bool expand_pseudo_locations,
        FILE* sam_ofp
    );

#endif // SAM_HEADER
