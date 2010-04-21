#ifndef SNP
#define SNP

#include "config.h"
#include "genome.h"

typedef struct {
    char alt_bp;
    float alt_frac;
    GENOME_LOC_TYPE loc;   
    /* The number of reference, and alternate bps observed */
    float alt_cnt;
    float ref_cnt;
} snp_t;


int 
cmp_snp( const void* void_snp1, const void* void_snp2 );

void
fprintf_snp( FILE* fp, snp_t* snp );

/********************************************************************************
 * Methods for dealing with the snp database ( for now just a sorted array )
 *
 */

struct snp_db_t {
    int num_snps;
    snp_t* snps;    
};

void
init_snp_db( struct snp_db_t** snp_db  );

void
free_snp_db( struct snp_db_t* snp_db );

void
sort_snp_db( struct snp_db_t* snp_db );

void
find_snps_in_snp_db( struct snp_db_t* snp_db, int chr, int start, 
                     int stop, /* NOT including stop */ 
                     int* num_snps, snp_t** snps );

void
fprintf_snp_db( FILE* fp, struct snp_db_t* snp_db );

void
add_snp_to_snp_db( struct snp_db_t* snp_db, snp_t* snp );

void
build_snp_db_from_snpcov_file( FILE* fp, struct genome_data* genome );

struct mapped_reads_db;

void
update_snp_estimates_from_candidate_mappings( 
    struct mapped_reads_db* rdb, struct genome_data* genome );

void
write_snps_to_file( FILE* ofp, struct genome_data* genome );

#endif
