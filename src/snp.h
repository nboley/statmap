#ifndef SNP
#define SNP

typedef struct {
    char alt_bp;
    float alt_frac;
    GENOME_LOC_TYPE loc;    
} snp_t;


int 
cmp_snp( const void* void_snp1, const void* void_snp2 );

void
fprintf_snp( FILE* fp, snp_t* snp );

struct snp_db_t;

void
init_snp_db( struct snp_db_t** snp_db  );

void
free_snp_db( struct snp_db_t* snp_db );

void
sort_snp_db( struct snp_db_t* snp_db );

void
find_snps_in_snp_db( struct snp_db_t* snp_db, int chr, int start, 
                     int stop, /* NOT including stop */ 
                     int* num_snps, int** snp_indexes );

void
fprintf_snp_db( FILE* fp, struct snp_db_t* snp_db );

void
add_snp_to_snp_db( struct snp_db_t* snp_db, snp_t* snp );

void
build_snp_db_from_snpcov_file( FILE* fp, genome_data* genome );

void
index_snp_sequences( struct snp_db_t* snp_db, genome_data* genome );

#endif
