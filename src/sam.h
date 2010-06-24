
struct mapped_read_t;
struct genome_data;
struct rawread;
struct rawread_db_t;
struct mapped_reads_db;
struct genome_data;
enum bool; 

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* pkd_rd,
    struct genome_data* genome,
    struct rawread* rr1,
    struct rawread* rr2
);

void
write_mapped_reads_to_sam( 
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mappings_db,
    struct genome_data* genome,
    enum bool reset_cond_read_prbs,
    FILE* sam_ofp );

