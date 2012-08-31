#ifndef SAM_HEADER
#define SAM_HEADER

#include "config.h"
#include "read.h"
#include "mapped_read.h"

/* defines copied from sam tools - these go into the 'flag' field */
/* the read is paired in sequencing, no matter 
   whether it is mapped in a pair */
#define BAM_FPAIRED        1
/* the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/* the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/* the mate is unmapped */
#define BAM_FMUNMAP        8
/* the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/* the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/* this is read1 */
#define BAM_FREAD1        64
/* this is read2 */
#define BAM_FREAD2       128
/* not primary alignment */
#define BAM_FSECONDARY   256
/* QC failure */
#define BAM_FQCFAIL      512
/* optical or PCR duplicate */
#define BAM_FDUP        1024

struct genome_data;
struct rawread;
struct rawread_db_t;
struct genome_data;

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    mapped_read_t* mpd_rd,
    struct cond_prbs_db_t* cond_prbs_db,    
    struct genome_data* genome,
    struct read* r );

void
write_nonmapping_reads_to_fastq( 
        struct rawread_db_t* rdb,
        struct mapped_reads_db* mappings_db );

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
        FILE* sam_ofp );

#endif // SAM_HEADER
