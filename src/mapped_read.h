/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef MAPPED_READ
#define MAPPED_READ

#include <stdio.h>

#include "config.h"
#include "candidate_mapping.h"
#include "rawread.h"


struct fragment_length_dist_t;
struct genome_data;
struct rawread_db_t;

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


/*************************************************************************
 *
 *  Mapped Reads
 * 
 *  Candidate Mappings that have been joined.
 *
 */

#define ML_PRB_TYPE float

/* TODO - make this a bitfield */

/* Set if the read that contribute to this mapped location are paired */
#define IS_PAIRED 1
/* Set if the first read ( ie read_name/1 ) in the pair was rev 
   complemented to map */
#define FIRST_READ_WAS_REV_COMPLEMENTED 2
/* Set if the first read ( ie read_name/1 ) maps to start_pos */
#define FIRST_PAIR_IS_FIRST_IN_GENOME 4
/* set if the first read covers a snp */
#define FIRST_READ_COVERS_SNP 8
/* set if the second read covers a snp */
#define SECOND_READ_COVERS_SNP 16

struct mapped_read_location {
    // if stop_pos > start_pos, the strand is pos
    // else, the strand is negative. For rna-seq, 
    // the strand will be determined by the exon
    // that it maps into
    // enum STRAND strand;
    unsigned char flag;
    unsigned snps_bm_r1 :MAX_NUM_SNPS;
    unsigned snps_bm_r2 :MAX_NUM_SNPS;
    unsigned char chr;
    unsigned int start_pos;
    /* THIS IS EXCLUSIVE, ie NOT including stop */
    unsigned int stop_pos;
    ML_PRB_TYPE seq_error;
    /* the probability of observing this fragment length */
    ML_PRB_TYPE fl_prob;
    ML_PRB_TYPE cond_prob;
} __attribute__((__packed__)) ;

struct mapped_read_t {
    unsigned long read_id;
    unsigned short num_mappings;
    struct mapped_read_location* locations;
};

typedef struct {
    size_t size;
    size_t allocated_size;
    struct mapped_read_t* reads;
} mapped_reads;

unsigned char
chr_index( char* chr_name );

void
init_mapped_read( struct mapped_read_t** rd );

void
free_mapped_read( struct mapped_read_t* rd );

void
add_location_to_mapped_read( 
    struct mapped_read_t* rd, struct mapped_read_location* loc );

void
reset_read_cond_probs( struct mapped_read_t* rd  );

void
fprintf_mapped_read( FILE* fp, struct mapped_read_t* r );

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* pkd_rd,
    struct genome_data* genome,
    struct rawread* rr1,
    struct rawread* rr2
);

void
build_mapped_read_from_candidate_mappings( 
    candidate_mappings* mappings, 
    struct mapped_read_t** mpd_rd,
    long read_id );


int 
write_mapped_read_to_file( struct mapped_read_t* read, FILE* of  );


/*
 *  END Mapped Reads
 *
 **************************************************************************/


/*****************************************************************************
 *
 * Mapped Reads DB Code
 *
 */

struct mapped_reads_db {
    FILE* fp;
    /* Set this to locked when we mmap it - 
       then forbid any new writes to the 
       file 
    */
    enum bool locked;
    /* mmap data */
    /* pointer to the mmapped data and its size in bytes */
    char* mmapped_data;
    size_t mmapped_data_size; 
    
    char** mmapped_reads_starts;
    unsigned long num_mmapped_reads;

    struct fragment_length_dist_t* fl_dist;

    unsigned long current_read;
};

typedef struct {
    /* the current position of the fp */
    long int fpos;
    struct mapped_reads_db* rdb;
} mapped_reads_db_cursor;

void
init_mapped_reads_db( struct mapped_reads_db** rdb, char* fname );

void
build_fl_dist_from_file( struct mapped_reads_db* rdb, FILE* fl_fp );

void
close_mapped_reads_db( struct mapped_reads_db* rdb );

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    struct mapped_read_t* rd
);

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb );

enum bool
mapped_reads_db_is_empty( struct mapped_reads_db* rdb );

/* returns EOF, and sets rd to NULL if we reach the end of the file */
int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* rdb, 
    struct mapped_read_t** rd );

void
set_all_read_fl_probs( struct mapped_reads_db* rdb );

void
reset_all_read_cond_probs( struct mapped_reads_db* rdb );

void
write_mapped_reads_to_sam( struct rawread_db_t* rdb,
                           struct mapped_reads_db* mappings_db,
                           struct genome_data* genome,
                           enum bool reset_cond_read_prbs,
                           FILE* sam_ofp );

/*
 * this requires code from iterative mapping to write out the
 * traces. It also rewquires the mapped_reads to be mmapped and 
 * indexed. 
 *
 * TODO - rewrite this to not require the mmapped mapped reads.
 * ( Why? Just 1 less dependency. Also, it confuses the order
 *   and makes it hard. Although, actually, we could probably
 *   just move the fp and read sicne we arent going to actually 
 *   change anything on disk. hmmm... )
 *
 */

void
write_marginal_mapped_reads_to_stranded_wiggles( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    FILE* fwd_wfp, FILE* bkwd_wfp );

void
mmap_mapped_reads_db( struct mapped_reads_db* rdb );

void
munmap_mapped_reads_db( struct mapped_reads_db* rdb );

void
index_mapped_reads_db( struct mapped_reads_db* rdb );


/*
 *  END Mapped Reads DB
 *
 **************************************************************************/

#endif // MAPPED_READ
