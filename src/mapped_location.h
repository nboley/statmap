/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef MAPPED_LOCATION
#define MAPPED_LOCATION

#include <stdio.h>

#include "config.h"
#include "genome.h"
#include "rawread.h"

/* The mapped location probability type */
#define ML_PRB_TYPE float

/* block size for newly allocated results memory */
#define RESULTS_GROWTH_FACTOR 100

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
    float penalty;
} mapped_location;

typedef struct {
    /* An array of mapped locations */
    mapped_location* locations;
    size_t length;
    size_t allocated_length;
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
                     float penalty );

void
print_mapped_locations( mapped_locations* results );

/*
 *  END Mapped Location
 *
 **************************************************************************/

/*************************************************************************
 *
 *  Candidate Mapping
 *
 *  These are after a read has gone to the index,  had the read
 *  information extracted, and probably gone through any necessary recheck.
 *  However, they have not been paired to make them mapped_read's.
 *
 */

typedef struct __attribute__((packed))__{
    /* 
     * Whether or not we need to go to the genome to update this 
     * information. 
     *
     * If this is set to NO, then the chr, start_bp,
     * and penalty are all assumed correct and subseq_offset = 0
     * and subseq_len = rd_len. ( which is not checked w/o asserts )
     * 
     * If this is YES, then we go back to chr in the genome at loc
     * start_bp. Then we recalculate the penalty from the genome
     * directly and, update the penalty for that location. 
     */
    enum RECHECK recheck;

    /*** Info relating to the location ***/
    /* the chromosome code */
    unsigned char chr;
    /* the bp in genomic coordinates */
    unsigned int start_bp;

    /*** Info related to the read ***/
    /* If this is a normal, or paired end read */
    enum READ_TYPE rd_type;
    /* the full length in bp's of the underlying read */
    READ_POSITION rd_len;

    /*** Info related to the mapped read ***/
    /* 
     * which strand direction this read is ( ie, it is fwd
     * if it maps to the fwd stranded genome, and rev if it's 
     * rev complement maps to the fwd genome )
     */
    enum STRAND rd_strnd;
    /* the penalty associated with this match */
    /* 
     * generally, this is the log10 probability of observing the sequence
     * from this genomic location given the sequenced sequence.
     */
    float penalty;

    /*** Info relating to properties of the subseq loc ***/
    /* 
     * The offset of the partial match. When we map junction reads
     * we are mapping the beggining k basepairs and the final N-k
     * bp's and then joining them together into a junction read if
     * the associated meta data suggests that this is possible. Thus,
     * subseq offset tells us where, in relation to the underling full 
     * read, that this junction read comes from. 
     *
     * For now, it should be 0 or N-k, but in the future if reads get 
     * much longer and could cross multiple junctions, this need not 
     * be the case 
     */
    READ_POSITION subseq_offset;
    /* The length of the underlying sub match */
    READ_POSITION subseq_len;     
} candidate_mapping;


#define CANDIDATE_MAPPINGS_GROWTH_FACTOR 10
typedef struct {
    unsigned int allocated_length;
    unsigned int length;
    candidate_mapping* mappings;
} candidate_mappings;

/* Deal with candidate mappings arrays */

/* Initialize the array */
void
init_candidate_mappings( candidate_mappings** mappings );

candidate_mapping
init_candidate_mapping_from_template( rawread* rp, 
                                      int indexed_seq_len, 
                                      float max_penalty_spread );

/* add a copy of a candidate mapping */
void
add_candidate_mapping( candidate_mappings* mappings,
                       candidate_mapping* mapping     );

void
free_candidate_mappings( candidate_mappings* mappings );


void
print_candidate_mapping( candidate_mapping* mapping );

void
print_candidate_mappings( candidate_mappings* mappings );

inline int
cmp_candidate_mappings( const candidate_mapping* m1, const candidate_mapping* m2 );

int 
sort_candidate_mappings( candidate_mappings* mappings );

/*
 *  END Candidate Mapping
 *
 **************************************************************************/

/*************************************************************************
 *
 *  Packed Mapped Reads
 * 
 *  Reads that have been joined, but unlike mapped reads proper have not
 *  had the 'extra' read inforamtion attached.
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

typedef struct  __attribute__((__packed__)) {
    // if stop_pos > start_pos, the strand is pos
    // else, the strand is negative. For rna-seq, 
    // the strand will be determined by the exon
    // that it maps into
    // enum STRAND strand;
    unsigned char flag;
    unsigned char chr;
    unsigned int start_pos;
    unsigned int stop_pos;
    ML_PRB_TYPE seq_error;
    ML_PRB_TYPE cond_prob;
} mapped_read_location;

typedef struct {
    unsigned long read_id;
    unsigned short num_mappings;
    mapped_read_location* locations;
} mapped_read;

typedef struct {
    size_t size;
    size_t allocated_size;
    mapped_read* reads;
} mapped_reads;

unsigned char
chr_index( char* chr_name );

void
init_mapped_read( mapped_read** rd );

void
free_mapped_read( mapped_read* rd );

void
add_location_to_mapped_read( 
    mapped_read* rd, mapped_read_location* loc );

void
fprintf_mapped_read( FILE* fp, mapped_read* r );

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    mapped_read* pkd_rd,
    genome_data* genome,
    rawread* rr1,
    rawread* rr2
);

void
build_mapped_read_from_candidate_mappings( 
    candidate_mappings* mappings, 
    mapped_read** mpd_rd,
    unsigned long read_id );

int 
write_mapped_read_to_file( mapped_read* read, FILE* of  );


/*
 *  END Mapped Reads
 *
 **************************************************************************/


/*****************************************************************************
 *
 * Mapped Reads DB Code
 *
 */

typedef struct {
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
    
} mapped_reads_db;

typedef struct {
    /* the current position of the fp */
    long int fpos;
    mapped_reads_db* rdb;
} mapped_reads_db_cursor;

void
init_mapped_reads_db( mapped_reads_db** rdb, char* fname );

void
close_mapped_reads_db( mapped_reads_db* rdb );

void
add_read_to_mapped_reads_db( 
    mapped_reads_db* rdb,
    mapped_read* rd
);

void
rewind_mapped_reads_db( mapped_reads_db* rdb );

enum bool
mapped_reads_db_is_empty( mapped_reads_db* rdb );

int
get_next_read_from_mapped_reads_db( 
    mapped_reads_db* rdb, 
    mapped_read** rd );

void
mmap_mapped_reads_db( mapped_reads_db* rdb );

void
munmap_mapped_reads_db( mapped_reads_db* rdb );

void
index_mapped_reads_db( mapped_reads_db* rdb );


/*
 *  END Mapped Reads DB
 *
 **************************************************************************/


#endif /* define MAPPED_LOCATION */
