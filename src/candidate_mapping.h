/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef CANDIDATE_MAPPING
#define CANDIDATE_MAPPING

#include "config.h"
#include "read.h"

/* FWD declaration for rawread */
// #include "rawread.h"
struct rawread;
// #include "genome.h"
struct genome_data;

/* block size for newly allocated results memory */
#define CAND_MAPPING_RESULTS_GROWTH_FACTOR 100

int
modify_mapped_read_location_for_index_probe_offset(  
    int read_location,
    const int chr,
    const enum STRAND strnd,
    const int subseq_offset,
    const int subseq_len,
    const int read_len,
    struct genome_data* genome
    ) ;

/*************************************************************************
 *
 *  Candidate Mapping
 *
 *  These are after a read has gone to the index,  had the read
 *  information extracted, and probably gone through any necessary recheck.
 *  However, they have not been paired to make them mapped_read's.
 *
 */

struct READ_TYPE {
    /* Whether this candidate mapping follows a gap in the reference genome
     * (such as an intron) */
    enum bool follows_ref_gap;
    /* Position of this candidate mapping in the underlying template.
     * This is equal to the index of the underlying read subtemplate */
    int pos;
};

struct CIGAR_ENTRY {
    char op;
    int len;
};

/* TODO - for ungapped assays, there will always be just one entry in the
 * cigar string */
#define MAX_CIGAR_STRING_ENTRIES 5

typedef struct __attribute__((packed))__{
    /*** Info relating to the location ***/
    /* the chromosome code */
    short chr;
    /* the bp in genomic coordinates */
    int start_bp;

    /*** Info related to the read ***/
    /* If this is a normal, or paired end read */
    struct READ_TYPE rd_type;
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

    /* cigar string */
    struct CIGAR_ENTRY cigar[MAX_CIGAR_STRING_ENTRIES];
    int cigar_len;

} candidate_mapping;


#define CANDIDATE_MAPPINGS_GROWTH_FACTOR 10
typedef struct {
    int allocated_length;
    int length;
    candidate_mapping* mappings;
} candidate_mappings;

/* Deal with candidate mappings arrays */

/* Initialize the array */
void
init_candidate_mappings( candidate_mappings** mappings );

candidate_mapping
init_candidate_mapping_from_read_subtemplate(
        struct read_subtemplate* rst );

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

void
sort_candidate_mappings( candidate_mappings* mappings );

candidate_mapping
convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
        struct genome_data* genome,
        candidate_mapping cm
    );
/*
 *  END Candidate Mapping
 *
 **************************************************************************/

// fwd declaration of diploid_map_data_t
struct diploid_map_data_t;

void
append_candidate_mappings(
        candidate_mappings* dest,
        candidate_mappings* src
    );

#endif // #ifdef CANDIDATE_MAPPING
