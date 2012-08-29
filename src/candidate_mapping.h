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
    /* whether this immediately follows a gap between read subtemplates */
    enum bool follows_ref_gap;
    /* pos of this candidate mapping relative to the subset of candidate
     * mappings it belongs to, which is derived from one subtemplate */
    int pos;
};

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

} candidate_mapping;


/* This is a list of pointers to candidate_mapping, where sets of pointers may
 * be separated by NULL pointers to indicate subsets of candidate mappings. */
typedef struct {
    /* the number of allocated pointers, including NULL pointers */
    int length;
    /* the number of candidate mappings */
    int num_candidate_mappings;

    candidate_mapping** mappings;
} candidate_mappings;

/* Deal with candidate mappings arrays */

/* Initialize the array */
void
init_candidate_mappings( candidate_mappings** mappings );

void
add_candidate_mapping( candidate_mappings* mappings,
                       candidate_mapping* mapping     );

void
add_null_separator_to_candidate_mappings( candidate_mappings* mappings );

void
free_candidate_mappings( candidate_mappings* mappings,
                         enum bool free_mappings );

candidate_mapping*
init_candidate_mapping_from_read_subtemplate(
        struct read_subtemplate* rst,
        float max_penalty_spread );

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
