/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef CONFIG
#define CONFIG

/****** configuration options                   ******/

// Candidate Mappings
#define DEFAULT_MIN_MATCH_PENALTY -7.0
#define DEFAULT_MAX_PENALTY_SPREAD 2.1

// Iterative Mapping
#define MAX_NUM_EM_ITERATIONS 500
#define MAX_PRB_CHANGE_FOR_CONVERGENCE 1e-2
#define EXPLORATION_PRIOR 1e-10

// Wiggle output
#define MAX_TRACE_FNAME "max_trace.wig"
#define MIN_TRACE_FNAME "min_trace.wig"

// Samples output 
// Bootstrap 
#define SAVE_BOOTSTRAP_SAMPLES true
#define SAVE_AGGREGATED_BOOTSTRAP_SAMPLES (false && SAVE_BOOTSTRAP_SAMPLES)
#define NUM_BOOTSTRAP_SAMPLES 25
#define BOOTSTRAP_SAMPLES_PATH "./bootstrap_samples/"
#define BOOTSTRAP_SAMPLES_MAX_PATH "./bootstrap_samples/max_traces/"
#define BOOTSTRAP_SAMPLES_MIN_PATH "./bootstrap_samples/min_traces/"
#define BOOTSTRAP_SAMPLES_ALL_PATH "./bootstrap_samples/all_traces/"


#define SAVE_STARTING_SAMPLES true
#define STARTING_SAMPLES_PATH "./starting_samples/"
#define STARTING_SAMPLES_META_INFO_FNAME "./starting_samples/meta_info.csv"

#define SAVE_SAMPLES true
#define RELAXED_SAMPLES_PATH "./samples/"
#define RELAXED_SAMPLES_META_INFO_FNAME "./samples/meta_info.csv"

#define MAPPED_READS_DB_FNAME "mapped_reads.db"

#define PSEUDO_LOCATIONS_FNAME "pseudo_locations.txt"

/**** determine how the letters are packed       ****/

#define PSEUDO_LOC_MIN_SIZE 100
#define EXPAND_UNPAIRED_PSEUDO_LOCATIONS false
#define EXPAND_PAIRED_PSEUDO_LOCATIONS true

typedef unsigned char LEVEL_TYPE;
typedef unsigned char READ_POSITION;

// the length of each index letter, in bps 
// this must correspond with LETTER_TYPE 
// 4 is the best option for conserving memory, 
// but for small genomes smaller values may
// provide a small increase in speed
#define LETTER_LEN 4
typedef unsigned char LETTER_TYPE;
// the number of sequences that fit into 
// a letter. this should always be
// 4^LETTER_LEN
#define ALPHABET_LENGTH 256

typedef unsigned char byte;

enum bool {
    false = 0,
    true = 1
};

/*
 *  A few commonly used type definitions
 */

/* For use in the below enums */
#define UNKNOWN 0

typedef char NODE_TYPE;

/* 
 * What strand a read came from. For raw reads, this will usually be unknown
 * ( unless the assay was designed to determine strand ). For mapped reads,
 * this is typically BKWD if we needed to take the reverse complement and 
 * FWD if not. Note that this does not indicate where it actually came from, 
 * just where it mapped to. For some reads, ( ie canonical junctions ) this is 
 * a known quantity ( ie, we know the gene was on the 5' strand )
 *
 */
enum STRAND {
    // UNKNOWN = 0,
    FWD = 1,
    BKWD = 2
};

/*
 *  If the read was paired end and, if so, which end it came from. Unknown 
 *  should never happen, but it's included for consistency 
 */

enum READ_END {
    // UNKNOWN = 0,
    NORMAL = 1,
    FIRST = 2,
    SECOND = 3
};

/* this should be used for NULL, in the below as commented */
#define UNKNOWN 0

enum JUNCTION_TYPE
{
    /*  0 = 'Not Set' - it should be the initial setting */
    PRE_INTRON = 1,
    POST_INTRON = 2,
    NON_CANONICAL = 3
};

enum READ_TYPE
{
    /*  0 = 'Not Set' - it should be the initial setting */
    SINGLE_END = 1,
    /* Paired end, with an indeterminate end */
    PAIRED_END = 2,
    /* Paired end, but this is end 1 ( ie, KEY/1 in the read name ) */
    PAIRED_END_1 = 3,
    /* Paired end, but this is end 2 ( ie, KEY/2 in the read name ) */
    PAIRED_END_2 = 4
};

enum RECHECK
{
    /* the default, it needs a full recheck */
    RECHECK_FULL = 0,
    /* recheck to make sure the penalties are above the min */
    RECHECK_PENALTY = 1,
    /* recheck against chr to make sure the penalty is correct */
    RECHECK_LOCATION = 2,
    /* we know this crosses a snp, so we test for one */
    COVERS_SNP = 3,
    /* all rechecks have been performed, and the loc is valid */
    VALID = 4,
    /* one or more rechecks failed */
    INVALID = 5
};


/*
 * THIS SECTION IS OBSOLETE *************************************
 
   Store potential junctions
   We just keep track of which strand the gene that could
   have been at this location could have come from.

   ie, if the genome location seq is CCCGT|AGTTT, then this read 
   could have been on the + strand ( assuming, of course, that 
   the genomic sequence was 5' -> 3' ) and if the genomic seq
   was AGTTT, then the read CCCTTT could have mapped in the
   positive direction. 
   
   So, when we observe the sequence CCCGT in chr 0 at bp 0, 
   the corresponding genome location is given by 
       gene_strnd = 1
       intron_start = 1 // before the intron
       chr = 0
       bp = 0
       seq = CCC
   
   In addition, when we observe AGTTT  we add the end  
       gene_strnd = 1
       intron_start = 0 // after the intron
       chr = 0
       bp = 0 + 3 + intron_length
       seq = TTT
    
    sequence is 
    
    5'+  GGC|CT|AC|CAA 3'
    3'-  CCG|GA|TG|GTT 5'

    then observing GGCCAA corresponds with observing a - strand gene 
    from a + read. So, we would add GGC here as a potential gene *end*
    for a *- strand gene* and CAA as a gene *start*
*/
#define CHR_BITS 15
#define CHR_NUM_MAX (32768 - 1) // 2**15 - 1
#define PSEUDO_LOC_CHR_INDEX 0
#define LOCATION_BITS 28
#define LOCATION_MAX (268435456 - 1) // 2**28 = 268435456
#define FRAGMENT_LENGTH_BITS 20
#define FRAGMENT_LENGTH_MIN ( -524288 + 1 ) // 2**20 = 1048576
#define FRAGMENT_LENGTH_MAX ( 524288 - 1 )
#define MAX_NUM_SNPS 3
typedef struct __attribute__((__packed__))
{
    /* 1 if it covers at least 1 snp */
    unsigned covers_snp     :1;
    unsigned snp_coverage   :MAX_NUM_SNPS;

    /* read_type 0 = normal, 1 = junction */
    unsigned read_type      :1;

    /* the chr that the read came from */
    unsigned chr            :CHR_BITS;

    /* 
     * Start of the sequence in the 5' direction - 
     * so a read of length L that maps here covers
     * bps (loc, loc+1, ..., loc+L-1)
     */
    unsigned loc            :LOCATION_BITS;

} GENOME_LOC_TYPE;

/*** Globally useful macros ***/
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))


#endif
