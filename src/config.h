/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef CONFIG
#define CONFIG

/****** version options                         *******/
#define GENOME_FILE_FORMAT_VERSION   0

/****** verbosity options                       ******/
/* how often we print out the mapping status */
#define MAPPING_STATUS_GRANULARITY 100000

/****** configuration options                   ******/

#define N_DEBUG
#define PROFILE_CANDIDATE_MAPPING true

/* Try to map 99% of the input reads with the estimated error model */
#define DEFAULT_ESTIMATED_ERROR_METAPARAMETER 0.99
/* Allow up to 10% of the bp's in a read to be mismatches in a mapping */
#define DEFAULT_MISMATCH_METAPARAMTER 0.1

/* Mismatch rates for the estimated error model bootstrap */
#define MAX_NUM_MM_RATE 0.10
#define MAX_NUM_MM_SPREAD_RATE 0.05

#define HIGH_QUALITY_BP_ERROR_PRB 1e-2

// the num of reads that we map before updating the error estimates
#define READS_STAT_UPDATE_STEP_SIZE 50000

// Iterative Mapping
#define MAX_NUM_EM_ITERATIONS 500
#define MAX_PRB_CHANGE_FOR_CONVERGENCE 1e-2
#define EXPLORATION_PRIOR 1e-10

// Wiggle output
#define MAX_TRACE_FNAME "max_trace.wig"
#define MIN_TRACE_FNAME "min_trace.wig"

// Samples output 
#define SAM_MARGINAL_OFNAME "mapped_reads.sam"
#define SAM_MARGINAL_NC_OFNAME "mapped_reads.nc.sam"

/********** Iterative Mapping **************************************************/
// this says that we will setup the iterative mapping framework,
// but not actually take any samples ( they can be taken later with 
// sample mappings )
#define DEFAULT_NUM_SAMPLES 0
#define SAVE_STARTING_SAMPLES true
#define STARTING_SAMPLES_PATH "./starting_samples/"
#define STARTING_SAMPLES_META_INFO_FNAME "./starting_samples/meta_info.csv"

#define SAVE_SAMPLES true
#define RELAXED_SAMPLES_PATH "./samples/"
#define RELAXED_SAMPLES_META_INFO_FNAME "./samples/meta_info.csv"

#define SAVE_BOOTSTRAP_SAMPLES true
#define SAVE_AGGREGATED_BOOTSTRAP_SAMPLES (false && SAVE_BOOTSTRAP_SAMPLES)
#define NUM_BOOTSTRAP_SAMPLES 25
#define BOOTSTRAP_SAMPLES_PATH "./bootstrap_samples/"
#define BOOTSTRAP_SAMPLES_MAX_PATH "./bootstrap_samples/max_traces/"
#define BOOTSTRAP_SAMPLES_MIN_PATH "./bootstrap_samples/min_traces/"
#define BOOTSTRAP_SAMPLES_ALL_PATH "./bootstrap_samples/all_traces/"

#define CALL_PEAKS true
#define CALLED_PEAKS_OUTPUT_DIRECTORY "./called_peaks/"
#define JOINED_CALLED_PEAKS_FNAME "./called_peaks/peaks.wig"
#define CALLED_PEAK_REGIONS_FNAME "./called_peaks/peaks.bed"

/*******************************************************************************
 *
 * default filenames 
 *
 *
 ***/

#define MAPPED_READS_DB_FNAME "mapped_reads.db"
#define MAPPED_NC_READS_DB_FNAME "mapped_NC_reads.db"

#define CANDIDATE_MAPPINGS_DB_FNAME "candidate_mappings"
#define CANDIDATE_MAPPINGS_NC_DB_FNAME "NC_candidate_mapings"

#define PSEUDO_LOCATIONS_FNAME "pseudo_locations.txt"

#define GENOME_FNAME "genome.bin"
#define GENOME_INDEX_FNAME "genome.bin.index"
#define GENOME_INDEX_PSLOCS_FNAME "genome.bin.index.pslocs"
#define GENOME_INDEX_DIPLOID_MAP_FNAME "genome.bin.index.dmap"

#define ERROR_STATS_LOG "error_stats.log"
#define BOOTSTRAP_ERROR_STATS_LOG "bootstrap_error_stats.log"


/**** set global constantsa for maximum read length, etc.       ****/

#define MAX_READ_LEN 255
#define MAX_READNAME_LEN 255
#define MAX_QUAL_SCORE 255

/**** determine how the letters are packed       ****/

#define PSEUDO_LOC_MIN_SIZE 50
#define EXPAND_UNPAIRED_PSEUDO_LOCATIONS true
#define EXPAND_PAIRED_PSEUDO_LOCATIONS true

typedef unsigned char LEVEL_TYPE;
/* TODO this probably isn't big enough - max read len of 255? */
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
 * The assay that generated these reads.
 * 
 * This only matters if we want to do 
 * iterative mapping, and for the wiggle
 * output type.
 *
 */
enum assay_type_t {
    // UNKNOWN = 0,
    CAGE = 1,
    CHIP_SEQ = 2,
    RNA_SEQ = 3
};

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

/*
 * if CHR_BITS or LOCATION_BITS are changed, INDEX_LOC_TYPE and the
 * mapped_read_t pseudo structure should be updated so they remain byte-aligned
 */
#define CHR_BITS 16
#define CHR_NUM_MAX (65536 - 1) // 2**16 - 1
#define PSEUDO_LOC_CHR_INDEX 0
#define LOCATION_BITS 29
#define LOCATION_MIN 0
#define LOCATION_MAX (536870912 - 1) // 2**29 = 536870912

#define FRAGMENT_LENGTH_BITS 20
#define FRAGMENT_LENGTH_MIN  0
#define FRAGMENT_LENGTH_MAX ( 524288 - 1 ) // 2**20 = 1048576

/* Reasonable default (TODO make user-adjustable?) */
#define REFERENCE_INSERT_LENGTH_MAX 500000

/* this needs to always be able to store up to LOCATION_MAX */
#define SIGNED_LOC int

/* if this structure is changed, GENOME_FILE_FORMAT_VERSION should be
 * incremented to reflect backwards incompatible changes in the binary
 * genome/index format */
typedef struct __attribute__((__packed__))
{
    unsigned is_paternal    :1;
    unsigned is_maternal    :1;
    
    unsigned unused_space   :1;
    
    /* the chr that the read came from */
    unsigned chr            :CHR_BITS;

    /* 
     * Start of the sequence in the 5' direction - 
     * so a read of length L that maps here covers
     * bps (loc, loc+1, ..., loc+L-1)
     */
    unsigned loc            :LOCATION_BITS;

} INDEX_LOC_TYPE;

enum CHR_SOURCE
{
    // UNKNOWN = 0,
    REFERENCE = 1,
    PATERNAL  = 2,
    MATERNAL  = 3
};

enum input_file_type_t {
    // We have not specified the file type
    // UNKNOWN = 0,
    // Fastq file. Quality scores are PHRED 
    // ( offset 40, log10 prob )
    SANGER_FQ = 1,
    // Illumina version 1.3
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    ILLUMINA_v13_FQ = 2,
    // Solexa version 1.4+
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    SOLEXA_v14_FQ = 2,
    // Solexa version 1.5+
    // Fastq file. Quality scores are PHRED 
    // ( offset 64, log10 prob )
    // In addition, scores of 0 and 1 are never used
    // Furthermore, a score of 2 is special: 
    /* 
     *  "A quality score of 2, encoded as a "B", is used as a special indicator. A
     *   quality score of 2 does not imply a specific error rate, but rather implies 
     *   that the marked region of the read should not be used for downstream analysis."
     *   http://docs.google.com/fileview?id=0B-lLYVUOliJFYjlkNjAwZjgtNDg4ZC00MTIyLTljNjgtMmUzN2M0NTUyNDE3&hl=en&pli=1
     */
    ILLUMINA_v15_FQ = 2,
    // Solexa version 1.0
    // Fastq file. Quality scores are 
    // ( offset 64, log10 log odds )
    SOLEXA_LOG_ODDS_FQ = 3,
    // Test suite error format
    // All scores are perfect h's, 
    // assume it's the test code, and the 
    // scores are ILLUMINA_v13_FQ
    TEST_SUITE_FORMAT = 2,
    // Marks SOLEXA ( from his ChIP-seq )
    // Maximum is ASCII V
    MARKS_SOLEXA = 4,
    // Illumina 1.8+
    ILLUMINA_v18_FQ = 1,
};

/*** Globally useful macros ***/
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))


#endif
