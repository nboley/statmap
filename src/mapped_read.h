/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef MAPPED_READ
#define MAPPED_READ

#define DONT_MALLOC_READS_DB

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "config.h"
#include "statmap.h"
#include "candidate_mapping.h"
#include "rawread.h"

typedef int MPD_RD_ID_T;

struct fragment_length_dist_t;
struct genome_data;
struct rawread_db_t;
struct cond_prbs_db_t;

float
calc_candidate_mapping_penalty(
        candidate_mapping* mapping,
        struct read_subtemplate* rst,
        struct genome_data* genome
    );


/*************************************************************************
 *
 *  Mapped Reads
 * 
 *  Candidate Mappings that have been joined.
 *
 */

/* 
 * Determine whether or not to use the compact representation of ML_PRB_TYPE.
 */

#define NO_COMPACT_ML_PRB_TYPE
#ifdef COMPACT_ML_PRB_TYPE
typedef struct __attribute__((__packed__)) {
    unsigned mantissa :18;
    unsigned exponent :6;
} ML_PRB_TYPE;
#define ML_PRB_MIN 1e-19
#else 
typedef float ML_PRB_TYPE;
#define ML_PRB_MIN FLT_MIN
#define ML_PRB_MAX FLT_MAX
#endif 

static inline float 
ML_PRB_TYPE_to_float( 
    const ML_PRB_TYPE value )
{
#ifdef COMPACT_ML_PRB_TYPE
    if( value.mantissa > 1000 )
        return -1.0;
    
    float rv = value.mantissa;
    rv /= 1000;
    rv *= pow(2, -value.exponent );
    
    assert( rv >= 0 );
    return rv;
#else
    return value;
#endif
}

static inline ML_PRB_TYPE 
ML_PRB_TYPE_from_float( float value  )
{
#ifdef COMPACT_ML_PRB_TYPE
    assert( value == -1 || value == 0 
            || ( value >=  ML_PRB_MIN && value <= 1.0 ) ); 
    
    ML_PRB_TYPE rv = {0,0};

    if( value < -0.5 )
    {
        rv.mantissa = 1024;
        assert( abs( value - ML_PRB_TYPE_to_float(rv) ) < 1e-3  );
        return rv;
    }
    
    int exponent = 0;
    float mantissa = 0;
    
    mantissa = frexpf( value, &exponent );
    assert( mantissa >= 0 );
    assert( mantissa <= 1 );
    
    /* Fix the exponent between -63 and 0 */
    rv.exponent = -(MAX( -64, MIN( exponent, 0  ) ));
    rv.mantissa = MAX( 0, (int) (mantissa*1000) );

    assert( abs( value - ML_PRB_TYPE_to_float(rv) ) == 0  );

    return rv;
#else
    return value;
#endif
}

/*
 * Compact representation of a full mapped read.
 *
 * We assume that each mapped_read_t has at least
 * 1) 1 read_id_node (for now, it has *exactly* 1 read_id_node)
 * 2) 1 mapped_read_location, which has at least
 *    a) 1 subread
 *
 * Note that we use signed bitfields for except single-fit fields. Signed
 * bitfields are much better for debugging, but signed 1-bit fields don't make
 * any sense (two's complement means the only two values it can have are 0 and
 * -1). 
 *
 *  Depending on the assay type, the sublocations may represent fragments in a
 *  paired end assay (originating from read subtemplates) or junctions in
 *  a gapped alignment (originating from indexable subtemplates). The
 *  next_subread_is flags handle all possibilities.
 *
 */

/*
mapped_read_t {
    // 4 bytes
    struct {
        signed read_id    :READ_ID_BITS;          // 31 bits
        unsigned are_more :1;
    } read_id_node;
    // potentially more read ids
    
    // MIN 16 bytes ( 4 per extra read id, 6 per extra subtemplate )
    struct mapped_read_location {
        // 6 bytes
        struct mapped_read_location_prologue{
            unsigned chr            :CHR_BITS;      // 15 bits
            unsigned strand         :1;
            unsigned are_more       :1;
            unsigned trimmed_length :6;
            unsigned unused_bits    :2;             // total: 2 bytes
        
            ML_PRB_TYPE seq_error[NUM_READ_IDS];    // 4 bytes
            // 1 sequencing error for each read id
        }
        
        // 6 bytes
        struct mapped_read_sublocation {
            signed start_pos    :LOCATION_BITS; ( 29 )
            unsigned length     :SUBTEMPLATE_LENGTH_BITS; ( 15 )

            unsigned rev_comp                 :1; // 0 = no, 1 = yes
            unsigned is_full_contig           :1; // unused for now
            unsigned next_subread_is_gapped   :1; // gapped
            unsigned next_subread_is_ungapped :1; // paired end
        }
        // potentially more sublocations
    }
    
    // potentially more mapped read locations

} __attribute__((__packed__));
*/

typedef void mapped_read_t;
typedef void mapped_read_location;

#define MRL_CHR_TYPE unsigned short
#define MRL_START_POS_TYPE unsigned int
#define MRL_STOP_POS_TYPE unsigned int
#define MRL_FL_TYPE unsigned short
#define MRL_TRIM_TYPE unsigned char

#define SUBTEMPLATE_LENGTH_BITS 15
#define MAX_SUBTEMPLATE_LENGTH 32767 // 2**15 - 1

#define READ_ID_BITS 31
#define MAX_READ_ID 1073741823 // 2**31 / 2 - 1 = 1073741823

#define TRIMMED_LENGTH_MAX 63 // 2**6 - 1 = 63
/*
 * Structures used to access and manipulate the mapped_read_t and
 * mapped_read_location pseudo structs
 *
 * Note - the following structures must all be aligned on byte boundaries so we
 * can use sizeof and pointer arithmetic to work with the mapped_read_t pseudo
 * structure
 */

// 4 bytes
typedef struct {
    unsigned read_id  :READ_ID_BITS;
    unsigned are_more :1;
} read_id_node;

// 6 bytes
typedef struct {
    unsigned chr            :CHR_BITS;
    unsigned strand         :1;
    unsigned are_more       :1;
    unsigned trimmed_length :6;
    //unsigned unused_bits    :0;
    /* For now, we assume that there is only 1 readid and therefore only
     * 1 corresponding seq_error */
    ML_PRB_TYPE log_seq_error;      
} mapped_read_location_prologue;

// 6 bytes
typedef struct {
    signed start_pos    :LOCATION_BITS;
    unsigned length     :SUBTEMPLATE_LENGTH_BITS;
    
    unsigned rev_comp                 :1;
    unsigned is_full_contig           :1;
    unsigned next_subread_is_gapped   :1;
    unsigned next_subread_is_ungapped :1;
} mapped_read_sublocation;

/* Utility functions for accessing components of mapped_read_t */

mapped_read_sublocation*
get_start_of_sublocations_in_mapped_read_location(
        const mapped_read_location* loc );

int
get_start_for_sublocations_group(
        mapped_read_sublocation* group );

int
get_stop_for_sublocations_group(
        mapped_read_sublocation* group );

enum bool
last_sublocation_in_mapped_read_location( mapped_read_sublocation* sub );

enum bool
end_of_sublocations_group( mapped_read_sublocation* sub_loc );

/* 
 * small, inline functions for accessing the items in mapped_read_location.
 * These are how elements should be accessed - there is no guarantee that
 * mapped_read_location will remain the same 
 *
 * These are defined in the header so that they can be inlined.
 *
 */

/** CHR **/

static inline MRL_CHR_TYPE
get_chr_from_mapped_read_location( const mapped_read_location* loc)
{
    mapped_read_location_prologue* prologue
        = (mapped_read_location_prologue*) loc;
    return (MRL_CHR_TYPE) prologue->chr;
}

/** FRAGMENT/READ COVERAGE **/

static inline MRL_START_POS_TYPE
get_start_from_mapped_read_location( const mapped_read_location* loc )
{
    /* The index of the first sublocation is the start of the mapped read
     * location */
    mapped_read_sublocation* first_subloc
        = get_start_of_sublocations_in_mapped_read_location( loc );

    MRL_START_POS_TYPE start = first_subloc->start_pos;

    assert( start <= LOCATION_MAX );

    return start;
}

static inline MRL_STOP_POS_TYPE
get_stop_from_mapped_read_location( const mapped_read_location* loc )
{
    mapped_read_sublocation* sl_start
        = get_start_of_sublocations_in_mapped_read_location( loc );

    MRL_STOP_POS_TYPE stop = get_stop_for_sublocations_group( sl_start );
    
    assert( stop <= LOCATION_MAX );

    return stop;
}

static inline MRL_FL_TYPE
get_fl_from_mapped_read_location( const mapped_read_location* loc )
{
    MRL_START_POS_TYPE start = get_start_from_mapped_read_location( loc );
    MRL_STOP_POS_TYPE stop = get_stop_from_mapped_read_location( loc );
    
    assert( stop - start <= FRAGMENT_LENGTH_MAX );

    return (MRL_FL_TYPE) stop - start;
}

/** SEQ_ERROR **/

static inline ML_PRB_TYPE
get_log_seq_error_from_mapped_read_location( const mapped_read_location* loc )
{
    mapped_read_location_prologue* prologue
        = (mapped_read_location_prologue*) loc;
    return (ML_PRB_TYPE) prologue->log_seq_error;
}

static inline ML_PRB_TYPE
get_seq_error_from_mapped_read_location( const mapped_read_location* loc )
{
    return pow(10, get_log_seq_error_from_mapped_read_location(loc));
}

static inline enum bool
mapped_read_location_is_paired( const mapped_read_location* loc )
{
    /* if one sublocation has next_subread_is_ungapped set, then the
     * mapped_read_location is from a paired end read */

    mapped_read_sublocation* sl_start
        = get_start_of_sublocations_in_mapped_read_location( loc );

    enum bool is_paired = false;

    char* ptr = (char*) sl_start;

    while( !last_sublocation_in_mapped_read_location(
                (mapped_read_sublocation*) ptr) )
    {
        mapped_read_sublocation* current_sl
            = (mapped_read_sublocation *) ptr;

        if( current_sl->next_subread_is_ungapped )
        {
            is_paired = true;
            break;
        }

        ptr += sizeof( mapped_read_sublocation );
    }

    return is_paired;
}

/* TODO - make sure we're setting rev_comp on the mapped_read_sublocation */
static inline enum bool
first_read_in_mapped_read_location_is_rev_comp(
        const mapped_read_location* loc )
{
    mapped_read_sublocation* sl_start
        = get_start_of_sublocations_in_mapped_read_location( loc );

    /* if the first sublocation is reverse complemented, then the first read is
     * reverse complement (whether it was gapped or not - strand matches for all
     * sublocations in a gapped read) */
    if( sl_start->rev_comp == 0 ) {
        return false;
    }

    return true;
}

void
init_mapped_read( mapped_read_t** rd );

void
free_mapped_read( mapped_read_t* rd );

void
add_location_to_mapped_read( mapped_read_location* loc,
                             mapped_read_t** rd );

MPD_RD_ID_T
get_read_id_from_mapped_read( mapped_read_t* rd );

void
fprintf_mapped_read( FILE* fp, mapped_read_t* r );

void
advance_pointer_to_start_of_next_joined_candidate_mappings(
        candidate_mapping** current_mapping,
        int current_set_index,
        int num_sets );

void
join_candidate_mappings( candidate_mappings* mappings, 
                         struct read *r,
                         candidate_mapping** joined_mappings, 
                         float* penalties,
                         int* joined_mappings_len );

int
filter_joined_candidate_mappings( candidate_mapping** joined_mappings, 
                                  float* penalties,
                                  int joined_mappings_len,
                                  
                                  struct genome_data* genome,
                                  struct read* r,
                                  struct fragment_length_dist_t* fl_dist,

                                  struct mapping_params* mapping_params
    );

mapped_read_t*
build_mapped_read_from_joined_candidate_mappings(
        MPD_RD_ID_T read_id,
        struct joined_candidate_mappings* joined_mappings,
        int joined_mappings_len );

size_t 
write_mapped_read_to_file( mapped_read_t* read, FILE* of  );

/***** mapped read index *****/

typedef struct {
    mapped_read_t* rd;
    MPD_RD_ID_T read_id;
    MPD_RD_ID_T num_mappings;
    mapped_read_location** mappings;
} mapped_read_index;

void
init_mapped_read_index( mapped_read_index** index,
                        mapped_read_t* rd );

void
free_mapped_read_index( mapped_read_index* index );

/*
 *  END Mapped Reads
 *
 **************************************************************************/


/*****************************************************************************
 *
 * Mapped Reads DB Code
 *
 */

struct mapped_reads_db_index_t {
    MPD_RD_ID_T read_id;
    char* ptr;
};

struct mapped_reads_db {

    /* A read may map, be considered unmappable, or be considered mappable but
     * fail to map (nonmapping). These files store data about each type of
     * read. */
    FILE* mapped_fp;
    MPD_RD_ID_T num_mapped_reads;
    pthread_mutex_t* mapped_mutex;

    FILE* unmappable_fp;
    MPD_RD_ID_T num_unmappable_reads;
    pthread_mutex_t* unmappable_mutex;

    FILE* nonmapping_fp;
    MPD_RD_ID_T num_nonmapping_reads;
    pthread_mutex_t* nonmapping_mutex;

    char mode; // 'r' or 'w'

    /* mmap data */
    /* pointer to the mmapped data and its size in bytes */
    char* mmapped_data;
    size_t mmapped_data_size; 
    
    /* mmap index */
    struct mapped_reads_db_index_t* index;
    
    /* store the number of times that each read has been
       iterated over, and found to be below the update threshold */
    char* num_succ_iterations;
    
    struct fragment_length_dist_t* fl_dist;

    MPD_RD_ID_T current_read;
};

void
open_mapped_reads_db_for_reading( struct mapped_reads_db** rdb, char* fname );

void
open_mapped_reads_db_for_writing( struct mapped_reads_db** rdb, char* fname );

void
close_mapped_reads_db( struct mapped_reads_db** rdb );

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    mapped_read_t* rd );

void
add_unmappable_read_to_mapped_reads_db(
        struct read* r,
        int error,
        struct mapped_reads_db* db );

void
add_nonmapping_read_to_mapped_reads_db(
        struct read* r,
        int error,
        struct mapped_reads_db* db );

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb );

/* returns EOF, and sets rd to NULL if we reach the end of the file */
int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* rdb, 
    mapped_read_t** rd
);

void
reset_read_cond_probs( struct cond_prbs_db_t* cond_prbs_db,
                       mapped_read_t* rd,
                       struct mapped_reads_db* mpd_rds_db );
                       
void
reset_all_read_cond_probs( struct mapped_reads_db* rdb,
                           struct cond_prbs_db_t* cond_prbs_db );

/*
 * this requires code from iterative mapping to write out the
 * traces. It also requires the mapped_reads to be mmapped and 
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

void
print_mapped_reads_db_index(
        struct mapped_reads_db* rdb
    );
/*
 *  END Mapped Reads DB
 *
 **************************************************************************/

/*****************************************************************************
 *
 * Mapped Reads Conditional Probabilities Code
 *
 */


struct cond_prbs_db_t {
    /* this should be the same length as the number of reads in the mapped read db */
    /* each pointer should contain num_mapped_location entries */
    ML_PRB_TYPE** cond_read_prbs;
    MPD_RD_ID_T max_rd_id;
};

void
init_cond_prbs_db_from_mpd_rdb( 
    struct cond_prbs_db_t** cond_prbs_db,
    struct mapped_reads_db* mpd_rdb
);

void
free_cond_prbs_db( struct cond_prbs_db_t* cond_prbs_db );

static inline void
set_cond_prb( 
    struct cond_prbs_db_t* cond_prbs_db,
    MPD_RD_ID_T mpd_rd_id,
    int mapping_index,
    float prb
)          
{
    cond_prbs_db->cond_read_prbs[ mpd_rd_id ][ mapping_index ] 
        = ML_PRB_TYPE_from_float( prb );
    return;
}

static inline float
get_cond_prb( 
    struct cond_prbs_db_t* cond_prbs_db,
    MPD_RD_ID_T mpd_rd_id,
    int mapping_index
)          
{
    return ML_PRB_TYPE_from_float( 
        cond_prbs_db->cond_read_prbs[ mpd_rd_id ][ mapping_index ] );
}

#endif // MAPPED_READ
