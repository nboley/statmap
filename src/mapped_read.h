/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef MAPPED_READ
#define MAPPED_READ

#define DONT_MALLOC_READS_DB

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "config.h"
#include "candidate_mapping.h"
#include "rawread.h"

typedef int MPD_RD_ID_T;

struct fragment_length_dist_t;
struct genome_data;
struct rawread_db_t;
struct cond_prbs_db_t;

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

/* TODO - make this a bitfield */

#define MRL_FLAG_TYPE unsigned char
#define MRL_CHR_TYPE unsigned short
#define MRL_START_POS_TYPE unsigned int
#define MRL_STOP_POS_TYPE unsigned int
#define MRL_FL_TYPE unsigned short

/* Set if the read that contribute to this mapped location are paired */
#define IS_PAIRED 1
/* Set if the first read ( ie read_name/1 ) in the pair was rev 
   complemented to map */
#define FIRST_READ_WAS_REV_COMPLEMENTED 2
/* Set if the first read ( ie read_name/1 ) maps to start_pos */
#define FIRST_PAIR_IS_FIRST_IN_GENOME 4

/*
 * A full mapped read.
 *
 */

/*
mapped_read_t {

    size_t num_allocated_bytes;

    // 4 bytes
    struct {
        read_id :31;
        are_more :1;
    } read_id_node;
    // potentially more read ids

    struct mapped_read_location {

        size_t num_allocated_bytes;

        // 2 bytes
        unsigned chr :CHR_BITS; // 15 bits
        unsigned flag :FLAG_BITS; // 8 bits TODO we may be able to get rid of this
        ML_PRB_TYPE seq_error;  // 4 bytes
        unsigned are_more :1    // 1 bit

        // 6 bytes
        struct {
            signed start_pos :LOCATION_BITS; ( 29 )
            signed length    :SUBTEMPLATE_LENGTH_BITS; ( 16 )
            signed strand    :1;
            signed are_more  :1;
        }
        // potentially more read subtemplates
    }

    // potentially more mapped read locations

} __attribute__((__packed__));
*/

typedef void mapped_read_location;

#define SUBTEMPLATE_LENGTH_BITS 16
#define MAX_SUBTEMPLATE_LENGTH 65535 // 2**16 - 1

#define FLAG_BITS 8
#define MAX_FLAG_LENGTH 255 // 2**8 - 1

// 11 bytes
typedef struct {

    size_t num_allocated_bytes; // 4 bytes

    unsigned chr :CHR_BITS; // 15 bits
    unsigned flag :FLAG_BITS; // 8 bits
    ML_PRB_TYPE seq_error; // 4 bytes
    unsigned are_more :1;   // 1 bits

} mapped_read_location_prologue;

// 6 bytes
typedef struct {
    unsigned start_pos :LOCATION_BITS; // 29 bits
    unsigned length :SUBTEMPLATE_LENGTH_BITS; // 16 bits
    unsigned strand :1;
    unsigned are_more :1;
} mapped_read_subtemplate_location;

/* 
 * small, inline functions for setting and accessing the items
 * in mapped_read_location. These are how elements should be 
 * accessed - there is no guarantee that mapped_read_location 
 * will remain the same 
 *
 * These are defined in the header so that they can be inlined.
 *
 */

void
init_mapped_read_location( mapped_read_location** loc,
                           MRL_CHR_TYPE chr,
                           MRL_FLAG_TYPE flag,
                           ML_PRB_TYPE seq_error,
                           enum bool are_more );

void
free_mapped_read_location( mapped_read_location* loc );

/** FLAG **/

static inline MRL_FLAG_TYPE
get_flag_from_mapped_read_location( const mapped_read_location* loc )
{
    /* cast loc to mapped_read_location_prologue to access the flag */
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    return (MRL_FLAG_TYPE) prologue->flag;
}

static inline void
set_flag_in_mapped_read_location( mapped_read_location* loc, MRL_FLAG_TYPE flag )
{
    assert( flag < MAX_FLAG_LENGTH );
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    prologue->flag = flag;
}

/** CHR **/

static inline MRL_CHR_TYPE
get_chr_from_mapped_read_location( const mapped_read_location* loc)
{
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    return (MRL_CHR_TYPE) prologue->chr;
}

static inline void
set_chr_in_mapped_read_location( const mapped_read_location* loc, MRL_CHR_TYPE value )
{
    assert( value <= CHR_NUM_MAX );
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    prologue->chr = value;
}

/** FRAGMENT/READ COVERAGE **/

static mapped_read_subtemplate_location*
get_start_of_subtemplate_locations( const mapped_read_location* loc )
{
    /* Given a pointer to a mapped_read_location, return a pointer to the start
     * of the first read subtemplate mapped location in the
     * mapped_read_location */

    // skip the location prologue
    char* ptr = (char*) loc;
    ptr += sizeof(mapped_read_location_prologue);

    return (mapped_read_subtemplate_location*) ptr;
}

static inline MRL_START_POS_TYPE
get_start_from_mapped_read_location( const mapped_read_location* loc )
{
    /* TODO for now, just return the start location of the first subtemplate
     * location */
    mapped_read_subtemplate_location* st_loc
        = get_start_of_subtemplate_locations( loc );
    return (MRL_STOP_POS_TYPE) st_loc->start_pos;
}

static inline MRL_STOP_POS_TYPE
get_stop_from_mapped_read_location( const mapped_read_location* loc )
{
    /* TODO for now, just return the stop of the first subtemplate location */
    mapped_read_subtemplate_location* st_loc
        = get_start_of_subtemplate_locations( loc );
    return (MRL_STOP_POS_TYPE) st_loc->start_pos + st_loc->length;
}

static inline MRL_FL_TYPE
get_fl_from_mapped_read_location( const mapped_read_location* loc )
{
    /* TODO for now, just return the stop of the first subtemplate location */
    mapped_read_subtemplate_location* st_loc
        = get_start_of_subtemplate_locations( loc );
    return (MRL_FLAG_TYPE) st_loc->length;
}

static inline void
set_start_and_stop_in_mapped_read_location(
        mapped_read_location* loc,
        int start,
        int stop )
{
    // TODO for now, just set start and stop of the first subtemplate location
    mapped_read_subtemplate_location* st_loc
        = get_start_of_subtemplate_locations( loc );

    assert( start >= 0 );
    assert( start < LOCATION_MAX );
    st_loc->start_pos = start;

    assert( stop >= 0 );
    assert( stop < LOCATION_MAX );

    assert( (stop - start) > FRAGMENT_LENGTH_MIN );
    assert( (stop - start) < FRAGMENT_LENGTH_MAX );

    st_loc->length = stop - start;
}

static inline ML_PRB_TYPE
get_seq_error_from_mapped_read_location( const mapped_read_location* loc )
{
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    return (ML_PRB_TYPE) prologue->seq_error;
}

static inline void
set_seq_error_in_mapped_read_location( 
        const mapped_read_location* loc,
        ML_PRB_TYPE value )
{
    mapped_read_location_prologue* prologue = 
        (mapped_read_location_prologue*) loc;
    prologue->seq_error = value;
}

typedef void mapped_read_t;

/* Note - the following structures must all be aligned on byte boundaries so we
 * can use sizeof and pointer arithmetic to work with the mapped_read_t pseudo
 * structure */

/* Structures used to access and manipulate the mapped_read_t and
 * mapped_read_location pseudo structs */

#define MAX_READ_ID 2147483648 // 2**31 = 2147483648

// 4 bytes
typedef struct {
    unsigned read_id :31;
    unsigned are_more :1;
} read_id_node;

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
build_mapped_read_from_candidate_mappings( 
        mapped_read_t** mpd_rd,
        candidate_mappings* mappings, 
        MPD_RD_ID_T read_id );

size_t 
write_mapped_read_to_file( mapped_read_t* read, FILE* of  );

/***** mapped read index *****/

typedef struct {
    mapped_read_t* rd;
    MPD_RD_ID_T read_id;
    MPD_RD_ID_T num_mappings; // TODO change to num_mappings
    mapped_read_location** mappings; // TODO change to mappings
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
    FILE* fp;

    char mode; // 'r' or 'w'
    MPD_RD_ID_T num_mapped_reads;

    pthread_mutex_t* mutex;

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
build_fl_dist_from_file( struct mapped_reads_db* rdb, FILE* fl_fp );

void
build_fl_dist_from_filename( struct mapped_reads_db* rdb, char* filename );

void
close_mapped_reads_db( struct mapped_reads_db** rdb );

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    mapped_read_t* rd
);

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
