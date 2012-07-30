/* Copyright (c) 2009-2012 Nathan Boley */

#ifndef QUALITY
#define QUALITY

/*
 * quality.h
 *
 * Functions and data types for converting fastq quality strings into
 * mutation probability arrays, and vice versa.
 *
 */

#include "dna_sequence.h"
#include "rawread.h"

extern enum bool ARE_LOG_ODDS;
extern int QUAL_SHIFT;

int 
bp_code( const char letter );

char
code_bp( int code );

struct error_model_t;

float
error_prb(
    char ref,
    char obs,
    char error_score,
    int pos,
    struct error_model_t* error_model
);

/* 
   the probability of observing obs given that the ref basepair is unobservable. 
   this is just log10( 0.25 ) - the N gives us 0 information. 
*/
static const double N_penalty = -0.60206;

/* Stores the penalty data for a single character at a specific position
   in a read */
struct penalty_t {
    /* 4x4 array of penalties - ref bp x obs bp */
    float penalties[4][4];
};

/*
   3d array of penalties:
   length of read * ref bp * obs bp
 */
struct penalty_array_t {
    /* length of read */
    int len;

    /* array of penalty structs, for each position in the read */
    struct penalty_t* array;
};

void
init_penalty_array( int len, struct penalty_array_t* pa );

void
free_penalty_array( struct penalty_array_t* pa );

void
build_penalty_array(
        struct rawread* rd,
        struct error_model_t* error_model,
        struct penalty_array_t* pa
);

void
build_reverse_penalty_array(
        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa
    );

/* convert a solexa quality score into a mutation probability */
double mutation_probability( int );

/* convert a mutation probability into a solexa quality character */
unsigned char quality_score( float );

/* 
 * populate a char array with solexa quality scores from an array of mutation
 * probabilities.
 */
void convert_into_quality_string( float* mutation_probs, char* quality, int seq_len );

/* fwd declarations */
struct rawread;
struct error_data_t;


/* Since log10(0) is undefined, we add a tiny fudge factor to
   the error probabilities (in case p=0) */
#define LOG_FFACTOR 1e-6

/* calculate the penalty of a packed sequence with respect to another */
float 
multiple_letter_penalty(
        const LETTER_TYPE* const reference,
        const LETTER_TYPE* const observed,

        const int start_position,
        const int seq_length,
        const int num_letters,
        const float min_penalty,

        struct penalty_array_t* pa
    );

/* Compute the penalty from char sequences */
float
recheck_penalty(
        char* reference,
        char* observed,
        const int seq_length,

        struct penalty_array_t* pa
    );

/* Compute the penalty form LETTER_TYPE sequences */
float
compute_penalty(
        /* 
         * these aren't const because they are copies. I bit shift
         * them while calculating the penalty.
         */
        LETTER_TYPE ref,
        LETTER_TYPE obs,

        /* the position in the sequence - this should be
           zero indexed */
        const int position, 
        /* the length of a full sequence */
        const int seq_length,
        /* the minimum allowable penalty */
        const float min_penalty,

        /* the penalty array */
        struct penalty_array_t* pa
    );

#endif /* #define QUALITY */
