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
//#include "read.h"

struct read_subtemplate; // fwd declaration (?)

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
    struct error_model_t* error_model,
    int read_subtemplate_index,
    enum STRAND strand
);

/* 
   the probability of observing obs given that the ref basepair is unobservable.
   this is just log10( 0.25 ) - the N gives us 0 information. 
*/
static const double N_penalty = -0.60206;

/* Stores the penalty data for a single character at a specific position
   in a read. This is the penalty addition to make for each possible bp if it
   matches to the given bp at this position in the read. */
struct penalty_t {
    float penalties[4];
};

struct penalty_array_t {
    int length;
    struct penalty_t* array;
};

void
init_penalty_array( struct penalty_array_t* pa, int length );

void
free_penalty_array( struct penalty_array_t* pa );

void
build_penalty_array(
        struct penalty_array_t* pa,
        struct read_subtemplate* rst,
        struct error_model_t* error_model
    );

void
build_reverse_penalty_array(
        struct penalty_array_t* rev_pa,
        struct read_subtemplate* rst,
        struct error_model_t* error_model
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
struct read;
struct error_data_t;


/* Since log10(0) is undefined, we add a tiny fudge factor to
   the error probabilities (in case p=0) */
#define LOG_FFACTOR 1e-6

/* calculate the penalty of a packed sequence with respect to another */
float 
multiple_letter_penalty(
        const LETTER_TYPE* const reference,

        const int start_position,
        
        const int num_letters,
        const float min_penalty,

        struct penalty_t* pa
    );

/* Compute the penalty from char sequences */
float
recheck_penalty(
        char* reference,
        /* Pointer into an array of penalty_t */
        struct penalty_t* pa,
        const int seq_length
    );

/* Compute the penalty form LETTER_TYPE sequences */
float
compute_penalty(
        /* 
         * not const because it is a copy. I bit shift
         * them while calculating the penalty.
         */
        LETTER_TYPE ref,

        /* the position in the sequence - this should be
           zero indexed */
        const int position, 
        
        /* the minimum allowable penalty */
        const float min_penalty,

        /* the penalty array */
        struct penalty_t* pa
    );

#endif /* #define QUALITY */
