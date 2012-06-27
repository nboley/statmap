/* Copyright (c) 2009-2010 Nathan Boley */

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

extern enum bool ARE_LOG_ODDS;
extern int QUAL_SHIFT;

/*
   3d array of penalties:
   length of read * ref bp * obs bp
 */
struct penalty_array {
    /* Array dimensions */
    int x, y, z;

    float*** array;
};

void
init_penalty_array( int x, int y, int z, struct penalty_array* pa );

void
free_penalty_array( struct penalty_array* pa );

// const int solexa_qual_start = 64;
// const int num_possible_quality_chars = 85;

/* convert a solexa quality score into a mutation probability */
double mutation_probability( int );

/* convert a mutation probability into a solexa quality character */
unsigned char quality_score( float );

/* 
 * populate a char array with solexa quality scores from an array of mutation
 * probabilities.
 */
void convert_into_quality_string( float* mutation_probs, char* quality, int seq_len );

/* FIXME - change the name */
extern void determine_bp_mut_rates( float** );
extern void determine_bp_mut_rates_for_bootstrap( float** );
extern void determine_bp_mut_rates_for_mismatches( float** );

/* fwd declarations */
struct rawread;
struct error_data_t;

/* populate float arrays with mutation probabilities */
void
build_lookup_table_from_rawread( struct rawread* rd,
                                 struct error_data_t* error_data,
                                 float* lookuptable_position,
                                 float* inverse_lookuptable_position,
                                 float* reverse_lookuptable_position,
                                 float* reverse_inverse_lookuptable_position
);

void
build_mismatch_lookup_table( float** lookuptable_position,
                             float** inverse_lookuptable_position,
                             float** lookuptable_bp,
                             int seq_len
);

/* print the float arrays */
void print_lookup_table( float*, float*, int seq_len );

float
recheck_penalty( char* reference, 
                 char* observed, 
                 const int seq_length,
                 const float* n_lookuptable_position,
                 const float* n_inverse_lookuptable_position,
                 const float* n_lookuptable_bp );


/* calculate the penalty of a packed sequence with respect to another */
float 
multiple_letter_penalty( const LETTER_TYPE* reference,
                         const LETTER_TYPE* observed,
                         const int start_position,
                         const int seq_length,
                         const int num_letters,
                         const float min_penalty,
                         const float* lookuptable_position,
                         const float* inverse_lookuptable_position,
                         const float* lookuptable_bp
);

float 
penalty_func( LETTER_TYPE reference, 
              LETTER_TYPE observed, 
              /* the position in the sequence - this should be
                 zero indexed */
              const int position, 
              const int seq_len,
              const float min_penalty,
              const float* lookuptable_position,
              const float* inverse_lookuptable_position,
              const float* lookuptable_bp 
);

float
est_error_prb( char bp, char error_score, enum bool inverse, 
               int pos, struct error_data_t* error_data );

#endif /* #define QUALITY */
