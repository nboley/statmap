/* Copyright (c) 2009-2010 Nathan Boley */

#ifndef DNA_SEQUENCE
#define DNA_SEQUENCE

#include "config.h"

/*
 * packed_sequence.h
 *
 * Functions and data types for manipulating packed sequences.
 *
 */

inline int calc_num_letters( const int seq_len );

LETTER_TYPE*
copy_first_k_basepairs( 
    LETTER_TYPE* read, LETTER_TYPE** new_read, 
    int read_len, int num_letters, int k
);

LETTER_TYPE*
copy_last_k_basepairs( 
    LETTER_TYPE* read, LETTER_TYPE** new_read, 
    int read_len, int k
);

/* returns 0 on success, 5 on error */
int
complement_read( char* read, char* mut_read, int seq_len );

int
rev_complement_read( char* read, char* mut_read, int seq_len );

/* returns 0 on success, 5 on error */
int
reverse_read( char* read, char* mut_read, int seq_len );

void convert_packed_sequence( LETTER_TYPE*, char*, int );

void
print_packed_sequence( LETTER_TYPE* seq, int );

/* This allocates memory for the new sequence */
inline LETTER_TYPE* 
translate_seq(char*, int num_letters, LETTER_TYPE** );

inline int 
cmp_letters( LETTER_TYPE, LETTER_TYPE );

inline int 
cmp_words( LETTER_TYPE*, LETTER_TYPE*, int );

#endif // define DNA_SEQUENCE
