/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "quality.h"
#include "rawread.h"
#include "error_correction.h"

/*
 *  This determines whether we consider the quality scores to 
 *  be log odds ( as in Solexa v <= 1.3 ) or log probabilites 
 *  ( pretty much everything else ). This parameter will be set 
 *  by the argument parsing function. 
 *
 */
enum bool ARE_LOG_ODDS = true;

/*
 * This determines how we shift ascii quality scores
 * to get the actual numeric quality score. This parameter
 * must be set in the argument parsing function. 
 *
 *
 *
 */
int QUAL_SHIFT = -1;

/* Utility functions to convert bps (char) to their encoded numerical reprs */
int 
bp_code( const char letter )
{
    switch( letter )
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    }

    fprintf(stderr, "PANIC - Error converting '%c' in recheck\n", letter );
    assert( false );
    exit( -1 );
}

char
code_bp( int code )
{
    switch( code )
    {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
    }

    fprintf(stderr, "PANIC - Error converting bp code %i in "
                    "penalty array build.\n", code
        );
    exit( -1 );
}

/**** penalty functions ****/

float
error_prb_for_mismatch( char ref, char obs )
{
    if( ref == obs )
        return 0;
    
    return -1;
}


float
error_prb_for_estimated_model(
        char ref,
        char obs,
        char error_score,
        int pos,
        struct error_model_t* error_model
    )
{
    struct freqs_array* freqs = error_model->data;
    double prb = freqs->freqs[(unsigned char)error_score][pos];

    if( ref == obs )
        return log10(1 - prb);
    else
        return log10(prb);
    
    assert( false );
}

float
error_prb(
        char ref,
        char obs,
        char error_score,
        int pos,
        struct error_model_t* error_model
    )
{
    switch ( error_model->error_model_type ) 
    {
    case MISMATCH:
        return error_prb_for_mismatch( ref, obs );
    case FASTQ_MODEL:
        assert( false );
        return 0;
    case ESTIMATED:
        return error_prb_for_estimated_model( 
            ref, obs, error_score, pos, error_model );
    default:
        assert( false );
    }
}


/**** Penalty array functions ****/

void
init_penalty_array( struct penalty_array_t* pa, int length )
{
    /* Set dimensions of array */
    pa->len = length;

    /* Allocate memory for array of len penalty_t structs */
    pa->array = malloc( length * sizeof(struct penalty_t) );
}

void
free_penalty_array( struct penalty_array_t* pa )
{
    if( pa == NULL ) return;

    free( pa->array );
}

void
fprintf_penalty_array( FILE* fp, struct penalty_array_t* pa )
{
    int i;
    for( i = 0; i < pa->len; i++ )
    {
        fprintf(fp, "%i: ", i); // index
        int x, y;
        for( x = 0; x < 4; x++ )
        {
            for( y = 0; y < 4; y++ )
            {
                fprintf(fp, "%f ", pa->array[i].penalties[x][y]);
            }
            if( x < 3 )
                fprintf(fp, "| ");
        }
        fprintf(fp, "\n");
    }
}

/**** penalty array builders ****/
void
build_penalty_array(
        struct penalty_array_t* pa,
        struct error_model_t* error_model,
        char* error_str
    )
{
    int pos, ref_bp, seq_bp;
    /* for each position in the rawread */
    for( pos = 0; pos < pa->len; pos++ )
    {
        /* for each possible basepair (A,C,G,T) in the reference sequence */
        for( ref_bp = 0; ref_bp < 4; ref_bp++ )
        {
            /* for each possible basepair (A,C,G,T) in the observed sequence */
            for( seq_bp = 0; seq_bp < 4; seq_bp++ )
            {
                /* estimate error probability based on observed sequence and
                   error data */
                pa->array[pos].penalties[ref_bp][seq_bp] = error_prb(
                        code_bp(ref_bp),
                        code_bp(seq_bp),
                        error_str[pos],
                        pos,
                        error_model
                    );
            }
        }
    }
}

/**** build reverse penalty array ****/
/* simply reverses the forward penalty array */
void
build_reverse_penalty_array(
        struct penalty_array_t* fwd_pa,
        struct penalty_array_t* rev_pa
    )
{
    // the underlying penalty arrays must be the same length
    assert( fwd_pa->len == rev_pa->len );

    int i;
    for( i = 0; i < fwd_pa->len; i++ )
    {
        rev_pa->array[rev_pa->len-1-i] = fwd_pa->array[i];
    }
}

double
mutation_probability( int qscore )
{
    double I = pow(10, ((((double) qscore)-QUAL_SHIFT)/-10.0) );
    if( ARE_LOG_ODDS == true ) {
        return I/( 1 + I );
    } else {
        return I;
    }
}

unsigned char
quality_score( float log10_p )
{
    double p = pow(10, log10_p);
    if( ARE_LOG_ODDS == true ) {
        p = p/(1-p);
    }
        
    p = log10(p);
    p = -10*p;
    
    int c  = (int) floor(p + 0.5);
    
    return (unsigned char) ( c + QUAL_SHIFT );
}

void
convert_into_quality_string( float* mutation_probs, char* quality, int seq_len )
{
    int i;
    for( i = 0; i < seq_len; i++)
    {
        quality[i] = quality_score( (double) mutation_probs[i] );
    }
    quality[seq_len] = '\0';
    return;
}

float 
multiple_letter_penalty(
        const LETTER_TYPE* const reference,
        const LETTER_TYPE* const observed,

        const int start_position,
        const int seq_length,
        const int num_letters,
        const float min_penalty,

        struct penalty_t* pa
    )
{
    int j;
    float cum_penalty = 0;

    for( j = 0; j < num_letters - start_position; j++ )
    {
        /* determine the penalty contribution from letter j in seq i */ 
        float added_penalty = compute_penalty(
                reference[j], observed[j],
                start_position + j,
                seq_length,
                min_penalty - cum_penalty,
                pa
            );
        
        /* 
         * if the penalty is > 0, then that indicates that it exceeded the
         * minimum at the (curr_penalty - 1)'th basepair in letter j. See
         * the penalty function for details.
         */
        if( added_penalty > 0.5 )  {
            /* skip this sequence */
            return j*LETTER_LEN + added_penalty;
        } else {
            /* update the current penalty */
            cum_penalty += added_penalty;
        }
    }

    return cum_penalty;
}

float
recheck_penalty(
        char* reference,
        char* observed,
        const int seq_length,

        struct penalty_array_t* pa
    )
{
    int i;
    float penalty = 0;

    for( i = 0; i < seq_length; i++ )
    {
        char ref = toupper( reference[i] );
        char obs = toupper( observed[i] );

        /* if it's an N, put in a max penalty substitution */
        if( ref == 'N' || obs == 'N' )
        {
            penalty += N_penalty;
        } else {
            penalty += pa->array[i].penalties[bp_code(ref)][bp_code(obs)];
        }
    }

    return penalty;
}

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
        struct penalty_t* pa
    )
{
    int i;
    float penalty = 0;

    for( i = 0; i < LETTER_LEN; i++ )
    {
        /*
           make sure we haven't run off the string
           (happens if seq_length is not divisible by LETTER_LEN)
        */
        if( LETTER_LEN*position + i >= seq_length ) break;

        /* 
           add penalty
           NOTE - penalty_t float array and LETTER_TYPE must use same encoding
        */
        penalty += pa[LETTER_LEN*position+i].penalties[(ref&3)][(obs&3)];  

        /* 
         * If we have surpassed the minimum allowed penalty, there is no 
         * reason to continue. We return the current bp index as a hint for 
         * the calling function - we may be able to skip basepairs with the 
         * same minor bits. Since the penalty *must* be negative, we can 
         * always interpret a positive number as a failure.
         */
        if( penalty < min_penalty ) {
            return i+1;
        }

        /* Shift the bits down to consider the next basepair */
        ref = ref >> 2;
        obs = obs >> 2;
    }

    return penalty;
}
