/* Copyright (c) 2009-2010, Nathan Boley */

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
            return 'C':
        case 2:
            return 'G';
        case 3:
            return 'T';
    }

    fprintf(stderr, "PANIC - Error converting bp code %i in "
                    "penalty array build.\n"
        );
    exit( -1 );
}


/**** Penalty array functions ****/

void
init_penalty_array( int len, struct penalty_array_t* pa )
{
    /* Set dimensions of array */
    pa->len = len;

    /* Allocate memory for array of len penalty_t structs */
    pa->penalties = malloc( len * sizeof(penalty_t) );
}

void
free_penalty_array( struct penalty_array_t* pa )
{
    if( pa == NULL ) return;

    free( pa->penalties );
}

/**** penalty array builders ****/
void
build_error_data_bootstrap_penalty_array_from_rawread(
        struct rawread* rd,
        struct error_data_t* error_data,
        struct penalty_array_t* pa,
    )
{
    int i, j, k;
    /* for each position in the rawread */
    for( i = 0; i < pa->len; i++ )
    {
        /* for each possible basepair (A,C,G,T) in the reference sequence */
        for( j = 0; j < 4; j++ )
        {
            /* for each possible basepair (A,C,G,T) in the observed sequence */
            for( k = 0; k < 4; k++ )
            {
                if( j == k ) // if bp's match
                    pa->penalties[i].array[j][k] = 0;
                else
                    /*
                       we only want perfect mappers - if the bp's mismatch,
                       return 1, which is an invalid penalty score and will
                       cause the search to terminate immediately
                    */
                    pa->penalties[i].array[j][k] = 1;
            }
        }
    }
}

void
build_error_data_penalty_array_from_rawread(
        struct rawread* rd,
        struct error_data_t* error_data,
        struct penalty_array_t* pa,
    )
{
    int i, j, k;
    /* for each position in the rawread */
    for( i = 0; i < pa->len; i++ )
    {
        /* for each possible basepair (A,C,G,T) in the reference sequence */
        for( j = 0; j < 4; j++ )
        {
            /* for each possible basepair (A,C,G,T) in the observed sequence */
            for( k = 0; k < 4; k++ )
            {
                /* estimate error probability based on observed sequence and
                   error data */
                pa->penalties[i].array[j][k] = est_error_prb(
                        code_bp(j),
                        code_bp(k),
                        rd->error_str[i],
                        i,
                        error_data
                    );

                /* add mismatch penalty */
                if( j != k )
                {
                    /* log10(1/3) = -0.4771213 */
                    pa->penalties[i].array[j][k] += -0.4771213;
                }
            }
        }
    }
}

void
build_mismatch_penalty_array_from_rawread(
        struct rawread* rd,
        struct error_data_t* error_data,
        struct penalty_array_t* pa,
    )
{
    int i, j, k;
    /* for each position in the rawread */
    for( i = 0; i < pa->len; i++ )
    {
        /* for each possible basepair (A,C,G,T) in the reference sequence */
        for( j = 0; j < 4; j++ )
        {
            /* for each possible basepair (A,C,G,T) in the observed sequence */
            for( k = 0; k < 4; k++ )
            {
                if( j == k ) // if bp's match
                    pa->penalties[i].array[j][k] = 0;
                else
                    /*
                       each mismatch is treated as -1 penalty score, so we
                       can intutively control the number of mismatches allowed
                       with existing parameters
                    */
                    pa->penalties[i].array[j][k] = -1;
            }
        }
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
error_prb(
        char ref,
        char obs,
        char error_score,
        int pos,
        struct error_data_t* error_data
    )
{
    /*
       compute error score from observed error data
     */

    double loc_component = ( error_data->position_mismatch_cnts[pos] /
                             error_data->num_unique_reads );

    double mismatch_component;
    // avoid divide-by-0
    if( error_data->qual_score_cnts[(unsigned char) error_score] == 0 )
    {
        mismatch_component = 0;
    }
    else {
        mismatch_component =
            error_data->qual_score_mismatch_cnts[(unsigned char) error_score] /
            error_data->qual_score_cnts[(unsigned char) error_score];
    }

    float prb = ( loc_component + mismatch_component ) / 2 + LOG_FFACTOR;

    assert( prb > 0 && prb <= 1 ); // make sure prb is, in fact, a probability

    /* convert to log prb */
    float log_prb;
    /* if the bp's match, we return the inverse of the probability of error */
    if( ref == obs )
        log_prb = log10(1 - prb);
    else
        log_prb = log10(prb);

    return log_prb;
}

float 
multiple_letter_penalty(
        const LETTER_TYPE* const reference,
        const LETTER_TYPE* const observed,

        const int start_position,
        const int seq_length,
        const int num_letters,
        const float min_penalty,

        struct penalty_array_t* pa
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
        struct penalty_array_t* pa
    )
{
    int i;
    float penalty = 0;

    for( i = 0; i < LETTER_LEN; i++ )
    {
        // NOTE - penalty_t float array and LETTER_TYPE must use same encoding
        penalty += pa->array[LETTER_LEN*position+i].penalties[(ref&3)][(obs&3)];

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
        reference = reference >> 2;
        observed = observed >> 2;
    }

    return penalty;
}
