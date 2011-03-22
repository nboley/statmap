/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "quality.h"
#include "rawread.h"

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


/*
 * quality.c
 *
 * Functions and data types for converting fastq quality strings into
 * mutation probability arrays, and vice versa.
 *
 */

/* 
 * Quality Lookup Tables
 * 
 * This converts the integer quality score returned by
 * solexa sequencing machines from the GA pipeline v <1.3.
 *
 * The sanger standard is that the Quality score is
 * Q_sanger = -10*log_10( p ), where p is the probability 
 * of a mutation. 
 *
 * Solexa calculates the quality as
 * Q_solexa = -10*log_10( p / 1-p )
 *
 * ( see http://en.wikipedia.org/wiki/FASTQ_format#Variations )
 *
 * The table, solexa_lookuptable_score stores the value of p
 * for quality scores ranging from -5 to 62. That is, there are
 * 68 entries total. 
 *
 *
 * here is the python code that generated these table:

from math import log10
npc = 4 # number per column

# generate the logodd q's
qs = [ pow(10, s/-10.0)/(1+pow(10, s/-10.0)) \
     for s in range(-5, 63) ]



# generate the log-odds mutation table
for i in range(0, len(qs), npc ):
    print ",\t".join( "%.3e" % log10(qs[j])
                  for j in xrange(i, i+npc) 
                            if j < len(qs) ) + ","

# generate the log-odds inverse mutation table
for i in range(0, len(qs), npc ):
    print ",\t".join( "%.3e" % log10(1-qs[j])
                  for j in xrange(i, i+npc) 
                            if j < len(qs) ) + ","

# generate the logprb q's
qs = [ min( 1-1e-15, pow(10, s/-10.0) ) for s in range(0, 63) ]
     
# generate the log-odds mutation table
for i in range(0, len(qs), npc ):
    print ",\t".join( "%.3e" % log10(qs[j])
                  for j in xrange(i, i+npc) 
                            if j < len(qs) ) + ","

# generate the log-odds inverse mutation table
for i in range(0, len(qs), npc ):
    print ",\t".join( "%.3e" % log10(1-qs[j])
                  for j in xrange(i, i+npc) 
                            if j < len(qs) ) + ","

 *
 *
 */

static const double 
logodds_lookuptable_score[ 68 ] = {
    -1.193e-01, -1.455e-01, -1.764e-01, -2.124e-01,
    -2.539e-01, -3.010e-01, -3.539e-01, -4.124e-01,
    -4.764e-01, -5.455e-01, -6.193e-01, -6.973e-01,
    -7.790e-01, -8.639e-01, -9.515e-01, -1.041e+00,
    -1.133e+00, -1.227e+00, -1.321e+00, -1.417e+00,
    -1.514e+00, -1.611e+00, -1.709e+00, -1.807e+00,
    -1.905e+00, -2.004e+00, -2.103e+00, -2.203e+00,
    -2.302e+00, -2.402e+00, -2.501e+00, -2.601e+00,
    -2.701e+00, -2.801e+00, -2.901e+00, -3.000e+00,
    -3.100e+00, -3.200e+00, -3.300e+00, -3.400e+00,
    -3.500e+00, -3.600e+00, -3.700e+00, -3.800e+00,
    -3.900e+00, -4.000e+00, -4.100e+00, -4.200e+00,
    -4.300e+00, -4.400e+00, -4.500e+00, -4.600e+00,
    -4.700e+00, -4.800e+00, -4.900e+00, -5.000e+00,
    -5.100e+00, -5.200e+00, -5.300e+00, -5.400e+00,
    -5.500e+00, -5.600e+00, -5.700e+00, -5.800e+00,
    -5.900e+00, -6.000e+00, -6.100e+00, -6.200e+00
};

static const double 
logodds_inverse_lookuptable_score[ 68 ] = {
    -6.193e-01, -5.455e-01, -4.764e-01, -4.124e-01,
    -3.539e-01, -3.010e-01, -2.539e-01, -2.124e-01,
    -1.764e-01, -1.455e-01, -1.193e-01, -9.732e-02,
    -7.901e-02, -6.389e-02, -5.150e-02, -4.139e-02,
    -3.320e-02, -2.657e-02, -2.124e-02, -1.695e-02,
    -1.352e-02, -1.077e-02, -8.580e-03, -6.829e-03,
    -5.433e-03, -4.321e-03, -3.436e-03, -2.732e-03,
    -2.171e-03, -1.726e-03, -1.371e-03, -1.090e-03,
    -8.657e-04, -6.878e-04, -5.464e-04, -4.341e-04,
    -3.448e-04, -2.739e-04, -2.176e-04, -1.729e-04,
    -1.373e-04, -1.091e-04, -8.664e-05, -6.883e-05,
    -5.467e-05, -4.343e-05, -3.450e-05, -2.740e-05,
    -2.177e-05, -1.729e-05, -1.373e-05, -1.091e-05,
    -8.665e-06, -6.883e-06, -5.467e-06, -4.343e-06,
    -3.450e-06, -2.740e-06, -2.177e-06, -1.729e-06,
    -1.373e-06, -1.091e-06, -8.665e-07, -6.883e-07,
    -5.467e-07, -4.343e-07, -3.450e-07, -2.740e-07
};

static const double 
logprb_lookuptable_score[ 63 ] = {
     0.000e+00, -1.000e-01, -2.000e-01, -3.000e-01,
    -4.000e-01, -5.000e-01, -6.000e-01, -7.000e-01,
    -8.000e-01, -9.000e-01, -1.000e+00, -1.100e+00,
    -1.200e+00, -1.300e+00, -1.400e+00, -1.500e+00,
    -1.600e+00, -1.700e+00, -1.800e+00, -1.900e+00,
    -2.000e+00, -2.100e+00, -2.200e+00, -2.300e+00,
    -2.400e+00, -2.500e+00, -2.600e+00, -2.700e+00,
    -2.800e+00, -2.900e+00, -3.000e+00, -3.100e+00,
    -3.200e+00, -3.300e+00, -3.400e+00, -3.500e+00,
    -3.600e+00, -3.700e+00, -3.800e+00, -3.900e+00,
    -4.000e+00, -4.100e+00, -4.200e+00, -4.300e+00,
    -4.400e+00, -4.500e+00, -4.600e+00, -4.700e+00,
    -4.800e+00, -4.900e+00, -5.000e+00, -5.100e+00,
    -5.200e+00, -5.300e+00, -5.400e+00, -5.500e+00,
    -5.600e+00, -5.700e+00, -5.800e+00, -5.900e+00,
    -6.000e+00, -6.100e+00, -6.200e+00
};

static const double 
logprb_inverse_lookuptable_score[ 63 ] = {
    -1.000e-01, -6.868e-01, -4.329e-01, -3.021e-01,
    -2.205e-01, -1.651e-01, -1.256e-01, -9.665e-02,
    -7.494e-02, -5.844e-02, -4.576e-02, -3.594e-02,
    -2.830e-02, -2.233e-02, -1.764e-02, -1.396e-02,
    -1.105e-02, -8.753e-03, -6.938e-03, -5.502e-03,
    -4.365e-03, -3.463e-03, -2.749e-03, -2.182e-03,
    -1.732e-03, -1.376e-03, -1.092e-03, -8.674e-04,
    -6.889e-04, -5.471e-04, -4.345e-04, -3.451e-04,
    -2.741e-04, -2.177e-04, -1.729e-04, -1.374e-04,
    -1.091e-04, -8.666e-05, -6.884e-05, -5.468e-05,
    -4.343e-05, -3.450e-05, -2.740e-05, -2.177e-05,
    -1.729e-05, -1.373e-05, -1.091e-05, -8.665e-06,
    -6.883e-06, -5.467e-06, -4.343e-06, -3.450e-06,
    -2.740e-06, -2.177e-06, -1.729e-06, -1.373e-06,
    -1.091e-06, -8.665e-07, -6.883e-07, -5.467e-07,
    -4.343e-07, -3.450e-07, -2.740e-07,
};


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

void
build_lookup_table_from_rawread ( struct rawread* rd,
                                  float* lookuptable_position,
                                  float* inverse_lookuptable_position,
                                  float* reverse_lookuptable_position,
                                  float* reverse_inverse_lookuptable_position
    )
{
    int i;
    for( i = 0; i < rd->length; i++ )
    {
        unsigned char quality_char = ((unsigned char) rd->error_str[i]) - QUAL_SHIFT;

        /* check to see if the read is an 'N'. If it is, set the qual to the min */
        if( rd->char_seq[i] == 'N' || rd->char_seq[i] == 'n' )
        {
            /* set the probability that this is incorrect to 0.75 */
            lookuptable_position[i] = -0.1249387;
        } else {
            if( ARE_LOG_ODDS == false )
            {
                lookuptable_position[i] = logodds_lookuptable_score[ quality_char ];
            } else {
                lookuptable_position[i] = logprb_lookuptable_score[ quality_char ];
            }
        }
        /* set the reverse position */
        reverse_lookuptable_position[rd->length-1-i] = lookuptable_position[i];
        
        /* check to see if the read is an 'N'. If it is, set the qual to the min */
        if( rd->char_seq[i] == 'N' || rd->char_seq[i] == 'n' )
        {
            /* set the probability that this is incorrect to 0.75 */
            inverse_lookuptable_position[i] = -0.60206;
        } else {        
            if( ARE_LOG_ODDS == false ) {
                inverse_lookuptable_position[i] = 
                    logodds_inverse_lookuptable_score[ quality_char ];
            } else {
                inverse_lookuptable_position[i] = 
                    logprb_inverse_lookuptable_score[ quality_char ];
            }
        }
        reverse_inverse_lookuptable_position[rd->length-1-i] 
            = inverse_lookuptable_position[i];
    }

    return;
}

void
build_mismatch_lookup_table( float** lookuptable_position,
                             float** inverse_lookuptable_position,
                             float** lookuptable_bp,
                             int seq_len
    )
{
    *lookuptable_position = malloc( sizeof(float)*seq_len );
    *inverse_lookuptable_position = malloc( sizeof(float)*seq_len );
    *lookuptable_bp = malloc( sizeof(float)*16 );

    int i;
    for( i = 0; i < seq_len; i++ )
    {
        (*lookuptable_position)[i] = -1.0;
        (*inverse_lookuptable_position)[i] = 0.0;
    }

    for( i = 0; i < 16; i++ )
    {
        (*lookuptable_bp)[i] = 0;
    }

    return;
}


void
print_lookup_table( float* lookuptable_position, 
                    float* inverse_lookuptable_position,
                    int seq_len )
{
    int i;
    for( i = 0; i < seq_len; i++ )
    {        
        fprintf(stdout, "%i\t%e\t%e\n", 
                i+1, lookuptable_position[i], inverse_lookuptable_position[i]);
    }
}

inline float 
multiple_letter_penalty( const LETTER_TYPE* const reference,
                         const LETTER_TYPE* const observed,
                         const int start_position,
                         const int seq_length,
                         const int num_letters,
                         const float min_penalty,
                         const float* const lookuptable_position,
                         const float* const inverse_lookuptable_position,
                         const float* const lookuptable_bp
    )
{
    int j;
    float cum_penalty = 0;
    for( j = 0; j < num_letters - start_position; j++ )
    {
        /* determine the penalty contribution from letter j in seq i */ 
        float added_penalty = penalty_func( 
            reference[j], observed[j],
            start_position + j, 
            seq_length,
            min_penalty - cum_penalty,
            lookuptable_position,
            inverse_lookuptable_position,
            lookuptable_bp
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

extern void 
determine_bp_mut_rates( float** lookuptable_bp_ref )
{
    *lookuptable_bp_ref = malloc(16*sizeof(float));
    float* lookuptable_bp = *lookuptable_bp_ref;

    lookuptable_bp[4*0+0] = -1000;                /*A->A*/
    lookuptable_bp[4*0+1] = -0.4771213;           /*A->C*/
    lookuptable_bp[4*0+2] = -0.4771213;           /*A->G*/
    lookuptable_bp[4*0+3] = -0.4771213;           /*A->T*/

    lookuptable_bp[4*1+0] = -0.4771213;           /*C->A*/   
    lookuptable_bp[4*1+1] = -1000;                /*C->C*/
    lookuptable_bp[4*1+2] = -0.4771213;           /*C->G*/
    lookuptable_bp[4*1+3] = -0.4771213;           /*C->T*/

    lookuptable_bp[4*2+0] = -0.4771213;           /*G->A*/
    lookuptable_bp[4*2+1] = -0.4771213;           /*G->C*/
    lookuptable_bp[4*2+2] = -1000;                /*G->G*/
    lookuptable_bp[4*2+3] = -0.4771213;           /*G->T*/

    lookuptable_bp[4*3+0] = -0.4771213;           /*T->A*/
    lookuptable_bp[4*3+1] = -0.4771213;           /*T->C*/    
    lookuptable_bp[4*3+2] = -0.4771213;           /*T->G*/
    lookuptable_bp[4*3+3] = -1000;                /*T->T*/   

    return;         
}

inline int 
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

float 
recheck_penalty( char* reference, 
                 char* observed, 
                 const int seq_length,
                 const float* n_lookuptable_position,
                 const float* n_inverse_lookuptable_position,
                 const float* n_lookuptable_bp )
{
    int i = 0;
    float penalty = 0;
    for( i = 0; i < seq_length; i++ )
    {
        char ref = toupper( reference[i] );
        char obs = toupper( observed[i] );
        if( ref == obs )
        {
            penalty = penalty 
                + n_inverse_lookuptable_position[ i ] ;
        } else {
            penalty = penalty 
                + n_lookuptable_position[ i ];
            
            int bp_lookup = bp_code( reference[i] );
            bp_lookup = (bp_lookup << 2);
            bp_lookup += bp_code( observed[i] );
            
            penalty += n_lookuptable_bp[ bp_lookup ];
        }
    }
    
    return penalty;
}


float 
penalty_func( /* 
               * these aren't const because they are copies. I bit shift
               * them while calculating the penalty.
               */
              LETTER_TYPE reference, 
              LETTER_TYPE observed, 
              /* the position in the sequence - this should be
                 zero indexed */
              const int position, 
              /* the length of a full sequence */
              const int seq_length,
              /* the minimum allowable penalty */
              const float min_penalty,
              const float* const n_lookuptable_position,
              const float* const n_inverse_lookuptable_position,
              const float* const n_lookuptable_bp )
{ 
    int i;    
    float penalty = 0; 

    assert( position >= 0 );

    // it doesnt matter that the data type sapce is too big, givem
    // the alignment concerns
    for ( i = 0; i < LETTER_LEN; i++)
    {    
        /* 
            The naive way of doing this...        
        
        unsigned short ref_bp = (unsigned short) reference;
        ref_bp = ref_bp >> 2*i;      // ref_bp /= pow(4, i);
        ref_bp = ref_bp&3;  // ref_bp = ref_bp%4;

        unsigned short seq_bp = (unsigned  short) observed;
        seq_bp = seq_bp >> 2*i;      // seq_bp /= pow(4, i);
        seq_bp = seq_bp&3;  // seq_bp = seq_bp%4;
        */
        
        /* 
         * If the basepairs are the different, find out what the penalty is.
         * TODO make this an XOR for a ( very ) minor speed up
         */
        if ( (reference&3) != (observed&3) ) 
        {
            /* penalty due to the position */
            /* 
             * position is the position of the node in the tree, ie level. Thus, 
             * since we are considering the lowest order bits first ( we bit shift 
             * right at each iteration ) and we bitmask by 3, the bp_position
             * is given by level*LETTER_LEN - i. ie, if the letter length is 4
             * and we are in the first loop iteration( i = 0 ) then for TCGT we have
             * 10110110, so 10110110&00000011 = 00000010 which corresponds to the 
             * 4th letter in 1 indexing, or LETTER_LEN - 1 in 0 indexing.
             */

            /* TODO a potential optimization may be to cache the lookup table */
            penalty = penalty 
                    + n_lookuptable_position[ LETTER_LEN*position + i ]
                    + n_lookuptable_bp[ ((reference&3)<<2) + (observed&3) ];
            
        } 
        /* otherwise, calulate the penalty if they do match */
        else {
            if( LETTER_LEN*position + i < seq_length  ) {
                penalty += n_inverse_lookuptable_position[ LETTER_LEN*position + i ];
            }
        }

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

