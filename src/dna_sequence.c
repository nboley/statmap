/* Copyright (c) 2009-2012, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "dna_sequence.h"
#include "log.h"

/*
 *
 * Functions and data types for manipulating packed sequences.
 *
 */

/*
 * Change all N's to A's
 * 
 * If we are trying to map sequences that contain N's,
 * then we do a lowest quality sequence subsitution and
 * then change the actual basepair to 'A' ( arbitrarily ).
 * This is equivalent to us knowing nothing about the BP.
 *
 */
void
replace_ns_inplace( char* read, int seq_len )
{
    /* take the read complement */
    int i;
    for( i = 0; i < seq_len; i++ )
    {
        if( read[i] == 'N' || read[i] == 'n' )
        {
            read[i] = 'A';
        }
    }

    return;
}

int
calc_num_letters( const int seq_len )
{
    if( seq_len % LETTER_LEN == 0 )
        return seq_len/LETTER_LEN;

    return seq_len/LETTER_LEN + 1;
}


/* returns 0 on success, 5 on error */
int
complement_read( char* read, char* mut_read, int seq_len )
{
    /* take the read complement */
    int i;
    for( i = 0; i < seq_len; i++ )
    {
        char bp = read[i];
        switch(bp) {
        case 'A':
        case 'a':
            mut_read[i] = 'T';
            break;
        case 'C':
        case 'c':
            mut_read[i] = 'G';
            break;
        case 'G':
        case 'g':
            mut_read[i] = 'C';
            break;
        case 'T':
        case 't':
            mut_read[i] = 'A';
            break;
        case 'N':
        case 'n':
            return 5;
        default:
            statmap_log( LOG_ERROR, "Unable to translate %c", bp );
            return 5;
        }
    }
    
    mut_read[seq_len] = '\0';

    return 0;
};

int
rev_complement_read( char* read, char* mut_read, int seq_len )
{
    /* take the read complement */
    int i;
    for( i = 0; i < seq_len; i++ )
    {
        char bp = read[i];
        switch(bp) {
        case 'A':
        case 'a':
            mut_read[seq_len-1-i] = 'T';
            break;
        case 'C':
        case 'c':
            mut_read[seq_len-1-i] = 'G';
            break;
        case 'G':
        case 'g':
            mut_read[seq_len-1-i] = 'C';
            break;
        case 'T':
        case 't':
            mut_read[seq_len-1-i] = 'A';
            break;
        case 'N':
        case 'n':
            mut_read[seq_len-1-i] = 'N';
            break;
        default:
            statmap_log( LOG_ERROR, "Unable to translate %c", bp );
            return 5;
        }
    }
    
    mut_read[seq_len] = '\0';

    return 0;
};


/* returns 0 on success, 5 on error */
int
reverse_read( char* read, char* mut_read, int seq_len )
{
    /* take the read complement */
    int i;
    for( i = 0; i < seq_len; i++ )
    {
        mut_read[i] = read[seq_len-1-i];
    }
    mut_read[seq_len] = '\0';
    
    return 0;
};


void 
convert_packed_sequence( LETTER_TYPE* seq, char* new_seq, int seq_len )
{
    int i;
    for( i = 0; i < seq_len; i++)
    {
        /* TODO - MAKE Valgrind not complain */
        LETTER_TYPE letter = seq[(i/LETTER_LEN)];
        int shift_amt = 6 - 2*(i%LETTER_LEN);        
        letter = (letter >> shift_amt);
        letter = letter&3;
        char bp = 0;
        
        switch( letter ) 
        {
        case 0:
            bp = 'A';
            break;
        case 1:
            bp = 'C';
            break;
        case 2:
            bp = 'G';
            break;
        case 3:
            bp = 'T';
            break;
        default:
            bp = 'N';
            break;
        }
        new_seq[i] = bp;
    }
    new_seq[seq_len] = '\0';
    return;
};

void
print_packed_sequence( LETTER_TYPE* seq, int seq_length )
{
    char conv_seq[100];
    convert_packed_sequence( seq, conv_seq, seq_length );
    printf( "%s\n", conv_seq );
}

LETTER_TYPE* 
translate_seq(char* seq, int seq_len)
{
    unsigned short result_len, trans_val;
    unsigned char bp;

    /* determine the length of the result sequence */
    result_len = calc_num_letters( seq_len );

    /* calloc initializes the memory to 0 */
    LETTER_TYPE* result = calloc(result_len, sizeof(LETTER_TYPE));

    /* Memory check */
    if( result == NULL ) {
        statmap_log( LOG_FATAL, "Out of memory in translate_seq(). ( there is probably a memory leak )" );
    }

    /*
     * This loop is the core of the translate sequence, and a little 
     * confusing because it bitpacks the letters. 
     *
     * Basically, the loop goes through each character in the passed string,
     * sequence. Noting that the string is read in the forward direction, ie
     * if the sequence is "ACTAAT" then we take the left to be the 5' end. 
     * Also, we bit pack sequences in the "reverse" order in the sense that the
     * higher power bits store the *last* sequences. That is, for the example
     * above, if the letter length is 2, then we have 'AC' 'TA' 'AT' which 
     * correspond with ( for 'AC' ) 4**1 * 'C' + 4**0 * 'A' = 4 + 0 = 00000010
     * So, the full sequence is (00000010)(00000011)(00001100), noting that the
     * leftmost bits do nothing. If we bitpack this into a byte ( which makes 
     * much more sense, we have 'ACTA''ATAA' where the last two A's are ignored 
     * by the global sequence length. That is, we just ignore the trailing 
     * values. Thus, we get (00110100)(00001100), which is 
     * 4**3 * 'A' + 4**2 * 'C' + 4**1 * 'T' + 4**0 * 'A' ( you get the idea, if 
     * the power notation is confusing, think about the base ten case )
     *
     * We bitpack the letters in reverse order to facilitate a quicker penalty 
     * function. Since we generally care about the higher order terms first, 
     * ( because they are less likely to show mutation, which allows us to 
     *  abort the loop for a given translation value ) we want to be able to
     * bitwise and with 3, index the penalty array directly, and then bitshift
     * the high bits off of the end. If the lower order terms were stored first,
     * we would need to bitshift left and then re-bitshift right to  get the 
     * correct values. See the penalty function for details.
     */

    /* 
     * TODO - this function can be optimized if we consider the letters 
     * in order and then use bitwise operations. Of course, this would lower the 
     * potential generality, since then we could only use LETTER_LEN in powers
     * of 2, but im not convinced that this matters.  If we only let letters 
     * to be in powers of 2, then the best algorithm would probably be just to 
     * add the highest order term into the rightmost bits, and then bitshift 2 
     * every single iteration.
     */

    int index;
    for ( index=0; index < result_len; index++ )
    {
        /*
         * There is some subtlty in the following when sequence lengths
         * are not evenly divisible by letter lengths. If not, ( consider 
         * 'TCG' with LETTER_LEN = 2. Then we want (0111)(0010). Thus, we 
         * really just want 'TCGA', which is how we implement it
         */
        LETTER_TYPE current_letter = 0;
        int bit_index;
        for( bit_index = 0; bit_index < LETTER_LEN; bit_index++ )
        {
            /* 
             * note that, for bit index = 0, we take the *lowest* order bp 
             * since we are bitshifting, this will finish as the highest. Also,
             * see above for how we are dealing with phantom bits.
             * TODO optimize this with a macro TODO
             */
            unsigned short bp_index = LETTER_LEN*index + bit_index;
            if( bp_index >= seq_len ) {
                bp = 'A'; 
            } else {
                bp = seq[bp_index];
            }

            switch( bp )
            {
                case 'A':
                case 'a':
                    trans_val = 0;
                    break;
                case 'C':
                case 'c':
                    trans_val = 1;
                    break;
                case 'G':
                case 'g':
                    trans_val = 2;
                    break;
                case 'T':
                case 't':
                    trans_val = 3;
                    break;
                case 'N':
                case 'n':
                    free(result);
                    return NULL;
                default:          
                    free(result);
                    fprintf (stderr, "Error in the sequence - read '%c'\n", bp);
                    return NULL;
            }
            current_letter = current_letter << 2;
            current_letter = (current_letter | trans_val);
        }

        /* set the current letter in the array */
        result[index] = current_letter;
    }
    
    return result;
}

 LETTER_TYPE*
copy_first_k_basepairs( LETTER_TYPE* read, 
                        LETTER_TYPE** new_read, 
                        int read_len, int num_letters,
                        int k)
{
    if( k > read_len ) {
        *new_read = NULL;
        return NULL;
    }
    
    const int new_num_letters = calc_num_letters( k );

    *new_read = calloc(new_num_letters, sizeof(LETTER_TYPE));
    
    /* copy the first new_num_letters letter */
    memcpy( *new_read, read, new_num_letters );

    /* if k%num_letters == 0, we are done */
    if( k%num_letters == 0 )
        return *new_read;

    /* otherwise, zero out the final basepairs */
    /* generate the final letters mask - 
       noting that higher order letters are first */
    /* int num_bps_in_last_letter = k%num_letters */
    LETTER_TYPE mask = ( 1 << ( 2*(k%LETTER_LEN) ) ) - 1;
    (*new_read)[ new_num_letters - 1 ] &= mask;

    return *new_read;
}

 LETTER_TYPE*
copy_last_k_basepairs(  LETTER_TYPE* read, 
                        LETTER_TYPE** new_read, 
                        int read_len, int k)
{
    if( k > read_len ) {
        *new_read = NULL;
        return NULL;
    }
    
    // const int num_letters = num_letters_needed( read_len );
    const int new_num_letters = calc_num_letters( k );
    
    *new_read = calloc(new_num_letters, sizeof(LETTER_TYPE));
    
    int i;
    for( i=0; i < k; i++ )
    {
        /* find the index of the read's bp */
        int read_bp_index = (read_len - k)+i;
        int read_letter_index = read_bp_index/LETTER_LEN;
        int read_bp_offset = read_bp_index%LETTER_LEN;
        int read_bp = (((read[read_letter_index])&(3<<(2*read_bp_offset)))>>(2*read_bp_offset));


        /* find the letter index of the i'th bp */
        int new_letter_index = i/LETTER_LEN;
        int new_bp_offset = i%LETTER_LEN;

        /* add the new letter */
        ((*new_read)[new_letter_index]) |= (read_bp  << (2*new_bp_offset));
    }

    return *new_read;
}


 int 
cmp_letters( LETTER_TYPE letter1, LETTER_TYPE letter2 )
{
    /* 
     * compare two letters in *reverse* bit order. ie, pay attention to
     * the lower order bits before the higher.
     */
    int i;
    /* check the lowest order bits. If they are identical, shift down... */
    for( i = 0; i < LETTER_LEN; i++ )
    {
        if( (letter1&3) > (letter2&3)) {
            return -1;
        } 
        if( (letter2&3) > (letter1&3) ) {
            return 1;
        } 

        letter1 >>= 2;
        letter2 >>= 2;
    }

    return 0;
}

 int 
cmp_words(   LETTER_TYPE* seq1, 
             LETTER_TYPE* seq2, 
             const int num_letters )
{
    int i;
    for( i = 0; i < num_letters; i++ ) 
    {
        int value = cmp_letters( seq1[i], seq2[i] );
        if( value != 0 )    {
            return value;
        }
    }
    return 0;

}
