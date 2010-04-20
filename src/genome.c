/* Copyright (c) 2009-2010, Nathan Boley */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "config.h"
#include "index_genome.h"
#include "genome.h"
#include "snp.h"

void
init_genome( struct genome_data** gen )
{
    *gen = malloc( sizeof( struct genome_data ) );
    (*gen)->index = NULL;
    (*gen)->num_chrs = 0;
    (*gen)->chr_names = NULL;
    (*gen)->chrs = NULL;
    (*gen)->chr_lens = NULL;
    (*gen)->snp_db = NULL;
}

void
free_genome( struct genome_data* gen )
{
    int i;
    
    if( gen->index != NULL )
        free_tree( gen->index );

    for( i = 0; i < gen->num_chrs; i++ )
    {
        free( (gen->chrs)[i] );
        free( (gen->chr_names)[i] );
    }
    free( gen->chr_names );

    
    free( gen->chrs );
    free( gen->chr_lens );

    if( gen->snp_db != NULL )
        free_snp_db( gen->snp_db );

    free( gen );
}

int
find_chr_index( struct genome_data* genome, const char* const chr_name )
{
    int i = 0;
    for( i = 0; i < genome->num_chrs; i++ )
    {
        if( 0 == strcmp( chr_name, genome->chr_names[i] ) )
        {
            return i;
        }
    }
    
    return -1;
}

void
add_chr_to_genome( char* chr_name, 
                   char* chr_str, 
                   unsigned int chr_len,
                   struct genome_data* gen )
{
    /* find the length of the chromosome name  */
    int chr_name_len = strlen( chr_name );
    
    /* increment the number of chromosomes */
    gen->num_chrs += 1;

    /* copy the chr name into the array */
    gen->chr_names = 
        realloc( gen->chr_names, sizeof(char*)*(gen->num_chrs) );
    (gen->chr_names)[gen->num_chrs - 1] = malloc( sizeof(char)*(chr_name_len+1) );
    memcpy( (gen->chr_names)[gen->num_chrs - 1], chr_name, sizeof(char)*chr_name_len );
    (gen->chr_names)[gen->num_chrs - 1][chr_name_len] = '\0';

    /* put in the chr str len */
    gen->chr_lens = realloc( gen->chr_lens, sizeof(size_t)*(gen->num_chrs) );
    (gen->chr_lens)[gen->num_chrs - 1] = chr_len;
    
    /* copy the chr string into memory */
    gen->chrs = 
        realloc( gen->chrs, sizeof(char*)*(gen->num_chrs) );
    (gen->chrs)[gen->num_chrs - 1] = malloc( sizeof(char)*(chr_len+1) );
    memcpy( (gen->chrs)[gen->num_chrs - 1], chr_str, sizeof(char)*chr_len );
    (gen->chrs)[gen->num_chrs - 1][chr_len] = '\0';
}

extern void 
add_chrs_from_fasta_file( 
    struct genome_data* gen, FILE* f   )
{
    int error;
    
    /* FIXME - ge rid of this via mmap */
    // int max_num_bytes = MAX_GENOME_LOC;
    unsigned int allcd_size = 300000;

    /* BUG OVERFLOW */
    /* store the chromosome name - limit it to 255 characters */
    char chr_name[255];
    chr_name[0] = '\0';

    char* chr;
    chr = malloc( (allcd_size+1)*sizeof(char) );
    if( chr == NULL ) {
        fprintf(stderr, "FATAL       :  Could not allocate enough memory for genome\n");
        assert( 0 );
        exit(-1);
    }
    
    /* store the index of the current chr */
    int chr_index = -1;

    unsigned int i = 0;
    while( !feof(f) )
    {
        char bp = (char) fgetc(f);
        if( !isalpha(bp) || (chr_index == -1 && i==0) ) 
        {
            /* 
             *  if we read a new line, then the next line could be a label. 
             * Therefore, we read the next character and if it is not a >
             * then we continue as usual. If it is, we read until the next 
             * chr.
             */
            if( isspace(bp) ) {
                bp = (char) fgetc(f);
            }
            /* check to see if we are at the end of the file */
            if( feof(f) )
                break;
            
            if( bp == '>' )
            {
                /* if we are past the first chr, add the previous chr */
                if ( chr_index >= 0 )
                {
                    /* 
                     * put a trailing null in to make it a proper string, 
                     * and add it to the genome
                     */
                    chr[i] = '\0';    
                    /* shrink the size of chr so that we arent wasting memory */
                    chr = realloc( chr, (i+1)*sizeof(char) );

                    assert( chr_name[0] != '\0' );
                    fprintf(stderr, "NOTICE      :  Added '%s' with %i basepairs\n", 
                            chr_name, i);
                    add_chr_to_genome( chr_name, chr, i, gen );
                    
                    if( !feof(f)) {
                        /* reset the bp index to 0 */
                        i = 0;
                    }
                }

                /* move to the next chr index */
                chr_index++;

                /* OVERFLOW BUG - check chr_name */
                /* read in the next string as the chromosome */
                error = fscanf(f, "%s\n", chr_name );

                /* if we are at the end of the file, break */
                if( feof(f) ) {
                    break;
                } else {
                    bp = (char) fgetc(f);
                }                    
            } else if ( i == 0 && chr_index == -1 )
            {
                /* move to the next chr index */
                chr_index++;
                sprintf(chr_name, "Default");
            }
        }
        
        assert(isalpha(bp));
        chr[i] = bp;
        i++;

        if( i == allcd_size ) {
            allcd_size += allcd_size;
            chr = realloc( chr, (allcd_size+1)*sizeof(char) );
            if( chr == NULL ) {
                fprintf(stderr, "FATAL       :  Could not allocate enough memory for genome\n");
                assert( 0 );
                exit(-1);
            }
        }
    }
    
    /* put a trailing null in to make it a proper string, and add the chr */
    chr[i] = '\0';    

    /* shrink the size of chr so that we arent wasting memory */
    chr = realloc( chr, (i+1)*sizeof(char) );

    assert( chr_name[0] != '\0' );
    add_chr_to_genome( chr_name, chr, i, gen );
    
    fprintf(stderr, "NOTICE      :  Added '%s' with %i basepairs\n", chr_name, i);
    free(chr);


    /* close the input file */
    fclose(f);

    return;
}


/*
 * Forward declaration for translate_sequence. I wanted to avoid including
 * DNA sequence.h for this single declaration 
 */
inline LETTER_TYPE* 
translate_seq(char*, int num_letters, LETTER_TYPE** );

extern void
index_genome( struct genome_data* genome, int indexed_seq_len )
{
    /* initialize the tree structure */
    init_tree( &(genome->index), indexed_seq_len );

    /* TODO - use snps directly */
    struct snp_db_t* snp_db = genome->snp_db;
    /* 
     * If there is no initialized snp DB, then we set the
     * number of snps to zero so that we dont try and 
     * add any downstream. Maybe, in the future, I will 
     * just always the db but set the number of snps to 0.
     *
     */
    int num_snps;
    if( snp_db == NULL )
    {
        num_snps = 0;
    } else {
        num_snps = snp_db->num_snps;
    }
    
    /* initialize the constant loc elements */
    GENOME_LOC_TYPE loc;
    /* Not a junction read */
    loc.read_type = 0;
    
    int seq_len = genome->index->seq_length;
    char* tmp_seq = malloc(seq_len*sizeof(char));

    int chr_index;
    unsigned int bp_index;
    
    int snp_lb = 0;
    int snp_ub = 0;

    /* 
     * Iterate through each indexable sequence in the genome. If the 
     * sequence covers snps, then add them.
     */
    for( chr_index = 0; chr_index < genome->num_chrs; chr_index++ )
    {
        if( genome->chr_lens[chr_index] > LOCATION_MAX )
        {
            fprintf( stderr, 
                     "FATAL ERROR       :  Max indexable chr size is '%i'\n", 
                     LOCATION_MAX );
            exit( -1 );
        }

        fprintf(stderr, "NOTICE      :  Indexing '%s'\n", (genome->chr_names)[chr_index] );

        /* Set the chr index in the soon to be added location */
        loc.chr = chr_index;

        for( bp_index = 0; bp_index < genome->chr_lens[chr_index] - seq_len; bp_index += 1 )
        {            
            /* Set the basepair index */
            loc.loc = bp_index;
            
            /* 
             * Update the snp_lb. While the snps are strictly below bp_index, the
             * lower bound, then the snps can not fall into the sequence and there
             * is no need to consider these
             */
            while( snp_lb < num_snps
                   && snp_db->snps[snp_lb].loc.chr == chr_index 
                   && snp_db->snps[snp_lb].loc.loc < bp_index  )
            {
                snp_lb += 1;
            }
                        
            /* 
             * Update the snp_ub. While the snps are strictly below bp_index, the
             * lower bound, then the snps can not fall into the sequence and there
             * is no need to consider these
             */
            snp_ub = MAX( snp_lb, snp_ub );
            while( snp_ub < num_snps
                   && snp_db->snps[snp_ub].loc.chr == chr_index 
                   && snp_db->snps[snp_ub].loc.loc < bp_index + seq_len )
            {
                snp_ub += 1;
            }
                   
            /* 
             * If lb == ub, then we know this sequence doesnt cover any snps. 
             */
            if( snp_lb == snp_ub )
            {
                memcpy( tmp_seq, genome->chrs[chr_index] + bp_index, sizeof(char)*seq_len );
                
                /* Add the normal sequence */
                LETTER_TYPE *translation;
                translate_seq( tmp_seq, seq_len, &(translation));
                
                /* if we cant add this sequence ( probably an N ) continue */
                if( translation == NULL ) {
                    continue;
                }
                
                loc.covers_snp = 0;
                loc.snp_coverage = 0;
                
                /* Add the sequence into the tree */
                add_sequence(genome->index, translation, seq_len, loc);
                
                free( translation );                                
            } else {
                /* If there are too many snps in this sequence, print a warning */
                if( snp_ub - snp_lb > MAX_NUM_SNPS )
                {
                    fprintf(stderr, "ERROR       :  Can not add the sequence at 'chr %i: bp %i' - it contains %i SNPs ( %i max )\n", 
                            chr_index, bp_index, snp_ub - snp_lb, MAX_NUM_SNPS );
                    continue;
                }
                
                /* 
                 * We need to iterate through every combination of snps in the following
                 * list. Ie, if there are 2 snps, then we need to add all of the 
                 * subsequences with s1 on, s2 on, sn1 on, s2 off, s1 off, s2 off, 
                 * and both off.  This might be hard in general, but there
                 * is an easy way of dealing with this. Since there are 2^NS -1 total
                 * combinations, we just count from 1 to 2^NS. Then, if the ith bit
                 * is set, we know that the ith snp should be on. 
                 */
                
                /* Make sure there is room in the bitmap to store all of the snips */
                assert( sizeof(unsigned int)*8 > MAX_NUM_SNPS  );
                /* bm stands for bitmap */
                unsigned int bm;
                assert( snp_ub > snp_lb );
                for( bm = 0; bm < (((unsigned int)1) << (snp_ub - snp_lb)); bm++ )
                {
                    /* Make a copy of the sequence that we will be mutating */
                    /* 
                     * TODO - Make this more efficient. We shouldnt need to recopy the sequence
                     * every single time. Although, since snps should typically be pretty 
                     * sparse, this may not matter in practice.
                     */
                    memcpy( tmp_seq, genome->chrs[chr_index] + bp_index, sizeof(char)*seq_len );
                    
                    /* 
                     * Loop through each possible snp. If the bit is set, then 
                     * set the bp location in the seq to the alternate.
                     */
                    int snp_index; 
                    for( snp_index = 0; snp_index < snp_ub - snp_lb; snp_index++ )
                    {
                        /* If the correct bit is set */
                        if( (bm&(1<<snp_index)) > 0 )
                            tmp_seq[ snp_db->snps[snp_lb+snp_index].loc.loc - bp_index  ] 
                                = snp_db->snps[snp_lb+snp_index].alt_bp;                        
                    }
                    
                    LETTER_TYPE *translation;
                    translate_seq( tmp_seq, seq_len, &(translation));
                    
                    /* if we cant add this sequence ( probably an N ) continue */
                    if( translation == NULL ) {
                        continue;
                    }
                    
                    loc.covers_snp = 1;
                    loc.snp_coverage = bm;
                    
                    /* Add the sequence into the tree */
                    add_sequence(genome->index, translation, seq_len, loc);
                    
                    free( translation );                
                }
            }
        }
    }
    free( tmp_seq );

    return;
}

