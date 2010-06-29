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
    
    /* Add the pseduo chromosome */
    add_chr_to_genome( "Pseudo", "", 0, *gen );
}

void
free_genome( struct genome_data* gen )
{
    int i;
    
    if( gen->index != NULL )
        free_tree( gen->index );
    
    for( i = 0; i < gen->num_chrs; i++ )
    {
        if( (gen->chr_lens)[i] > 0 )
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
    assert( chr_len <= LOCATION_MAX );

    /* copy the chr string into memory */
    gen->chrs = 
        realloc( gen->chrs, sizeof(char*)*(gen->num_chrs) );
    (gen->chrs)[gen->num_chrs - 1] = malloc( sizeof(char)*(chr_len+1) );
    
    if( chr_len > 0 ) 
        memcpy( (gen->chrs)[gen->num_chrs - 1], chr_str, sizeof(char)*chr_len );

    (gen->chrs)[gen->num_chrs - 1][chr_len] = '\0';
}

extern void 
add_chrs_from_fasta_file( 
    struct genome_data* gen, FILE* f   )
{
    int error;
    
    /* FIXME - get rid of this via mmap */
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
    /* we set it to 0 to skip the pseudo chr */
    int chr_index = 0;

    unsigned int i = 0;
    while( !feof(f) )
    {
        char bp = (char) fgetc(f);
        if( !isalpha(bp) || (chr_index == 0 && i==0) ) 
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
                /* if we are past the first non pseudo chr, add the previous chr */
                if ( chr_index > 0 )
                {
                    /* 
                     * put a trailing null in to make it a proper string, 
                     * and add it to the genome
                     */
                    chr[i] = '\0';    
                    /* shrink the size of chr so that we arent wasting memory */
                    chr = realloc( chr, (i+1)*sizeof(char) );
                    allcd_size = i;

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
            } else if ( i == 0 && chr_index == 0 )
            {
                /* move to the next chr index */
                chr_index++;
                sprintf(chr_name, "Default");
            }
        }

        if( i == allcd_size ) {
            allcd_size += allcd_size;
            chr = realloc( chr, (allcd_size+1)*sizeof(char) );
            if( chr == NULL ) {
                fprintf(stderr, "FATAL       :  Could not allocate enough memory for genome\n");
                assert( 0 );
                exit(-1);
            }
        }        

        assert(isalpha(bp));
        chr[i] = bp;
        i++;        
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



