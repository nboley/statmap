/* Copyright (c) 2009-2010, Nathan Boley */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <fcntl.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <sys/stat.h> /* permission bit defines */

#include "config.h"
#include "index_genome.h"
#include "genome.h"
#include "snp.h"
#include "pseudo_location.h"

/**** ON DISK Code **********************************************************************/

/* returns size written to disk */
static size_t
write_reference_data_header_to_disk( struct genome_header* header, FILE* fp )
{
    int rv;
    fprintf( fp, "SM_OD_GEN" );
    
    rv = fwrite( &(header->size), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fwrite( &(header->genome_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fwrite( &(header->pseudo_locs_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fwrite( &(header->snp_db_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );
    
    /* // keep the index in a separate file
    rv = fwrite( &(header->index_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );
    */
    
    return ( 9*sizeof(char) + 4*sizeof(size_t) );
}

static void
read_reference_data_header_from_disk( struct genome_header* header, FILE* fp )
{
    int rv;
    
    char magic_number[9];
    rv = fread( magic_number, sizeof(unsigned char), 9, fp );

    #ifdef DEBUG
    fprintf( stderr, "DEBUG       :  Genome Magic Number - %.9s\n", magic_number );
    #endif

    assert( rv == 9 );
    if( 0 != memcmp( magic_number, "SM_OD_GEN", 9 ) )
    {
        fprintf( stderr, "FATAL       :  Genome Magic Number ('%.9s') is incorrect ( it should be 'SM_OD_GEN' - cmp %i  ) \n", 
                 magic_number, memcmp( magic_number, "SM_OD_GEN", 9 ) );
        fprintf( stderr, "HINT        :  Is this a fasta file? Fasta files need to be converted with build_index.\n" );
        exit( 1 );
    }

    rv = fread( &(header->size), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fread( &(header->genome_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fread( &(header->pseudo_locs_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fread( &(header->snp_db_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    /* // keep the index in a separate file
    rv = fread( &(header->index_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );
    */
    
    return;
}

/* END on disk code **********************************************************/

void
init_genome( struct genome_data** gen )
{
    *gen = malloc( sizeof( struct genome_data ) );

    /* init the raw genome info */
    (*gen)->num_chrs = 0;
    (*gen)->chr_names = NULL;
    (*gen)->chrs = NULL;
    (*gen)->chr_lens = NULL;

    /* the genome index */
    (*gen)->index = NULL;
    
    /* snps - we dont init because they may not exist */
    (*gen)->snp_db = NULL;

    (*gen)->is_mmapped = false;

    /** Init the pseudo locations **/
    (*gen)->ps_locs = NULL;    
    init_pseudo_locations( &((*gen)->ps_locs)  );
    
    /* Add the pseudo chromosome */
    add_chr_to_genome( "Pseudo", "", 0, *gen );
}

/* returns the size written */
size_t
write_standard_genome_to_file( struct genome_data* gen, FILE* ofp  )
{
    size_t size_written = 0;
    int rv;
    int i;
    
    /* write the number of chromosomes */
    rv = fwrite( &(gen->num_chrs), sizeof( int ), 1, ofp );
    size_written += sizeof( int );
    assert( rv == 1 );

    /* write the chr names */
    for( i = 0; i < gen->num_chrs; i++ )
    {
        /* DONT include the terminating null */
        int chr_name_len = strlen( gen->chr_names[i] ) + 1;
        
        /* write the number of characters in the chr name */
        rv = fwrite( &chr_name_len, sizeof( int ), 1, ofp );
        size_written += sizeof( int );
        assert( rv == 1 );
        
        /* write the chr name */
        rv = fwrite( gen->chr_names[i], sizeof( char ), chr_name_len, ofp );
        assert( rv == chr_name_len );
        size_written += chr_name_len*sizeof( char );
    }

    /* write the chr lengths */
    rv = fwrite( gen->chr_lens, sizeof( unsigned int ), gen->num_chrs, ofp );
    size_written += (gen->num_chrs)*sizeof( unsigned int );
    assert( rv == gen->num_chrs);
    
    
    /* write the actual chromosomes */
    for( i = 0; i < gen->num_chrs; i++ )
    {
        /* write the number of characters in the chr name */
        rv = fwrite( gen->chrs[i], sizeof( char ), gen->chr_lens[i], ofp );
        size_written += sizeof( char )*gen->chr_lens[i];
        assert( rv == (int) gen->chr_lens[i] );        
    }
    
    return size_written;
}


void
populate_standard_genome_from_mmapped_file( struct genome_data* gen, char* data  )
{    
    int i;

    /* get the number of chromosomes */
    gen->num_chrs = *((int*) data);
    data += sizeof( int );
    
    /* read the chr names */
    gen->chr_names = malloc( gen->num_chrs*sizeof(char*) );
    for( i = 0; i < gen->num_chrs; i++ )
    {
        /* include the terminating null */
        int chr_name_len = *((int*) data);
        data += sizeof( int );
        
        gen->chr_names[i] = data;
        data += chr_name_len*sizeof( char );
    }

    /* read the chr lengths */
    gen->chr_lens = (unsigned int*) data;
    data += sizeof(unsigned int)*gen->num_chrs;
    
    /* read the actual chromosomes */
    gen->chrs = malloc( gen->num_chrs*sizeof(char**) );
    for( i = 0; i < gen->num_chrs; i++ )
    {
        gen->chrs[i] = data;
        data += sizeof(char)*gen->chr_lens[i];
    }
    
    return;
}


void
write_genome_to_disk( struct genome_data* gen, char* fname )
{
    FILE* ofp;
    ofp = fopen( fname, "w+" );
    if( ofp == NULL )
    {
        fprintf( stderr, "Can not open '%s' for writing.", fname );
        exit( 1 );
    }
      
    size_t size_written;
    
    /** write the header  */
    /* 
       we dont actually know what any of these values are, so we write 
       them to allocate the sapce, then we will go back and make them 
       correct 
    */
    struct genome_header header;
    
    header.size = 0;
    header.genome_offset = 0;
    header.snp_db_offset = 0;
    header.pseudo_locs_offset = 0;
    // header.index_offset = 0;

    size_written = write_reference_data_header_to_disk( &header, ofp );
    header.size += size_written;
    header.genome_offset = size_written;
    
    size_written = write_standard_genome_to_file( gen, ofp  );
    header.size += size_written;
    header.snp_db_offset = header.genome_offset + size_written;
    
    size_written = write_snp_db_to_binary_file( gen->snp_db, ofp );
    header.size += size_written;
    header.pseudo_locs_offset = header.snp_db_offset + size_written;

    /* write the updated header */
    fseek( ofp, 0, SEEK_SET );
    size_written = write_reference_data_header_to_disk( &header, ofp );
    
    fclose( ofp );
    
    return;
}

void
load_genome_from_disk( struct genome_data** gen, char* fname )
{
    /* 
       first, 
       open the file containing the index to ensure that
       the magic number is correct and to get the size.
    */
    
    FILE* genome_fp = fopen( fname, "r" );
    if( genome_fp == NULL )
    {
        fprintf( stderr, "Cannot open the binary genome '%s' for reading", fname );
        exit( -1 );
    }

    struct genome_header header;
    read_reference_data_header_from_disk( &header, genome_fp );

    #ifdef DEBUG
    fprintf( stderr, "DEBUG       :  HEADER DATA: %zu %zu %zu %zu\n", 
            header.size, header.genome_offset, 
            header.snp_db_offset, header.pseudo_locs_offset  );
    #endif
    
    *gen = malloc( sizeof( struct genome_data ) );
    
    fclose( genome_fp );

    (*gen)->is_mmapped = true;
    
    int fd;
    if ((fd = open(fname, O_RDONLY )) < 0)
        fprintf(stderr, "FATAL     : can't create %s for reading", fname);
    
    char* OD_genome;
    assert( sizeof(char) == 1 );
    if ((OD_genome = mmap (0, header.size,
                          PROT_READ,
                          MAP_SHARED, fd, 0)) == (caddr_t) -1)
        fprintf(stderr, "FATAL     : mmap error for genome file");

    char* curr_pos = OD_genome + header.genome_offset;
    populate_standard_genome_from_mmapped_file( *gen, curr_pos );

    curr_pos = OD_genome + header.snp_db_offset;
    load_snp_db_from_mmapped_data( *gen, curr_pos );

    curr_pos = OD_genome + header.pseudo_locs_offset;
    load_pseudo_locations_from_mmapped_data( &((*gen)->ps_locs), curr_pos );
        
    return;
}


void
free_genome( struct genome_data* gen )
{
    int i;
    
    if( NULL == gen )
        return;
    
    if( gen->is_mmapped == true )
    {
        free( gen->chr_names );
        free( gen->chrs );
        free( gen );
        return;
    }

    /* BUG - make this recognize the tree */
    if( false && gen->index != NULL )
    {
        free_tree( gen->index );
    }

    if( gen->index != NULL )
    {
        if( NULL != gen->index->ps_locs )
        {
            free_pseudo_locations( gen->index->ps_locs );
        }

        free( gen->index );
    }
    
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

char* 
find_seq_ptr( struct genome_data* genome, 
              int chr_index, unsigned int loc, 
              int read_len )
{
    assert( chr_index < genome->num_chrs );
    
    if( chr_index != PSEUDO_LOC_CHR_INDEX 
        && ( 
            loc + read_len > genome->chr_lens[chr_index]
            )
        ) 
    {
        return NULL;
    }

    assert( chr_index == PSEUDO_LOC_CHR_INDEX 
            || ( 
                loc < genome->chr_lens[chr_index]
            )
        );
    
    /* 
       the pseudo locations need to be treated differently,
       because we don't store the chromosome length ( because
       that doesn't really mean anything. ) However, we do need
       to make sure that the pseudo location *exists* 
    */
    
    assert( chr_index != (int) PSEUDO_LOC_CHR_INDEX 
            || loc < (unsigned long) genome->index->ps_locs->num );    

    /* if this is a pseudo chromosome, we need to 
       get the sequence associated with it. Since
       the sequence for every location associated with 
       a pseudo location is the same, it is sufficient to 
       just get the sequence of the first location.
    */
    
    if( chr_index == PSEUDO_LOC_CHR_INDEX ) {
        GENOME_LOC_TYPE* ps_loc = &(genome->index->ps_locs->locs[loc].locs[0]);
        chr_index = ps_loc->chr;
        loc = ps_loc->loc;
        // fprintf( stderr, "%.20s\n", genome->chrs[chr_index] + loc );
    }

    return genome->chrs[chr_index] + loc;
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

        if( !isalpha(bp))
        {
            fprintf( stderr, "ERROR           : Unrecognized character '%c' - ignoring\n", bp );
            
        } else {
            chr[i] = bp;
        }
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

/*
 * Forward declaration for translate_sequence. I wanted to avoid including
 * DNA sequence.h for this single declaration 
 */
LETTER_TYPE* 
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
     * just always init the db but set the number of snps to 0.
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
     *
     * We skip chr 0 - the pseudo chromosome.
     */
    for( chr_index = 1; chr_index < genome->num_chrs; chr_index++ )
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
                
                /* Add the sequence into the tree */
                add_sequence(genome->index, genome->index->ps_locs, translation, seq_len, loc);
                
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
                    
                    /* Add the sequence into the tree */
                    add_sequence(genome->index, genome->index->ps_locs, 
                                 translation, seq_len, loc);
                    
                    free( translation );                
                }
            }
        }
    }
    free( tmp_seq );

    /* sort all of the pseudo locations */
    sort_pseudo_locations( genome->index->ps_locs );

    /* We stop doing this, in favor of writing them to 
       disk in the actual index file */
    #if 0 

    FILE* ps_locs_of = fopen(PSEUDO_LOCATIONS_FNAME, "w");
    if( NULL == ps_locs_of )
    {
        char buffer[500];
        sprintf( buffer, "Error opening '%s' for writing.", PSEUDO_LOCATIONS_FNAME );
        perror( buffer );
        exit( -1 );
    }
    fprint_pseudo_locations( ps_locs_of, genome->index->ps_locs );
    fclose(ps_locs_of);
    
    #endif
    
    return;
}



