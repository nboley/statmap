/* Copyright (c) 2009-2012, Nathan Boley */

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
#include "pseudo_location.h"
#include "diploid_map_data.h"

#include "log.h"

long
calc_genome_len( struct genome_data* genome )
{
    long genome_len = 0;
    
    int i;
    for( i = 0; i < genome->num_chrs; i++ )
    {
        genome_len += genome->chr_lens[i];
    }
    
    return genome_len;
}


/**** ON DISK Code **********************************************************************/

/* returns size written to disk */
static size_t
write_reference_data_header_to_disk( struct genome_header* header, FILE* fp )
{
    int rv;
    fprintf( fp, "SM_OD_GEN" );

    rv = fwrite( &(header->format_version), sizeof(size_t), 1, fp );
    assert( rv == 1 );
    
    rv = fwrite( &(header->size), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fwrite( &(header->genome_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    return ( 9*sizeof(char) + 3*sizeof(size_t) );
}

static void
read_reference_data_header_from_disk( struct genome_header* header, FILE* fp )
{
    int rv;
    
    char magic_number[9];
    rv = fread( magic_number, sizeof(unsigned char), 9, fp );

    assert( rv == 9 );
    if( 0 != memcmp( magic_number, "SM_OD_GEN", 9 ) )
    {
        statmap_log( LOG_FATAL,
                "Genome Magic Number ('%.9s') is incorrect ( it should be 'SM_OD_GEN' - cmp %i  )",
                magic_number,
                memcmp( magic_number, "SM_OD_GEN",  9 )
            );
    }

    rv = fread( &(header->format_version), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    if( header->format_version != GENOME_FILE_FORMAT_VERSION )
    {
        statmap_log( LOG_FATAL,
                "Genome file format version is incompatible with this release of Statmap (found %zu, require %zu)",
                header->format_version,
                GENOME_FILE_FORMAT_VERSION
            );
    }

    rv = fread( &(header->size), sizeof(size_t), 1, fp );
    assert( rv == 1 );

    rv = fread( &(header->genome_offset), sizeof(size_t), 1, fp );
    assert( rv == 1 );
    
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
    (*gen)->chr_sources = NULL;

    /* the genome index */
    (*gen)->index = NULL;
    
    (*gen)->is_mmapped = false;

    /* Add the pseudo chromosome */
    add_chr_to_genome( "Pseudo", "", 0, REFERENCE, *gen );
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

    /* write the chr sources */
    rv = fwrite( gen->chr_sources, sizeof( enum CHR_SOURCE ), gen->num_chrs, ofp );
    size_written += (gen->num_chrs)*sizeof( enum CHR_SOURCE );
    assert( rv == gen->num_chrs );
    
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

    /* read the chr sources */
    gen->chr_sources = (enum CHR_SOURCE *) data;
    data += sizeof(enum CHR_SOURCE)*gen->num_chrs;
    
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
    if( ofp == NULL ) {
        statmap_log( LOG_FATAL, "Could not open '%s' for writing", fname );
    }
      
    size_t size_written;
    
    /** write the header  */
    /* 
       we don't actually know what the values of size or genome_offset are yet,
       so we write them to allocate the space, then we will go back and make
       them correct 
    */
    struct genome_header header;
    
    header.format_version = GENOME_FILE_FORMAT_VERSION;
    header.size = 0;
    header.genome_offset = 0;

    size_written = write_reference_data_header_to_disk( &header, ofp );
    header.size += size_written;
    header.genome_offset = size_written;
    
    size_written = write_standard_genome_to_file( gen, ofp  );
    header.size += size_written;
    
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
    if( genome_fp == NULL ) {
        statmap_log( LOG_FATAL, "Cannot open the binary genome '%s' for reading", fname );
    }

    struct genome_header header;
    read_reference_data_header_from_disk( &header, genome_fp );
    
    *gen = malloc( sizeof( struct genome_data ) );
    
    fclose( genome_fp );

    (*gen)->is_mmapped = true;
    
    int fd;
    if ((fd = open(fname, O_RDONLY )) < 0)
        statmap_log( LOG_FATAL, "can't create %s for reading",  fname );
    
    char* OD_genome;
    assert( sizeof(char) == 1 );
    if ((OD_genome = mmap (0, header.size,
                          PROT_READ,
                          MAP_SHARED, fd, 0)) == (caddr_t) -1)
        statmap_log( LOG_FATAL, "mmap error for genome file" );

    char* curr_pos = OD_genome + header.genome_offset;
    populate_standard_genome_from_mmapped_file( *gen, curr_pos );

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

    if( gen->index != NULL )
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
    free( gen->chr_sources );

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

/*
 * Returns the index of the chr matching the given prefix and source
 * Returns -1 if no chr matches
 */
int
find_diploid_chr_index(
    struct genome_data* genome,
    const char* const prefix,
    enum CHR_SOURCE source
)
{
    int i;
    /* loop over chrs */
    for( i=0; i < genome->num_chrs; i++ )
    {
        if( genome->chr_sources[i] == source &&
            0 == strncmp( genome->chr_names[i], prefix, strlen(prefix) )
          )
            return i;
    }
    return -1;
}

char* 
find_seq_ptr( struct genome_data* genome, 
              int chr_index, int start_bp, 
              int read_len )
{
    assert( chr_index < genome->num_chrs );

    /* We can get negative locs from the assay specific correction code if the
     * original mapping is close to a contig boundary. */
    if( start_bp < 0 )
    {
        return NULL;
    }

    unsigned int loc = (unsigned int) start_bp;

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
        INDEX_LOC_TYPE* ps_loc = &(genome->index->ps_locs->locs[loc].locs[0]);
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
                   enum CHR_SOURCE chr_source,
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

    /* copy the chr source */
    gen->chr_sources = 
        realloc( gen->chr_sources, sizeof(enum CHR_SOURCE)*(gen->num_chrs));
    (gen->chr_sources)[gen->num_chrs - 1] = chr_source;

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
    struct genome_data* gen, char* fname, enum CHR_SOURCE chr_source
    )
{
    int error;

    /* open the input file */
    FILE* f = fopen( fname, "r" );
    if( NULL == f )
    {
        statmap_log( LOG_FATAL, "Error opening input file %s",  fname  );
        exit( -1 );
    }
    
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
        statmap_log( LOG_FATAL, "Could not allocate enough memory for genome" );
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
                    statmap_log( LOG_NOTICE, "Added '%s' with %i basepairs",  chr_name, i );
                    add_chr_to_genome( chr_name, chr, i, chr_source, gen );
                    
                    if( !feof(f)) {
                        /* reset the bp index to 0 */
                        i = 0;
                    }
                }

                /* move to the next chr index */
                chr_index++;

                /* OVERFLOW BUG - check chr_name */
                /* read in the next string as the chromosome */
                error = fscanf(f, "%254s", chr_name );
                assert( error > 0 );
                while( fgetc(f) != '\n' )
                    ;

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
                statmap_log( LOG_FATAL, "Could not allocate enough memory for genome" );
                assert( 0 );
                exit(-1);
            }
        }        

        if( !isalpha(bp))
        {
            statmap_log( LOG_ERROR, "Unrecognized character '%c' - ignoring",  bp  );
            
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
    add_chr_to_genome( chr_name, chr, i, chr_source, gen );
    
    statmap_log( LOG_NOTICE, "Added '%s' with %i basepairs",  chr_name, i );
    free(chr);

    /* close the input file */
    fclose(f);

    return;
}

/*
 * Iterate through each indexable sequence in the given chromosome
 * between bp's start and stop, adding it to the index.
 */
void
index_contig(
    struct genome_data* genome,
    int seq_len,
    SIGNED_LOC start,
    SIGNED_LOC stop,
    INDEX_LOC_TYPE loc
)
{
    SIGNED_LOC bp_index;

    for( bp_index = start; bp_index < stop; bp_index += 1 )
    {
        /* skip negative bp indices explicitly */
        if( bp_index < 0 )
            continue;

        loc.loc = bp_index;

        /* Add the normal sequence */
        LETTER_TYPE *translation;
        translation = translate_seq(genome->chrs[loc.chr] + bp_index, seq_len);

        /* If we cant add this sequence (probably an N ), continue */
        if( translation == NULL ) {
            continue;
        }

        /* Add the sequence into the tree */
        add_sequence( genome, translation, seq_len, loc );

        free( translation );
    }
}

/* 
 * Return the prefix of chr_name up to the first underscore.
 * NOTE: Caller must free memory allocated by this function
 */
char* get_chr_prefix( char* chr_name )
{
    char* first_underscore = strchr( chr_name, '_' );
    /* TODO: what's the best way to handle this? this is one of our two major assumptions */
    assert( first_underscore != NULL );
    if( first_underscore == NULL )
    {
        statmap_log( LOG_FATAL,
                "Diploid genome chr names must have a prefix delimited by an underscore. Found %s",
                chr_name
            );
    }

    int prefix_len = first_underscore - chr_name;
    char* prefix = calloc( (prefix_len + 1), sizeof(char) ); // calloc for terminating \0
    strncpy( prefix, chr_name, prefix_len );

    return prefix;
}

/*
 * Returns the index into genome->index->map_data for the diploid_map_data_t
 * structure corresponding to the chr specified by chr_index.
 */
int get_map_data_index_from_chr_index(
    struct genome_data* genome,
    int chr_index
)
{
    /* Get the prefix (e.g. chr1 from chr1_paternal) */
    char* prefix = get_chr_prefix( genome->chr_names[chr_index] );

    /* loop over the array of diploid_map_data_t, comparing prefixes */
    int map_data_index;
    for( map_data_index=0;
         map_data_index < genome->index->diploid_maps->count;
         map_data_index++ )
    {
        /* diploid_map_data_t->chr_name is the same as the prefix of a chr_name */
        if( 0 == strncmp( prefix,
                         genome->index->diploid_maps->maps[map_data_index].chr_name,
                         strlen(prefix)
                       ) )
            break;
    }
    free( prefix );

    /* This should never happen - an error would have been raised while building the genome */
    assert( map_data_index < genome->index->diploid_maps->count );
    return map_data_index;
}

void
index_diploid_chrs(
    struct genome_data* genome,
    int chr_index,
    int seq_len,
    INDEX_LOC_TYPE loc
)
{
    /* get paternal and maternal indexes for this chr */
    int paternal_chr_index = chr_index;
    char* prefix = get_chr_prefix( genome->chr_names[chr_index] );
    int maternal_chr_index = find_diploid_chr_index(
            genome, prefix, MATERNAL
        );
    // DEBUG
    char debug_fname[255];
    sprintf( debug_fname, "%s_sequence_segments", prefix );
    free( prefix );

    // DEBUG
    //fprintf( stderr, "Paternal index: %d\n", paternal_chr_index );
    //fprintf( stderr, "Maternal index: %d\n", maternal_chr_index );

    /* get corresponding diploid map data struct index */
    int map_data_index = get_map_data_index_from_chr_index(
            genome, paternal_chr_index
        );

    /* build unique sequence segments */
    struct chr_subregion_t* segments;
    int num_segments;
    build_unique_sequence_segments(
            &(genome->index->diploid_maps->maps[map_data_index]),
            seq_len,
            &segments,
            &num_segments
        );

    /* loop over sequence segments */
    int i;

    // DEBUG
    FILE* debugfp = fopen( debug_fname, "w" );
    for( i=0; i < num_segments; i++ )
    {
        fprintf( debugfp, "Segment #%i: %i, %i, %i\n", i,
                segments[i].paternal_start_pos,
                segments[i].maternal_start_pos,
                segments[i].segment_length
            );
    }
    fclose(debugfp);

    for( i=0; i < num_segments; i++ )
    {

        int start_pos; 

        /*** Set bit flags and chr_index ***/
        /* identical sequence on paternal and maternal */
        if( segments[i].paternal_start_pos > 0 &&
            segments[i].maternal_start_pos > 0 )
        {
            /* set flags */
            loc.is_paternal = 1; loc.is_maternal = 1;
            /* set start_pos */
            start_pos = segments[i].paternal_start_pos;
            /* set chr_index */
            loc.chr = paternal_chr_index;
        }
        /* unique sequence on the paternal */
        else if( segments[i].maternal_start_pos == 0 ) {
            /* set flags */
            loc.is_paternal = 1; loc.is_maternal = 0;
            /* set start_pos */
            start_pos = segments[i].paternal_start_pos;
            /* set chr_index */
            loc.chr = paternal_chr_index;
        }
        /* unique sequence on the maternal */
        else if( segments[i].paternal_start_pos == 0 ) {
            /* set flags */
            loc.is_paternal = 0; loc.is_maternal = 1;
            /* set start pos */
            start_pos = segments[i].maternal_start_pos;
            /* loc.chr defaults to paternal - set to maternal */
            loc.chr = maternal_chr_index;
        }
        else {
            // silence the compiler
            start_pos = -1;
            statmap_log( LOG_FATAL,
                    "Impossible sequence segment : %i, %i, %i",
                    segments[i].paternal_start_pos,
                    segments[i].maternal_start_pos,
                    segments[i].segment_length
                );
        }

        /* index sequence segment */
        index_contig( genome, seq_len,
                start_pos-1, // mappings are 1-indexed, but statmap is 0-indexed
                start_pos-1 + segments[i].segment_length,
                loc
            );
    }
}

void
index_haploid_chr(
    struct genome_data* genome,
    int chr_index,
    int seq_len,
    INDEX_LOC_TYPE loc
)
{
    /* set flags */
    loc.is_paternal = 0; loc.is_maternal = 0;

    /* set chr_index */
    loc.chr = chr_index;

    /* take the whole chr and index as a contig */
    index_contig(
        genome, seq_len,
        0,
        genome->chr_lens[chr_index]-seq_len,
        loc
    );
}

extern void
index_genome( struct genome_data* genome )
{
    /* initialize the constant loc elements */
    INDEX_LOC_TYPE loc;
    
    int seq_len = genome->index->seq_length;
    
    /* 
     * Iterate through each indexable sequence in the genome.
     * We skip chr 0 - the pseudo chromosome.
     */
    int chr_index;
    for( chr_index = 1; chr_index < genome->num_chrs; chr_index++ )
    {
        if( genome->chr_lens[chr_index] > LOCATION_MAX ) {
            statmap_log( LOG_FATAL, "Max indexable chr size is '%i'", LOCATION_MAX );
        }

        statmap_log( LOG_NOTICE, "Indexing '%s'", (genome->chr_names)[chr_index] );

        /* Set the chr index in the soon to be added location */
        loc.chr = chr_index;

        /*** Index chromosome, depending on whether it is diploid or haploid ***/

        /* If chr is diploid, we will index both the paternal chr and its 
         * corresponding maternal chr when we iterate over the paternal,
         * so we skip maternal chrs here. */
        if( genome->chr_sources[chr_index] == PATERNAL )
            index_diploid_chrs( genome, chr_index, seq_len, loc );
        else if ( genome->chr_sources[chr_index] == MATERNAL )
            continue;
        else
            index_haploid_chr( genome, chr_index, seq_len, loc );
    }

    /* sort all of the pseudo locations */
    sort_pseudo_locations( genome->index->ps_locs );

    return;
}
