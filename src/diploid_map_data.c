#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h> // for toupper()

#include "config.h"
#include "genome.h"
#include "diploid_map_data.h"

/*
 * Find containing indexed loc to given paternal loc (binary search)
 * Returns index loc where paternal (loc) matches the index (index.loc) exactly,
 * OR returns the previous index. This is the index of the start of the
 * "containing" contig.
 */
static int
search_index( struct loc_and_index_t* index, int len, SIGNED_LOC loc)
{
    int low = 0;
    int high = len;
    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        //printf("%i-%i-%i\t%i\t%i\n", low, mid, high, index[mid].index, len );
        if( index[mid].loc < loc ) {
            low = mid + 1; 
        } else {
            /* can't be high = mid-1: here A[mid] >= value, */
            /* so high can't be < mid if A[mid] == value    */
            high = mid; 
        }
    }

    /* make sure the binary search is working */
    assert( low <= high );
    assert( low <= len );
    assert( low >= 0 );
    
    if( index[low].loc == loc )
        return index[low].index;
    
    assert( low > 0 );
    return index[low-1].index;
}

void
init_diploid_map_data( 
    struct diploid_map_data_t* map_data, 
    char* chr_name, 
    unsigned int* chr_lens )
{
    /* NOTE: assumes memory was allocated in calling function */

    memcpy( map_data->chr_lens, chr_lens, sizeof(unsigned int)*3 );
    
    size_t chr_name_len = strlen(chr_name)+1;
    map_data->chr_name = malloc( chr_name_len  );
    memcpy( map_data->chr_name, chr_name, chr_name_len );
    
    map_data->num_mappings = 0;
    map_data->mappings = NULL;
    map_data->index_len = 0;
    map_data->index = 0;
}

void
free_diploid_map_data( struct diploid_map_data_t* map_data )
{
    free( map_data->mappings );
    free( map_data->chr_name );
    free( map_data->index );
    free( map_data );
}

void
add_diploid_mapping( 
    struct diploid_map_data_t* map_data,
    struct diploid_mapping_t* mapping 
)
{
    map_data->num_mappings += 1;
    map_data->mappings = realloc( 
        map_data->mappings, 
        sizeof(struct diploid_mapping_t)*(map_data->num_mappings)
    );
    
    map_data->mappings[ map_data->num_mappings - 1 ] = *mapping;
    
    return;
}

/* Writes one diploid_map_data_t struct to file */
void
write_diploid_map_data_t_to_file(
    struct diploid_map_data_t* map_data,
    FILE* fp
)
{
    int rv;

    /* write chr_name, prefixed with > */
    rv = fprintf( fp, ">%s\n", map_data->chr_name );
    assert( rv > 0 );

    /* write chr_lens */
    rv = fprintf( fp, "%u %u\n", map_data->chr_lens[0], map_data->chr_lens[1] );
    assert( rv > 0 );

    /* write the number of mappings */
    rv = fprintf( fp, "%zu\n", map_data->num_mappings );
    assert( rv > 0 );

    /* write the mappings */
    unsigned int i; // loop counter
    // num_mappings + 1 to write the final chr_lens 
    for( i = 0; i < map_data->num_mappings + 1; i++ )
    {
        rv = fprintf( fp, "%i %i\n",
                map_data->mappings[i].paternal,
                map_data->mappings[i].maternal );
        assert( rv > 0 );
    }

    /* write the length of the diploid map index */
    rv = fprintf( fp, "%zu\n", map_data->index_len );
    assert( rv > 0 );

    /* write the diploid map index */
    for( i = 0; i < map_data->index_len; i++ )
    {
        rv = fprintf( fp, "%i %i\n",
                map_data->index[i].loc,
                map_data->index[i].index );
        assert( rv > 0 );
    }
}

/* Writes an array of diploid_map_data_t structs to file */
void write_diploid_map_data_to_file(
    struct diploid_map_data_t* map_data,
    int num_chrs,
    FILE* fp
)
{
    int rv;

    /* Set first byte to zero if there's no data to write */
    if( map_data == NULL )
    {
        rv = fprintf( fp, "%i\n", 0 );
        assert( rv > 0 );
        return;
    } else {
        /* Set it to 1 if there is data */
        rv = fprintf( fp, "%i\n", 1 );
        assert( rv > 0 );
    }

    /* Write the number of structs that will be written */
    rv = fprintf( fp, "num_diploid_chrs : %i\n", num_chrs );
    assert( rv > 0 );

    /* Write the diploid_map_data_t structs */
    int i;
    for( i=0; i < num_chrs; i++ )
    {
        write_diploid_map_data_t_to_file( &(map_data[i]), fp );
    }
}

/* Read a single diploid_map_data_t struct from file
 * Leaves fp at the start of the next diploid_map_data_t or at the end of the file */
void
load_diploid_map_data_t_from_file(
    struct diploid_map_data_t* map_data,
    FILE* fp
)
{
    int rv;
    unsigned int i;

    /* load initial info needed to init diploid_map_data_t */
    char chr_name[500];
    rv = fscanf( fp, ">%s\n", chr_name );
    assert( rv == 1 );
    unsigned int chr_lens[2];
    rv = fscanf( fp, "%u %u\n", &chr_lens[0], &chr_lens[1] );
    assert( rv == 2 );

    /* initialize diploid_map_data struct */
    /* NOTE:  memory is already allocated by calling function */
    init_diploid_map_data( map_data, chr_name, chr_lens );

    /* load mappings */
    size_t num_mappings;
    rv = fscanf( fp, "%zu\n", &num_mappings );
    assert( rv == 1 );

    // NOTE: add_diploid_mappings incremenets num_mappings automatically
    // num_mappings + 1 to read the final chr_lens
    for( i=0; i < num_mappings + 1; i++ )
    {
        struct diploid_mapping_t mapping;
        rv = fscanf( fp, "%i %i\n",
                &(mapping.paternal),
                &(mapping.maternal) );
        assert( rv == 2 );
        add_diploid_mapping( map_data, &mapping );
    }

    map_data->num_mappings -= 1; // Corrects for "extra" mapping (the chr lens)

    /* load diploid map index */
    rv = fscanf( fp, "%zu\n", &map_data->index_len );
    assert( rv == 1 );
    assert( map_data->index_len > 0 );

    // Allocate memory for index
    map_data->index = calloc(
        map_data->index_len,
        sizeof( struct loc_and_index_t )
    );
    for( i=0; i < map_data->index_len; i++ )
    {
        rv = fscanf( fp, "%i %i\n",
                &(map_data->index[i].loc),
                &(map_data->index[i].index) );
        assert( rv == 2 );
    }
}

/* Read diploid_map_data_t structs from file, returns number of structs read */
int
load_diploid_map_data_from_file(
    struct diploid_map_data_t** map_data,
    FILE* fp
)
{
    // Initialize array of diploid map data
    *map_data = NULL;

    int rv;

    /* Check if there is diploid map data to read, marked by first byte */
    int marker;
    rv = fscanf( fp, "%i\n", &marker );
    assert( rv == 1 );
    if( marker == 0 )
        return 0;

    /* Read the number of structs in the file */
    int num_structs;
    rv = fscanf( fp, "num_diploid_chrs : %i\n", &num_structs );
    assert( rv == 1 );
    assert( num_structs > 0 );

    /* Allocate memory for the diploid_map_data_t structs */
    (*map_data) = malloc(
        num_structs * sizeof( struct diploid_map_data_t )
    );

    /* Read in the diploid_map_data_t structs */
    int i;
    for( i=0; i < num_structs; i++ )
    {
        load_diploid_map_data_t_from_file( &((*map_data)[i]), fp );
    }

    return num_structs;
}

/*
 * Loop through all mappings from .map file.
 *
 * If both paternal and maternal entries are non-zero, add this mapping
 * to the index. 
 *
 * Builds data->index, storing in each loc_and_index_t entry:
 *      
 *      data->index.loc = paternal loc (value from .map file)
 *      data->index.index = corresponding index in diploid_map_data_t->mappings
 */
void
index_diploid_map_data( struct diploid_map_data_t* data )
{
    size_t i; // index_len and num_mappings are size_t
    
    /* find the number of non zero reference positions */
    data->index_len = 0;
    for( i = 0; i < data->num_mappings; i++ )
    {
        /* we only index rows where both the paternal and maternal are non-zero. 
           If either entries is zero, then we will just put in the correct 
           coordinates for the diverging sequence. 
        */
        if( data->mappings[i].paternal > 0 && data->mappings[i].maternal > 0 )
            data->index_len++;
    }
    
    /* build and populate the index */
    int mapping_counter = 0;
    data->index = calloc( data->index_len, sizeof(struct loc_and_index_t) );
    for( i = 0; i < data->num_mappings; i++ )
    {
        if( data->mappings[i].paternal > 0 && data->mappings[i].maternal > 0 )
        {
            data->index[mapping_counter].loc = data->mappings[i].paternal;
            assert( data->mappings[i].maternal > 0 );
            data->index[mapping_counter].index = i;
            mapping_counter++;
        }
    }    
    
    return;
}

/*
 * Lookup the maternal loc corresponding to a loc on the paternal genome.
 *
 * If paternal_pos is in unique paternal sequence (where there is no corresponding
 * maternal loc), returns -1
 */
int
find_diploid_locations( struct diploid_map_data_t* data, 
                        SIGNED_LOC paternal_pos
)
{
    int index = search_index( data->index, data->index_len, paternal_pos  );
    assert( data->mappings[index].maternal > 0 );
    assert( data->mappings[index].paternal > 0 );

    int maternal_pos = data->mappings[index].maternal + 
                       ( paternal_pos - data->mappings[index].paternal );

    // look ahead to determine the paternal length
    // return -1 if in unique paternal sequence
    if( data->mappings[index+1].maternal == 0 &&
        data->mappings[index+1].paternal != 0 &&
        paternal_pos >= data->mappings[index+1].paternal )
    {
        //printf("Skipping %i\n", paternal_pos);
        return -1;
    }

    return maternal_pos;
}

void
fprintf_diploid_map_data( FILE* stream, struct diploid_map_data_t* map  )
{
    /* print the header */
    fprintf( stream, "#PAT     MAT\n" );
    
    size_t i;
    for( i = 0; i < map->num_mappings; i++ )
    {
        fprintf( stream, "%i\t%i\n", 
                 map->mappings[i].paternal, 
                 map->mappings[i].maternal );
    }
    
    return;
}

/*
 * File opening boilerplate
 */

FILE*
open( char* fname )
{
    FILE* fp = NULL;
    fp = fopen( fname, "r" );
    if( fp == NULL )
    {
        fprintf( stderr, "FATAL         : Can not open '%s'\n", fname );
        exit( 1 );
    }
    return fp;
}

/*
 * Loops through FASTA file; counts characters that are not newlines
 * Skips lines that start with > or ;
 */
int
count_bp_in_fasta_file( FILE* fp )
{
    char c;
    int bp_count = 0;
    rewind( fp ); // don't assume at the beginning
    while( ( c = fgetc( fp ) ) != EOF )
    {
        if( c == '>' || c == ';' )
        {
            // loop until newline; this is a header line or comment
            // and does not count towards the total bp's
            while( ( c = fgetc( fp ) ) != '\n') continue;
            continue;
        }
        if( c != '\n' )
            bp_count++;
    }

    return bp_count;
}

/*
 * Loop through .map file, storing mapping entries.
 * Needs a copy of the genome so it can add entries for SNPs
 */
void
parse_map_file( char* fname, 
                struct diploid_map_data_t* map_data,
                struct genome_data* genome,
                int paternal_chr_index,
                int maternal_chr_index
              )
{
    /* Open .map file */
    FILE* fp = NULL;
    fp = fopen( fname,"r" );
    if( fp == NULL )
    {
        fprintf( stderr, "FATAL          :  Can not open '%s'\n", fname );
        assert( 0 );
        exit( 1 );
    }
    
    /* verify that this is a '.map' file ( just search for the suffix ) */
    if( 0 != strcmp( ".map", fname + (strlen(fname) - 4) ) )
    {
        char* extension = fname + (strlen(fname) - 4);
        fprintf( stderr, "ERROR         :  '%s' doesn't appear to be a map file ( the extension is %.4s, not '.map' )\n", fname, extension );
        return;
    }
    
    /* Initialize the map data structure */
    /* Use string up to first underscore in chr name as diploid chr_name */
    char chr_name[500];
    sscanf( genome->chr_names[paternal_chr_index], "%[^_]", chr_name );

    unsigned int chr_lens[2] = {
        genome->chr_lens[paternal_chr_index],
        genome->chr_lens[maternal_chr_index]
    };
    
    init_diploid_map_data( map_data, chr_name, chr_lens );

    /* move to the beginning of the file */
    fseek( fp, 0, SEEK_SET );
    char buffer[500];    
    char* rv;
    int line_num = 0;
    while( !feof(fp) )
    {
        line_num++;
        rv = fgets( buffer, 500, fp );
        // read the next line
        if( NULL == rv )
            break;
        
        // check for a comment line
        if( buffer[0] == '#' )
            continue;
        
        // Add mapping
        struct diploid_mapping_t mapping;
        sscanf( buffer, "%*i\t%i\t%i", &(mapping.paternal), &(mapping.maternal) );
        add_diploid_mapping( map_data, &mapping );

        // Look ahead (for SNPs) if we're in a shared region
        if( mapping.paternal > 0 && mapping.maternal > 0 )
        {
            // Save current position in file
            long current_pos = ftell( fp );
            int nextp, nextm;
            // loop until the next paternal mapping
            while( 1 )
            {
                rv = fgets( buffer, 500, fp );
                if( NULL == rv ) {
                    // we were on the last line; use chr_lens
                    nextp = chr_lens[0]; nextm = chr_lens[1];
                    break;
                } else {
                    // read next .map entry; if paternal is non-zero, this is the upper bound
                    sscanf( buffer, "%*i\t%i\t%i", &nextp, &nextm );
                    if( nextp != 0 )
                        break;
                }
            }

            // loop through bps in genome
            int i;
            for( i = 0; i < nextp - mapping.paternal; i++ )
            {
                // - 1 because .map file is 1-indexed, but statmap is 0-indexed
                char *p_ptr = find_seq_ptr( genome, paternal_chr_index, mapping.paternal + i - 1, 1 );
                char *m_ptr = find_seq_ptr( genome, maternal_chr_index, mapping.maternal + i - 1, 1 );
                // test for mismatch (SNP); normalize case of bp letter
                if( toupper(*p_ptr) != toupper(*m_ptr) )
                {
                    // add mappings for SNP
                    struct diploid_mapping_t tmp_mapping;
                    tmp_mapping.paternal = mapping.paternal + i; tmp_mapping.maternal = 0;
                    add_diploid_mapping( map_data, &tmp_mapping );
                    tmp_mapping.paternal = 0; tmp_mapping.maternal = mapping.maternal + i;
                    // end SNP, resuming contig
                    tmp_mapping.paternal = mapping.paternal + i + 1;
                    tmp_mapping.maternal = mapping.maternal + i + 1;
                    add_diploid_mapping( map_data, &tmp_mapping );
                }
            }

            // restore position in file for next mapping
            fseek( fp, current_pos, SEEK_SET );
        }
    }

    /* finally, add a diploid mapping onto the end that includes the chr
       stops. This simplifies logic in that w can always 'peek ahead' */
    struct diploid_mapping_t mapping = { chr_lens[0], chr_lens[1] };
    add_diploid_mapping( map_data, &mapping );
    /* decrement the number of mappings */
    map_data->num_mappings -= 1;

    fclose( fp );
    return;
}

void
build_unique_sequence_segments( struct diploid_map_data_t* data, 
                                int seq_len,
                                struct chr_subregion_t** segments,
                                int* num_segments 
    )
{
    assert( seq_len > 0 );

    /* Allocate memory for segments; one segment for each mapping */
    *segments = calloc( sizeof( struct chr_subregion_t ), data->num_mappings );
    assert( *segments != NULL );
    
    *num_segments = 0;
    size_t i;
    for( i = 0; i < data->num_mappings; i++ )
    {
        int maternal_start = 0;
        int paternal_start = 0;

        int maternal_length = 0;
        int paternal_length = 0;        
        
        int case_code = 2*(int)( data->mappings[i].paternal > 0 );
        case_code += 1*(int)( data->mappings[i].maternal > 0 );

        switch ( case_code )
        {
        case 0:
            continue;
        case 1:
            maternal_start = MAX( 1, data->mappings[i].maternal - seq_len + 1 );
            break;
        case 2:
            paternal_start = MAX( 1, data->mappings[i].paternal - seq_len + 1 );
            break;
        case 3:
            maternal_start = MAX( 1, data->mappings[i].maternal - seq_len + 1 );
            paternal_start = MAX( 1, data->mappings[i].paternal - seq_len + 1 );
            break;
        default:
            fprintf( stderr, 
                     "FATAL         :  impossible case value '%i'\n", 
                     case_code );
            abort();
        }
        
        /* the continue in the case statment should have prevented this,
           but I'm worry because that's a pretty weird construct */
        assert( maternal_start > 0 || paternal_start > 0 );
        
        /* determine the region length */
        /* we can peak ahead one because we added the chr lengths to the end */
        size_t j;
        
        if( maternal_start > 0 )
        {
            for( j = i+1; data->mappings[j].maternal == 0; j++ )
                assert( j <= data->num_mappings );
            maternal_length = data->mappings[j].maternal - maternal_start - seq_len + 1;
        }

        if( paternal_start > 0 )
        {
            for( j = i+1; data->mappings[j].paternal == 0; j++ )
                assert( j <= data->num_mappings );
            paternal_length = data->mappings[j].paternal - paternal_start - seq_len + 1;
        }
        
        assert( maternal_length == 0 || paternal_length == 0
                || maternal_length == paternal_length || i == data->num_mappings-1 );
        
        /* the map files are 1 indexed, but statmap uses zero based locations. */
        paternal_start -= 1;
        maternal_start -= 1;
        
          // DEBUG
#if 0
        fprintf( stderr, "%i-%i\t%i-%i\t%i, %i\n", 
                 paternal_start, paternal_start+paternal_length, 
                 maternal_start, maternal_start+maternal_length,
                 paternal_length, maternal_length );
#endif

        struct chr_subregion_t subregion =  { 
            paternal_start, maternal_start, MAX(paternal_length, maternal_length)
        };
        (*segments)[*num_segments] = subregion;
        *num_segments += 1;
    }
    
    *segments = realloc( *segments, sizeof( struct chr_subregion_t )*(*num_segments) );
    /* This should never fail because I'm always making the allocation smaller */
    assert( *segments != NULL );
    
    
    return;
}

#if 0
int 
main( int argc, char** argv )
{
    struct genome_data* genome;
    init_genome( &genome );
    add_chrs_from_fasta_file( genome, fopen( argv[2], "r" ), PATERNAL );
    add_chrs_from_fasta_file( genome, fopen( argv[3], "r" ), MATERNAL );

    struct diploid_map_data_t* map_data;
    parse_map_file( argv[1], &map_data, genome, 1, 2 );
    if( NULL == map_data )
    {
        fprintf( stderr, "fatal         :  can not parse map file.\n" );
        return -1;
    }
    

    index_diploid_map_data( map_data );

    // test writing and reading from file
    // write map_data to file
    FILE* ofp = fopen( "test.dmap", "r" );
    write_diploid_map_data_to_file( map_data, 1, ofp );
    fclose( ofp );

    // read map_data from saved file
    struct diploid_map_data_t* input_map_data;
    FILE* ifp = fopen( "test.dmap", "r" );
    int num_structs = read_diploid_map_data_from_file( &input_map_data, ifp );
    fclose( ifp );

    // write it out again; then diff the output files
    FILE* check_fp = fopen( "check_test.dmap", "r" );
    write_diploid_map_data_to_file( input_map_data, 1, check_fp );
    fclose( check_fp );

    struct chr_subregion_t* segments;
    int num_segments;
    int seq_len = 10;
    build_unique_sequence_segments( map_data, seq_len, &segments, &num_segments );
    
    struct genome_data* genome;
    init_genome( &genome );
    add_chrs_from_fasta_file( genome, fopen( argv[2], "r" ) );
    add_chrs_from_fasta_file( genome, fopen( argv[3], "r" ) );

    return 0;
}
#endif
