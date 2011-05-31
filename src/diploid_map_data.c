#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "genome.h"

/*

TODO - get rid of the reference part. It doesn't make sense to
keep it because we will never be mapping to it. Instead, just
use the paternal as the refernce and modify the code to convert
to the maternal. A couple points:

1) We can skip any row that has zeros for both paternal and maternal
2) We will, by default, use the paternal chromsome coordinate in cases 
   with identical coordinates.







 */


struct loc_and_index_t {
    SIGNED_LOC loc;
    int index;
};

static int
search_index( struct loc_and_index_t* index, int len, SIGNED_LOC loc)
{
    int low = 0;
    int high = len;
    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        // printf("%i-%i-%i\t%i\t%i\n", low, mid, high, index[mid].index, len );
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

struct diploid_mapping_t {
    SIGNED_LOC paternal;
    SIGNED_LOC maternal;
};

struct diploid_map_data_t {
    char* chr_name;
    unsigned int chr_lens[2];
    size_t num_mappings;
    struct diploid_mapping_t* mappings;
    size_t index_len;
    struct loc_and_index_t* index;
};

void
init_diploid_map_data( 
    struct diploid_map_data_t** map_data, 
    char* chr_name, 
    unsigned int* chr_lens )
{
    (*map_data) = malloc( sizeof( struct diploid_map_data_t ) );
    
    memcpy( (*map_data)->chr_lens, chr_lens, sizeof(unsigned int)*3 );
    
    size_t chr_name_len = strlen(chr_name)+1;
    (*map_data)->chr_name = malloc( chr_name_len  );
    memcpy( (*map_data)->chr_name, chr_name, chr_name_len );
    
    (*map_data)->num_mappings = 0;
    (*map_data)->mappings = NULL;
    (*map_data)->index_len = 0;
    (*map_data)->index = 0;
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

void
index_diploid_map_data( struct diploid_map_data_t* data )
{
    /* generic loop counter */
    int i;
    
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

void
find_diploid_locations( struct diploid_map_data_t* data, 
                        SIGNED_LOC paternal_pos,
                        SIGNED_LOC* maternal_pos
)
{
    int index = search_index( data->index, data->index_len, paternal_pos  );
    assert( data->mappings[index].maternal > 0 );
    assert( data->mappings[index].paternal > 0 );
    *maternal_pos = data->mappings[index].maternal + 
                       ( paternal_pos - data->mappings[index].paternal );
    return;
}

void
fprintf_diploid_map_data( FILE* stream, struct diploid_map_data_t* map  )
{
    /* print the header */
    fprintf( stream, "#PAT     MAT\n" );
    
    int i;
    for( i = 0; i < map->num_mappings; i++ )
    {
        fprintf( stream, "%i\t%i\n", 
                 map->mappings[i].paternal, 
                 map->mappings[i].maternal );
    }
    
    return;
}

void
parse_map_file( char* fname, struct diploid_map_data_t** map_data  )
{
    /* Initialize map data */
    *map_data = NULL;
    
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
    
    /* find the chromosome name */
    /* we assume the chromosome name is the string characters until 
       the first underscore or period */
    
    unsigned int chr_lens[2] = {135262207+1000, 135275611+1000};
    
    init_diploid_map_data( map_data, "test_chr", chr_lens );

    /* move to the beginning of the file */
    fseek( fp, 0, SEEK_SET );
    char buffer[500];    
    int line_num = 0;
    while( !feof(fp) )
    {
        line_num++;
        char* rv = fgets( buffer, 500, fp );
        // read the next line
        if( NULL == rv )
            break;
        
        // check for a comment line
        if( buffer[0] == '#' )
            continue;
        
        struct diploid_mapping_t mapping;
        
        sscanf( buffer, "%*i\t%i\t%i", &(mapping.paternal), &(mapping.maternal) );
        
        add_diploid_mapping( *map_data, &mapping );
    }

    /* finally, add a diploid mapping onto the end that includes the chr
       stops. This simplifies logic in that w can always 'peak ahead' */
    struct diploid_mapping_t mapping = { chr_lens[0], chr_lens[1] };
    add_diploid_mapping( *map_data, &mapping );
    /* decrement the number of mappings */
    (*map_data)->num_mappings -= 1;

    return;
}

struct chr_subregion_t {
    SIGNED_LOC paternal_start_pos;
    SIGNED_LOC maternal_start_pos;
    SIGNED_LOC segment_length;
};

void
build_unique_sequence_segments( struct diploid_map_data_t* data, 
                                int seq_len,
                                struct chr_subregion_t** segments,
                                int* num_segments 
    )
{
    assert( seq_len > 0 );

    *segments = calloc( sizeof( struct chr_subregion_t ), data->num_mappings );
    assert( *segments != NULL );
    
    *num_segments = 0;
    int i;
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
        int j;
        
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
                || maternal_length == paternal_length        );
        
        /* the map files are 1 indexed, but statmap uses zero based locations. */
        paternal_start -= 1;
        maternal_start -= 1;
        
        /*
          // DEBUG
        fprintf( stderr, "%i-%i\t%i-%i\n", 
                 paternal_start, paternal_start+paternal_length, 
                 maternal_start, maternal_start+maternal_length );
        */
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


int 
main( int argc, char** argv )
{
    struct diploid_map_data_t* map_data;
    parse_map_file( argv[1], &map_data  );
    if( NULL == map_data )
    {
        fprintf( stderr, "FATAL         :  Can not parse map file.\n" );
        return -1;
    }

    index_diploid_map_data( map_data );

    struct chr_subregion_t* segments;
    int num_segments;
    build_unique_sequence_segments( map_data, 10, &segments, &num_segments );
    
    struct genome_data* genome;
    init_genome( &genome );
    add_chrs_from_fasta_file( genome, fopen( "chr2_paternal.fa", "r")  );
    add_chrs_from_fasta_file( genome, fopen( "chr2_maternal.fa", "r")  );
    
    /* make sure the inverse is working */
    #if 0
    int i;
    SIGNED_LOC prev_maternal = 0;
    for( i = 87670; i < 135372312+10000; i++ )
    {
        SIGNED_LOC mat, pat;
        pat = i;
        find_diploid_locations( map_data, pat, &mat );
        fprintf( stderr, "%i\t%i\n", mat, prev_maternal );
        prev_maternal = mat;
    }
    #endif
    
    return 0;
}

