#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"

struct loc_and_index_t {
    SIGNED_LOC loc;
    int index;
};

int
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
        return low;
    
    assert( low > 0 );
    return low - 1;
}

struct diploid_mapping_t {
    SIGNED_LOC ref;
    SIGNED_LOC paternal;
    SIGNED_LOC maternal;
};

struct diploid_map_data_t {
    char* chr_name;
    size_t num_mappings;
    struct diploid_mapping_t* mappings;
    size_t index_len;
    struct loc_and_index_t* index;
};

void
init_diploid_map_data( struct diploid_map_data_t** map_data, char* chr_name )
{
    (*map_data) = malloc( sizeof( struct diploid_map_data_t ) );
    
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
        if( data->mappings[i].ref > 0 )
            data->index_len++;
    }
    
    /* build and populate the index */
    int mapping_counter = 0;
    data->index = calloc( data->index_len, sizeof(struct loc_and_index_t) );
    for( i = 0; i < data->num_mappings; i++ )
    {
        if( data->mappings[i].ref > 0 )
        {
            data->index[mapping_counter].loc = data->mappings[i].ref;
            data->index[mapping_counter].index = i;
            mapping_counter++;
        }
    }    
    
    return;
}

struct diploid_mapping_t*
find_diploid_mapping( struct diploid_map_data_t* data, SIGNED_LOC ref_pos )
{
    /* since the comparison only acres about the ref position,
       we ignore the other fields */
    struct loc_and_index_t key = {0,0};
    key.loc = ref_pos;
    int index = search_index( data->index, data->index_len, ref_pos  );
    assert( data->mappings[index].ref <= ref_pos );
    return data->mappings + index;
}

void
find_diploid_locations( struct diploid_map_data_t* data, 
                        SIGNED_LOC ref_pos,
                        SIGNED_LOC* paternal_pos,
                        SIGNED_LOC* maternal_pos
)
{
    *paternal_pos = -1;
    *maternal_pos = -1;
    
    struct diploid_mapping_t* mapping
        = find_diploid_mapping( data, ref_pos );
    
    if( mapping->ref == 0 )
        return;
    
    if( mapping->paternal > 0 ) {
        *paternal_pos = mapping->paternal; 
        *paternal_pos += ( ref_pos - mapping->ref ); 
    }

    if( mapping->maternal > 0 ) {  
        *maternal_pos = mapping->maternal;
        *maternal_pos += ( ref_pos - mapping->ref ); 
    }
    
    return;
}

void
fprintf_diploid_map_data( FILE* stream, struct diploid_map_data_t* map  )
{
    /* print the header */
    fprintf( stream, "#REF    PAT     MAT\n" );
    
    int i;
    for( i = 0; i < map->num_mappings; i++ )
    {
        fprintf( stream, "%i\t%i\t%i\n", 
                 map->mappings[i].ref, 
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
    
    
    init_diploid_map_data( map_data, "test_chr" );

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
        
        sscanf( buffer, "%i\t%i\t%i", 
                &(mapping.ref), &(mapping.paternal), &(mapping.maternal) );
        
        add_diploid_mapping( *map_data, &mapping );
    }

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

    int i;
    SIGNED_LOC prev_paternal = 0;
    SIGNED_LOC prev_maternal = 0;
    for( i = 132483510; i < 135372312+10000; i++ )
    {
        SIGNED_LOC mat, pat;
        find_diploid_locations( map_data, i, &mat, &pat );
        fprintf( stderr, "%i\t%i:%i\t%i:%i\n", i, mat, prev_maternal, pat, prev_paternal );
        assert( pat == -1 || pat > prev_paternal );
        prev_paternal = pat;
        assert( mat == -1 || mat > prev_maternal );
        prev_maternal = mat;
    }
    
    return 0;
}

