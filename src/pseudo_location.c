#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "pseudo_location.h"

static void* 
realloc_CE( void* ptr, size_t size )
{
    ptr = realloc( ptr, size );
    if( NULL == ptr )
    {
        fprintf( stderr, "FATAL       : Failed to alloc space for pseudo location.\n");
        assert( false );
        exit( -1 );
    }
    
    return ptr;
}

int
cmp_genome_location( void* loc1, void* loc2 )
{
    /* first test the chromosome identifier */
    if( ((GENOME_LOC_TYPE*) loc1)->chr
        > ((GENOME_LOC_TYPE*) loc2)->chr )
    {
        return 1;
    }
  
    if( ((GENOME_LOC_TYPE*) loc1)->chr 
        < ((GENOME_LOC_TYPE*) loc2)->chr )
    {
        return -1;
    }

    /* since the chromosomes must be identical... */
    if( ((GENOME_LOC_TYPE*) loc1)->loc
        > ((GENOME_LOC_TYPE*) loc2)->loc )
    {
        return 1;
    }
  
    if( ((GENOME_LOC_TYPE*) loc1)->loc 
        < ((GENOME_LOC_TYPE*) loc2)->loc )
    {
        return -1;
    }

    
    return 0;
}

/*
 * Pseudo location
 *
 * Each pseudo location corresponds to multiple real locations
 * with the same sequence. For very heavily repeated sequence, 
 * we use this to cut down on the storage and processing costs.
 *
 */

/* We assume that space has already been allocated for the struct */
static  void
init_pseudo_location( struct pseudo_location_t* loc )
{
    loc->num = 0;
    loc->locs = NULL;
}

static  void 
free_pseudo_location( struct pseudo_location_t* loc )
{
    if( NULL != loc->locs )
        free( loc->locs );

    return;
}

static  void
sort_pseudo_location( struct pseudo_location_t* loc )
{
    qsort( loc->locs, 
           loc->num, 
           sizeof(GENOME_LOC_TYPE), 
           (int(*)( const void*, const void* ))cmp_genome_location 
    );

    return;
}

 void
add_loc_to_pseudo_location( 
    struct pseudo_location_t* ps_loc, const GENOME_LOC_TYPE* const loc)
{
    ps_loc->num += 1;
    ps_loc->locs = realloc_CE( ps_loc->locs, ps_loc->num*sizeof(GENOME_LOC_TYPE) );
    ps_loc->locs[ps_loc->num-1] = *loc;
};

void
fprint_pseudo_locations( FILE* of, struct pseudo_locations_t* ps_locs )
{
    int i;
    for( i = 0; i < ps_locs->num; i++ )
    {
        assert( ps_locs->locs[i].num > 0 );
        fprintf( of, "%i\t%i", i, ps_locs->locs[i].num );
        int j;
        for( j = 0; j < ps_locs->locs[i].num; j++ )
        {
            fprintf( of, 
                "\t%i,%i", 
                ps_locs->locs[i].locs[j].chr,
                ps_locs->locs[i].locs[j].loc 
            );
        }
        fprintf( of, "\n" );
    }
}

size_t
write_pseudo_locations_to_file( struct pseudo_locations_t* ps_locs, FILE* of )
{
    int rv;    
    size_t size_written = 0;
    
    if( NULL == ps_locs )
    {
        int zero = 0;
        rv = fwrite( &zero, sizeof(int), 1, of );
        assert( rv == 1 );
        return sizeof( int );
    }
    
    /* write the size of the file. For now, we just put in a place holder
       and, after we write everything, we will go back and fill it in */
    rv = fwrite( &(size_written), sizeof(size_t), 1, of );
    assert( rv == 1 );
    size_written += sizeof( size_t );

    /* write the number of pseudo locations */    
    rv = fwrite( &(ps_locs->num), sizeof(int), 1, of );
    assert( rv == 1 );
    size_written += sizeof( int );
    
    int i;
    for( i = 0; i < ps_locs->num; i++ )
    {
        rv = fwrite( &(ps_locs->locs[i].num), sizeof(int), 1, of );
        assert( rv == 1 );
        size_written += sizeof( int );
        
        rv = fwrite( ps_locs->locs[i].locs, sizeof(GENOME_LOC_TYPE), ps_locs->locs[i].num, of );
        assert( rv == ps_locs->locs[i].num );
        size_written += ps_locs->locs[i].num*sizeof( GENOME_LOC_TYPE );
    }

    /* seek back to the beggining and update the size written */
    fseek( of, 0, SEEK_SET );
    rv = fwrite( &(size_written), sizeof(size_t), 1, of );
    assert( rv == 1 );
    
    return size_written;
}

void
load_pseudo_locations_from_mmapped_data( 
    struct pseudo_locations_t** ps_locs, char* data )
{
    int num = *data;
    /* if there are no pseudo locations, dont even init */
    if( num == 0 )
    {
        *ps_locs = NULL;
        return;
    }

    init_pseudo_locations( ps_locs );

    data += sizeof(size_t);
    
    (*ps_locs)->num = *( (int*) data );
    data += sizeof(int);
        
    (*ps_locs)->locs = malloc( (*ps_locs)->num*sizeof(struct pseudo_location_t) );
    
    int i, j;
    for( i = 0; i < (*ps_locs)->num; i++ )
    {
        (*ps_locs)->locs[i].num = *( (int*) data );
        data += sizeof(int);

        // allocate memory for GENOME_LOC_TYPEs
        (*ps_locs)->locs[i].locs = malloc( (*ps_locs)->locs[i].num * sizeof( GENOME_LOC_TYPE ) );
        // copy GENOME_LOC_TYPEs
        for( j = 0; j < (*ps_locs)->locs[i].num; j++ )
        {
            memcpy( &((*ps_locs)->locs[i].locs[j]), data, sizeof( GENOME_LOC_TYPE ) );
            data += sizeof( GENOME_LOC_TYPE );
        }
    }

    return;
}

void
load_pseudo_locations( 
    FILE* fp, struct pseudo_locations_t** ps_locs )
{
    size_t size = 0;
    size_t rv = 0;
    fread( &size, sizeof(size_t), 1, fp );
    assert( rv == 0 );
    
    char* data = malloc( size  );
    assert( data != NULL );
    fseek( fp, 0, SEEK_SET );
    rv = fread( data, 1, size, fp );
    assert( rv == size );
    
    load_pseudo_locations_from_mmapped_data( ps_locs, data );

    free( data );

    return;
}


void
load_pseudo_locations_from_text( FILE* fp, struct pseudo_locations_t** ps_locs )
{
    int error;
    
    /* initialize the container */
    init_pseudo_locations( ps_locs );

    if( fp == NULL )
        return;
    
    while( !feof( fp ) )
    {
        
        /* initialize a new pseudo location */
        struct pseudo_location_t ps_loc;
        init_pseudo_location( &ps_loc );
        
        int chr_index = -1; 
        int num_locs = -1;
        
        /* read the line header ( chromosome and number ) */
        error = fscanf( fp, "%i\t%i", &chr_index, &num_locs );
        if( error != 2 )
        {
            if( feof( fp ) )
            {
                /* gracefuly alowing spurious trailing newlines */
                break;
            }
            fprintf( stderr, "FATAL         : Error reading line header in pseudo loc parsing code.\n" );
            exit( 1 );
        }
        
        int i;
        for( i = 0; i < num_locs; i++ )
        {
            int chr = -1;
            int pos = -1;
            error = fscanf( fp, "\t%i,%i", &chr, &pos );
            if( error != 2 )
            {
                fprintf( stderr, "FATAL         : Error reading data in pseudo loc parsing code.\n" );
                exit( 1 );
            }
            
            GENOME_LOC_TYPE loc;
            memset( &loc, 0, sizeof( GENOME_LOC_TYPE ) );
            loc.chr = chr;
            loc.loc = pos;
            add_loc_to_pseudo_location( &ps_loc, &loc );
        }
        
        int index = add_new_pseudo_location( *ps_locs );
        (*ps_locs)->locs[index] = ps_loc;

        /* move the pointer to the next line. If there is nothing, then check for eof. */
        char nl = getc(fp);
        if( nl != '\n' )
        {
            assert( feof( fp ) );
        }

    }

}

/*
 * Pseudo Locations 
 *
 * container for pseudo location's
 * 
 */

void
init_pseudo_locations( 
    struct pseudo_locations_t** locs )
{
    *locs = realloc_CE( NULL, sizeof(struct pseudo_locations_t) );
    (*locs)->num = 0;
    (*locs)->locs = NULL;
}
    
void
free_pseudo_locations( struct pseudo_locations_t* locs )
{
    if( NULL != locs->locs )
    {
        int i;
        for( i = 0; i < locs->num; i++ )
            free_pseudo_location( locs->locs + i );

        free( locs->locs );
    }
    
    free( locs );
}

void
sort_pseudo_locations( struct pseudo_locations_t* locs )
{
    int i;
    for( i = 0; i < locs->num; i++ )
    {
        sort_pseudo_location( locs->locs + i );
    }
    
    return;
}


/* returns the index of the new entry */
unsigned int 
add_new_pseudo_location( struct pseudo_locations_t* locs )
{
    locs->num += 1;
    locs->locs = realloc( locs->locs, locs->num*sizeof(struct pseudo_location_t) );
    init_pseudo_location( locs->locs + locs->num - 1 );

    return locs->num - 1;
}

