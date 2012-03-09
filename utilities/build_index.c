/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <libgen.h> // for basename()

#include "../src/index_genome.h"

/* 
   usually, we can just let the OS cleanup the in memory index. 
   However, if we are searching for memory leaks, it may be useful 
   to manually free the index. 
*/
#define DONT_EXPLICITLY_FREE_INDEX

void usage()
{
    fprintf( stderr, "Usage: ./build_index indexed_seq_len output_filename genome.fa(s) [diploid.map(s)]\n" );
    return;
}

struct input_file_group {
    char* prefix;
    int num_files;
    char** filenames;
};

void
init_input_file_group(
    struct input_file_group** ifg,
    char* prefix
)
{
    (*ifg) = malloc( sizeof( struct input_file_group ) );
    assert( (*ifg) != NULL );

    // Allocate memory for, and copy, prefix string
    (*ifg)->prefix = calloc( strlen(prefix) + 1, sizeof(char) );
    assert( (*ifg)->prefix != NULL );
    strcpy( (*ifg)->prefix, prefix );

    (*ifg)->num_files = 0;
    (*ifg)->filenames = NULL;
}

// TODO: test
void
free_input_file_group(
    struct input_file_group* ifg
)
{
    int i;
    for( i=0; i < ifg->num_files; i++ )
    {
        free( ifg->filenames[i] );
    }
    free( ifg->filenames );
    free( ifg->prefix );
    free( ifg );
}

void add_file_to_input_file_group(
    struct input_file_group* ifg,
    char* filename
)
{
    ifg->num_files += 1;

    ifg->filenames = realloc(
        ifg->filenames,
        ifg->num_files * sizeof(char *)
    );
    assert( ifg->filenames != NULL );

    ifg->filenames[ifg->num_files - 1] = calloc(strlen(filename) + 1, sizeof(char) );
    assert( ifg->filenames[ifg->num_files - 1] != NULL );
    strcpy( ifg->filenames[ifg->num_files - 1], filename );
}

void
print_input_file_group(
    struct input_file_group* ifg
)
{
    printf("prefix : %s\n", ifg->prefix);

    int i;
    for( i=0; i < ifg->num_files; i++ )
    {
        printf("%s\n", ifg->filenames[i]);
    }
}

void
group_input_files(
    char** filenames,
    int num_files,
    struct input_file_group** ifgs,
    int* num_ifgs
)
{
    /* loop counters */
    int i, j;
    /* build a list of input file groups based on file prefixes */
    for( i = 0; i < num_files; i++ ) // loop through filenames
    {
        /* find the prefix */
        /* make a copy of the original string, since basename and strtok are destructive */
        char* fncopy = calloc( strlen(filenames[i]) + 1, sizeof(char) );
        strcpy( fncopy, filenames[i] );
        char* bname = basename( fncopy );
        char* prefix = strtok( bname, "_" ); // non-reentrant: uses strtok()! (use strtok_r() instead?)
        assert( prefix != NULL );

        /* loop through ifgs */
        for( j = 0; j < *num_ifgs; j++ )
        {
            /* is the prefix already in the list? */
            if( 0 == strcmp( prefix, ifgs[j]->prefix ) )
            {
                // yes - add it to the ifg
                add_file_to_input_file_group( ifgs[j], filenames[i] );
                break; // continue on the next input filename
            }
        }
        // no - create a new ifg, and add this file to the new group
        if( j == *num_ifgs )
        {
            *num_ifgs += 1;
            // realloc appears to be failing here - not that if realloc fails, it returns NULL
            *ifgs = realloc(
                *ifgs,
                *num_ifgs * sizeof( struct input_file_group * )
            );
            assert( *ifgs != NULL );
            init_input_file_group( &ifgs[*num_ifgs - 1], prefix );
            add_file_to_input_file_group( ifgs[*num_ifgs - 1], filenames[i] );
        }

        free( fncopy );
    }
}

int 
main( int argc, char** argv )
{
    int i; // loop counter

    if( argc < 4 )
    {
        usage();
        exit( 1 );
    }

    int indexed_seq_len = atoi( argv[1] );

    char* output_fname = argv[2];
    char index_fname[500];
    sprintf( index_fname, "%s.index", output_fname );

    /* Are we indexing a diploid or a haploid genome? */
    enum bool diploid = false;
    // are there .map files in the input?
    for( i = 3; i < argc; i++ )
    {
        /* verify that this is a '.map' file ( just search for the suffix ) */
        if( 0 == strcmp( ".map", argv[i] + (strlen(argv[i]) - 4))) {
            diploid = true;
            break;
        }
    }

    /* sort input files into groups by their prefixes */
    int num_groups = 0;
    struct input_file_group* groups = NULL;
    group_input_files( argv+3, argc-3, &groups, &num_groups );
    // debug - print input file groups
    printf("num_groups : %i\n", num_groups );
    for( i = 0; i < num_groups; i++ )
    {
        print_input_file_group( &groups[i] );
    }

    exit(0); // testing

    /* Initialize the genome data structure */
    struct genome_data* genome;
    init_genome( &genome );

    /*** Load the genome ***/
    /* Add the fasta files to the genome */
    for( i = 3; i < argc; i++ )
    {
        char* genome_fname = argv[i];
        fprintf( stderr, "NOTICE      :  Adding '%s'\n", genome_fname );
        FILE* genome_fp = fopen( genome_fname, "r");
        add_chrs_from_fasta_file( genome, genome_fp, REFERENCE );
    }
    
    /* index the genome */
    // TODO: implement diploid genome indexing
    index_genome( genome, indexed_seq_len, NULL );

    /* write the genome to file */
    write_genome_to_disk( genome, output_fname  );
    
    /* write the index to disk */
    build_ondisk_index( genome->index, index_fname  );

    /* Test the code - make sure the two copies are identical
    struct genome_data* genome_copy;
    load_genome_from_disk( &genome_copy, "ODGenome.bin" );
    // write the genome to file
    write_genome_to_disk( genome_copy, "ODGenome_2.bin"  );
    */

#ifdef EXPLICITLY_FREE_INDEX
    free_genome( genome );
#endif

    return 0;
}
