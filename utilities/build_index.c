/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "../src/index_genome.h"

/* 
   usually, we can just let the OS cleanup the in memory index. 
   However, if we are searching for memory leaks, it may be useful 
   to manually free the index. 
*/
#define DONT_EXPLICITLY_FREE_INDEX

void usage()
{
    fprintf( stderr, "Usage: ./build_index indexed_seq_len output_filename genome.fa(s)\n" );
    return;
}

struct input_file_group {
    char* prefix;
    int num_files;
    char** files;
};

void
init_input_file_group(
    struct input_file_group** ifg,
    char* prefix
)
{
    (*ifg) = malloc( sizeof( struct input_file_group ) );
    (*ifg)->num_files = 0;
    (*ifg)->files = NULL;

    // Init and copy prefix string
    (*ifg)->prefix = calloc( strlen(prefix) + 1, sizeof(char) );
    assert( (*ifg)->prefix != NULL );
    strcpy( (*ifg)->prefix, prefix );
}

void add_file_to_input_file_group(
    struct input_file_group** ifg,
    char* file
)
{
    (*ifg)->num_files += 1;
    (*ifg)->files = realloc(
        (*ifg)->files,
        (*ifg)->num_files * sizeof(char *)
    );
    (*ifg)->files[(*ifg)->num_files - 1] = calloc( strlen(file) + 1, sizeof(char) );
    strcpy( (*ifg)->files[(*ifg)->num_files - 1], file );
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
        printf("%s\n", ifg->files[i]);
    }
}

void
group_input_files(
    int argc,
    char** argv,
    struct input_file_group** ifgs,
    int* num_ifgs
)
{
    int i, j;
    /* build a list of input file groups based on file prefixes */
    for( i = 3; i < argc; i++ ) // loop through filenames in argv
    {
        /* find the prefix */
        // first, pull out the filename
        // start at the end of the string, step back until slash
        assert( strlen( argv[i] ) > 1 );
        int filename_start;
        for( filename_start = strlen( argv[i] ) ;
             filename_start > 1 && argv[i][filename_start-1] != '/' ;
             filename_start-- )
            ;
        // move forward until we encounter an underscore
        int prefix_end;
        for( prefix_end = filename_start ;
             prefix_end < strlen( argv[i] ) && argv[i][prefix_end] != '_' ;
             prefix_end++ )
            ;
        int prefix_len = prefix_end - filename_start;
        char prefix[prefix_len + 1];
        strncpy( prefix, argv[i] + filename_start, prefix_len );
        prefix[prefix_len] = '\0'; // append null byte

        // is this prefix already in ifgs?
        for( j = 0; j < *num_ifgs; j++ )
        {
            if( 0 == strcmp( prefix, ifgs[j]->prefix ) )
            {
                // yes - add it to the ifg
                add_file_to_input_file_group( &ifgs[j], argv[i] );
                break; // continue on the next input filename
            }
        }
        if( j == *num_ifgs )
        {
            // no - create a new ifg, and add this file to the new group
            *num_ifgs += 1;
            ifgs = realloc(
                ifgs,
                *num_ifgs * sizeof( struct input_file_group * )
            );
            init_input_file_group( &ifgs[*num_ifgs - 1], prefix );
            add_file_to_input_file_group( &ifgs[*num_ifgs - 1], argv[i] );
        }
    }
}

int 
main( int argc, char** argv )
{
    int i, j; // loop counter

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
    group_input_files( argc, argv, &groups, &num_groups );
    // debug - print input file groups
    printf("num_groups : %i\n", num_groups );
    for( i = 0; i < num_groups; i++ )
    {
        print_input_file_group( &groups[i] );
    }

    exit(0);

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
