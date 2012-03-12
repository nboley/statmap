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

struct file_group {
    char* prefix;
    int num_files;
    char** filenames;
};

struct file_group_list {
    int num_groups;
    struct file_group* groups;
};

void
init_file_group (
    struct file_group* fg,
    char* prefix
)
{
    /* NOTE: assumes memory for file_group struct fg has been allocated
     * (by realloc in add_group_to_file_group_list) */

    // allocate memory for, and copy, prefix
    fg->prefix = calloc( strlen(prefix) + 1, sizeof(char) );
    strcpy( fg->prefix, prefix );

    fg->num_files = 0;
    fg->filenames = NULL;
}

void
free_file_group(
    struct file_group* fg
)
{
    // free all allocated filenames
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        free( fg->filenames[i] );
    }
    
    free( fg->filenames );
    free( fg->prefix );
}

void
init_file_group_list(
    struct file_group_list** fgl
)
{
    (*fgl) = malloc( sizeof( struct file_group_list ) );
    
    (*fgl)->num_groups = 0;
    (*fgl)->groups = NULL;
}

void
free_file_group_list(
    struct file_group_list* fgl
)
{
    // free memory allocated in file groups
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        free_file_group( &(fgl->groups[i]) );
    }

    // free contigious block of file_group structs
    free( fgl->groups );
    free( fgl );
}

void
add_filename_to_file_group(
    struct file_group* fg,
    char* filename
)
{
    fg->num_files += 1;

    // Allocate a pointer to char for filename
    fg->filenames = realloc(
        fg->filenames,
        fg->num_files * sizeof(char *)
    );
    // Allocate space for the filename
    fg->filenames[fg->num_files - 1] = calloc( strlen(filename) + 1, sizeof(char) );
    // Copy filename
    strcpy( fg->filenames[fg->num_files - 1], filename );
}

void
add_group_to_file_group_list(
    char* prefix,
    struct file_group_list* fgl
)
{
    fgl->num_groups += 1;

    // Allocate memory for new file group
    fgl->groups = realloc(
        fgl->groups,
        fgl->num_groups * sizeof( struct file_group )
    );

    // Initialize new file group
    init_file_group( &(fgl->groups[fgl->num_groups - 1]), prefix );
}

void
print_file_group(
    struct file_group* fg
)
{
    // print the prefix
    printf("prefix : %s\n", fg->prefix);
    // print the filenames
    int i;
    for( i=0; i < fg->num_files; i++ )
    {
        printf("%s\n", fg->filenames[i] );
    }
}

void
print_file_group_list(
    struct file_group_list* fgl
)
{
    // print the file groups
    int i;
    for( i=0; i < fgl->num_groups; i++ )
    {
        print_file_group( &(fgl->groups[i]) );
    }
}

int
find_prefix_index_in_file_group_list(
    struct file_group_list* fgl,
    char* prefix
)
{
    int prefix_index;
    /* loop through files in group */
    for( prefix_index=0; prefix_index < fgl->num_groups; prefix_index++ )
    {
        if( 0 == strcmp( fgl->groups[prefix_index].prefix, prefix ) )
            return prefix_index;
    }

    /* if prefix not found in list of file groups, return -1 */
    return -1;
}

void
group_files(
    char** filenames,
    int num_files,
    struct file_group_list* fgl
)
{
    int i;
    /* loop through filenames */
    for( i = 0; i < num_files; i++ )
    {
        /* find the prefix */
        /* make a copy of the original string, since basename and strtok are destructive */
        char* fncopy = calloc( strlen(filenames[i]) + 1, sizeof(char) );
        strcpy( fncopy, filenames[i] );
        char* bname = basename( fncopy );
        char* prefix = strtok( bname, "_" ); // non-reentrant: uses strtok()! (use strtok_r() instead?)

        int prefix_index = find_prefix_index_in_file_group_list( fgl, prefix );

        /* if we didn't find the prefix in the list, 
         * add a new file group to the list */
        if( prefix_index == -1 )
        {
            add_group_to_file_group_list( prefix, fgl );
            prefix_index = fgl->num_groups - 1; // set prefix_index to newly added file group
        }

        /* add the filename to the file group */
        add_filename_to_file_group( &(fgl->groups[prefix_index]), filenames[i] );

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
    struct file_group_list* fgl = NULL;
    init_file_group_list( &fgl );
    group_files( argv+3, argc-3, fgl );
    // debug - print file group list
    print_file_group_list( fgl );
    // free file group list
    free_file_group_list( fgl );

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
