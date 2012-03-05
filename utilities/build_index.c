/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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

int 
main( int argc, char** argv )
{
    if( argc < 4 )
    {
        usage();
        exit( 1 );
    }
    
    int indexed_seq_len = atoi( argv[1] );

    char* output_fname = argv[2];
    char index_fname[500];
    sprintf( index_fname, "%s.index", output_fname );

    /* Initialize the genome data structure */
    struct genome_data* genome;
    init_genome( &genome );

    /*** Load the genome ***/
    /* Add the fasta files to the genome */
    int i;
    for( i = 3; i < argc; i++ )
    {
        char* genome_fname = argv[i];
        fprintf( stderr, "NOTICE      :  Adding '%s'\n", genome_fname );
        FILE* genome_fp = fopen( genome_fname, "r");
        add_chrs_from_fasta_file( genome, genome_fp, REFERENCE );
    }
    
    /* index the genome */
    index_genome( genome, indexed_seq_len );

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
