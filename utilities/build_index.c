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
#define EXPLICITLY_FREE_INDEX

void usage()
{
    fprintf( stderr, "Usage: ./build_index genome.fa indexed_seq_len output_filename\n" );
    return;
}

int 
main( int argc, char** argv )
{
    if( argc != 4 )
    {
        usage();
        exit( 1 );
    }
    
    char* genome_fname = argv[1];
    int indexed_seq_len = atoi( argv[2] );
    
    /*** Load the genome ***/
    FILE* genome_fp = fopen( genome_fname, "r");
    struct genome_data* genome;
    init_genome( &genome );
    /* load the genome into memory */
    add_chrs_from_fasta_file( genome, genome_fp );
    
    /* index the genome */
    index_genome( genome, indexed_seq_len );

    /* write the genome to file */
    write_genome_to_disk( genome, argv[3]  );
    
    /* write the index to disk */
    char buffer[500];
    sprintf( buffer, "%s.index", argv[3] );
    build_ondisk_index( genome->index, buffer  );

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
