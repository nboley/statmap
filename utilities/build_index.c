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

int 
main( int argc, char** argv )
{
    assert( argc == 3 );
    
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
    
    /* write the index to disk */
    build_ondisk_index( genome->index, "ODIndex.bin"  );

#ifdef EXPLICITLY_FREE_INDEX
    /* free the old index */
    free_node_and_children( (static_node*) genome->index->index, 's' );
#endif
    
    free_genome( genome );

    return 0;
}
