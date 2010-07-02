#include <stdio.h>
#include <stdlib.h>

#include "../src/mapped_read.h"
#include "../src/rawread.h"
#include "../src/wiggle.h"

/* make it possible to link */
int min_num_hq_bps = -1;

/* This is to eliminate a link error with including genome.o. It's OK because this never
   indexes the genome, so it never needs to free it. But, its pretty hacky. I probably
   need to do some refactoring...  */

void free_tree( void )
{ return; }

static FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s\n'", fname );
        exit( -1 );
    }
    return tmp;
}

void usage()
{
    fprintf( stderr, "Usage: ./mapped_reads_to_sam genome.fa mapped_reads.db fwd_output.wig bkwd_output.wig\n" );
}

int main( int argc, char** argv )
{
    if( argc != 5 )
    {
        usage();
        exit(1);
    }

    /* Load the genome */
    struct genome_data* genome;
    init_genome( &genome );
    FILE* genome_fp = open_check_error( argv[1], "r" );
    add_chrs_from_fasta_file( genome,  genome_fp );
    // BUG Why does it segfault when I close the file handle?
    
    /* load the mapped read db */
    char* mpd_rd_fname = argv[2];
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );
    
    FILE* fwd_wig_fp = open_check_error( argv[3], "w" );
    FILE* bkwd_wig_fp = open_check_error( argv[4], "w" );
    write_marginal_mapped_reads_to_stranded_wiggles(
        mpd_rdb, genome, fwd_wig_fp, bkwd_wig_fp );
    fclose( fwd_wig_fp );
    fclose( bkwd_wig_fp );

cleanup:
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
}
