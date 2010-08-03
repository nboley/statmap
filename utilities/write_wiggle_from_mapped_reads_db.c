#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../src/mapped_read.h"
#include "../src/genome.h"
#include "../src/rawread.h"
#include "../src/wiggle.h"

/* make it possible to link */
int min_num_hq_bps = -1;


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
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory fwd_output.wig bkwd_output.wig\n" );
}

int main( int argc, char** argv )
{
    if( argc != 4 )
    {
        usage();
        exit(1);
    }

    /* Initialize variables to stored mapped info */
    struct genome_data* genome = NULL;
    struct mapped_reads_db* mpd_rdb = NULL;

    /* open the output files */
    FILE* fwd_wig_fp = NULL;
    fwd_wig_fp = open_check_error( argv[2], "w" );
    FILE* bkwd_wig_fp = NULL;
    bkwd_wig_fp = open_check_error( argv[3], "w" );

    /* change to the directory where the mapped output is stored */
    int error = chdir( argv[1] );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit(1);
    }
    
    /* Load the genome */
    init_genome( &genome );
    FILE* genome_fp = open_check_error( "genome.fa", "r" );
    add_chrs_from_fasta_file( genome,  genome_fp );
    // BUG Why does it segfault when I close the file handle?
    
    /* load the mapped read db */
    char* mpd_rd_fname = "mapped_reads.db";
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );
    
    write_marginal_mapped_reads_to_stranded_wiggles(
        mpd_rdb, genome, fwd_wig_fp, bkwd_wig_fp );
    
    goto cleanup;

cleanup:
    if( NULL != genome ) 
        free_genome( genome );

    if( NULL != mpd_rdb )
        close_mapped_reads_db( mpd_rdb );

    fclose( fwd_wig_fp );
    fclose( bkwd_wig_fp );
    
    return 0;
}
