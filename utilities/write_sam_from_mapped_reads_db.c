#include <stdio.h>
#include <stdlib.h>

#include "../src/mapped_read.h"
#include "../src/rawread.h"

/* make it possible to link */
int min_num_hq_bps = -1;

/* This is to silence a link error with including genome.o. It's OK because this never
   indexes the genome, so it never needs to free it */

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
    fprintf( stderr, "Usage: ./mapped_reads_to_sam genome.fa mapped_reads.db raw_reads.fastq(s)\n" );
}

int main( int argc, char** argv )
{
    if( argc != 4 && argc != 5 )
    {
        usage();
        exit(1);
    }

    /* load the genome */
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

    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    init_rawread_db( &raw_rdb );
    
    /* if the reads are single ended */
    if( argc == 4 )
    {
        add_single_end_reads_to_rawread_db(
            raw_rdb, argv[3], FASTQ 
        );
    } else if( argc == 5 )
    {
        add_paired_end_reads_to_rawread_db(
            raw_rdb, 
            argv[3],
            argv[4],
            FASTQ 
        );
    }
    
    write_mapped_reads_to_sam( raw_rdb, mpd_rdb, genome, true, false, stdout );

cleanup:
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );
}



