#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "../src/mapped_read.h"
#include "../src/genome.h"
#include "../src/sam.h"
#include "../src/rawread.h"

/* make it possible to link */
int min_num_hq_bps = -1;

/* This is to silence a link error with including genome.o. It's OK because this never
   indexes the genome, so it never needs to free it */
void free_tree( void )
{ return; }
\
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
    fprintf( stderr, "Usage: ./mapped_reads_to_sam output_directory\n" );
}

int main( int argc, char** argv )
{
    if( argc != 2 )
    {
        usage();
        exit(1);
    }

    /* change to the output directory */
    int error = chdir( argv[1] );
    if( -1 == error )
    {
        perror( "FATAL       :  Cannot move into output directory ");
        exit( 1 );
    }


    /* Load the genome */
    struct genome_data* genome;
    init_genome( &genome );
    FILE* genome_fp = open_check_error( "genome.fa", "r" );
    add_chrs_from_fasta_file( genome,  genome_fp );
    // BUG Why does it segfault when I close the file handle?
    
    /* load the mapped read db */
    char* mpd_rd_fname = "mapped_reads.db";
    struct mapped_reads_db* mpd_rdb;
    open_mapped_reads_db( &mpd_rdb, mpd_rd_fname );

    /* load the rawreads */
    struct rawread_db_t* raw_rdb;
    init_rawread_db( &raw_rdb );
    
    /* try to open the single end file to see if the reads are single ended */
    FILE* tmp = fopen( "reads.unpaired", "r" );
    if( tmp != NULL )
    {
        fclose( tmp );

        add_single_end_reads_to_rawread_db(
            raw_rdb, "reads.unpaired", FASTQ 
        );
    } 
    /* else, assume the reads are paired */
    else {
        add_paired_end_reads_to_rawread_db(
            raw_rdb, 
            "reads.pair1",
            "reads.pair2",
            FASTQ 
        );
    }
    
    write_mapped_reads_to_sam( raw_rdb, mpd_rdb, genome, true, false, stdout );
    
    goto cleanup;
    
cleanup:
    free_genome( genome );
    close_mapped_reads_db( mpd_rdb );
    close_rawread_db( raw_rdb );

    return 0;
}



