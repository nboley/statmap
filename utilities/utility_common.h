#ifndef UTILITY_COMMON_H
#define UTILITY_COMMON_H

/* reading directory entries */
#include <dirent.h>

#include "../src/config.h"

#include "../src/rawread.h"

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

static int
determine_next_sample_index()
{
    DIR *dp;
    struct dirent *ep;

    dp = opendir ( BOOTSTRAP_SAMPLES_ALL_PATH );
    if( NULL == dp )
    {
        perror( "FATAL       : Could not open bootstrap samples directory." );
        exit( 1 );
    }

    int sample_index = 0;
    while( 0 != ( ep = readdir (dp) ) )
    {
        /* skip the references */
        if( ep->d_name[0] == '.' )
            continue;     
        sample_index += 1;
    }

    return sample_index;
}

static void
populate_rawread_db( 
    struct rawread_db_t** raw_rdb,
    char* unpaired_fname,
    char* pair1_fname,
    char* pair2_fname )
{
    *raw_rdb = NULL;
    
    /* try to open the single end file to see if the reads are single ended */
    FILE* tmp = fopen( unpaired_fname, "r" );
    if( tmp != NULL )
    {
        fclose( tmp );
        init_rawread_db( raw_rdb );
        add_single_end_reads_to_rawread_db(
            *raw_rdb, unpaired_fname, FASTQ, UNKNOWN
        );
        
        return;
    } 

    tmp = fopen( pair1_fname, "r" );
    if( tmp != NULL )
    {
        fclose( tmp );
        init_rawread_db( raw_rdb );
        
        add_paired_end_reads_to_rawread_db(
            *raw_rdb, 
            pair1_fname,
            pair2_fname,
            FASTQ,
            UNKNOWN
        );
        
        return;
    }
}

#endif // UTILITY_COMMON_H
