/* Copyright (c) 2009-2010, Nathan Boley */

#include "config_parsing.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
fprintf_name_or_null( FILE* arg_fp, const char* header, char* value )
{
    fprintf( arg_fp, "%s:\t", header );
    if( NULL == value )
    {
        fprintf( arg_fp, "NULL\n" );
    } else {
        fprintf( arg_fp, "%s\n", value );
    }

    return;
}

void
write_config_file_to_stream( FILE* arg_fp, struct args_t* args  )
{
    fprintf_name_or_null( 
        arg_fp, "genome_fname", args->genome_fname );
    fprintf_name_or_null( 
        arg_fp, "genome_index_fname", args->genome_index_fname );

    fprintf_name_or_null( 
        arg_fp, "unpaired_reads_fnames", args->unpaired_reads_fnames );
    fprintf_name_or_null( 
        arg_fp, "pair1_reads_fnames", args->pair1_reads_fnames );
    fprintf_name_or_null( 
        arg_fp, "pair2_reads_fnames", args->pair2_reads_fnames );
    // struct rawread_db_t* rdb;

    fprintf_name_or_null( 
        arg_fp, "unpaired_NC_reads_fnames", args->unpaired_NC_reads_fnames );
    fprintf_name_or_null( 
        arg_fp, "pair1_NC_reads_fnames", args->pair1_NC_reads_fnames );
    fprintf_name_or_null( 
        arg_fp, "pair2_NC_reads_fnames", args->pair2_NC_reads_fnames );
    // struct rawread_db_t* NC_rdb;

    fprintf_name_or_null( 
        arg_fp, "frag_len_fname", args->frag_len_fname );
    // FILE* frag_len_fp;

    fprintf_name_or_null( 
        arg_fp, "output_directory", args->output_directory );

    fprintf_name_or_null( 
        arg_fp, "sam_output_fname", args->sam_output_fname );


    fprintf_name_or_null( 
        arg_fp, "log_fname", args->log_fname );
    // FILE* log_fp;

    fprintf_name_or_null( 
        arg_fp, "log_fname", args->log_fname );

    fprintf( arg_fp, "min_match_penalty:\t%.4f\n", args->min_match_penalty );
    fprintf( arg_fp, "max_penalty_spread:\t%.4f\n", args->max_penalty_spread );
    fprintf( arg_fp, "min_num_hq_bps:\t%i\n", args->min_num_hq_bps );

    fprintf( arg_fp, "num_starting_locations:\t%i\n", args->num_starting_locations );

    fprintf( arg_fp, "indexed_seq_len:\t%i\n", args->indexed_seq_len );

    fprintf( arg_fp, "num_threads:\t%i\n", args->num_threads );

    fprintf( arg_fp, "search_type:\t%i\n", args->search_type );

    fprintf( arg_fp, "input_file_type:\t%i\n", args->input_file_type );
    fprintf( arg_fp, "assay_type:\t%i\n", args->assay_type );

    return;
}

void
fscanf_name_or_null( FILE* arg_fp, const char* header, char** value ) 
{
    // first, build the formatting strinf
    char format_string[200];
    sprintf( format_string, "%s:\t%%s\n", header);
    
    // allocate space to store the read value
    char tmp_value[500];
    
    // scan for the value
    fscanf( arg_fp, format_string, &tmp_value );
    
    // if the value is equal to NULL, then set *value= NULL
    // this indicates there was no argument in the initial
    // argument list
    if( 0 == strcmp( tmp_value, "NULL" ) )
    {
        *value = NULL;
        return;
    } 
    // if it's not NULL, then allocate space for the new string,
    // and copy the option over
    else {
        *value = calloc(1, strlen(tmp_value)+1);
        strncpy( *value, tmp_value, strlen(tmp_value) );
        return;
    }
}

void
read_config_file_fname_from_disk( char* fname, struct args_t** args  )
{
    *args = malloc( sizeof(struct args_t)  );
    FILE* arg_fp = fopen( fname, "r" );

    fscanf_name_or_null( 
        arg_fp, "genome_fname", &((*args)->genome_fname) );
    fscanf_name_or_null( 
        arg_fp, "genome_index_fname", &((*args)->genome_index_fname) );

    fscanf_name_or_null( 
        arg_fp, "unpaired_reads_fnames", &((*args)->unpaired_reads_fnames) );
    fscanf_name_or_null( 
        arg_fp, "pair1_reads_fnames", &((*args)->pair1_reads_fnames) );
    fscanf_name_or_null( 
        arg_fp, "pair2_reads_fnames", &((*args)->pair2_reads_fnames) );
    // struct rawread_db_t* rdb;

    fscanf_name_or_null( 
        arg_fp, "unpaired_NC_reads_fnames", &((*args)->unpaired_NC_reads_fnames) );
    fscanf_name_or_null( 
        arg_fp, "pair1_NC_reads_fnames", &((*args)->pair1_NC_reads_fnames) );
    fscanf_name_or_null( 
        arg_fp, "pair2_NC_reads_fnames", &((*args)->pair2_NC_reads_fnames) );
    // struct rawread_db_t* NC_rdb;

    fscanf_name_or_null( 
        arg_fp, "frag_len_fname", &((*args)->frag_len_fname) );
    // FILE* frag_len_fp;

    fscanf_name_or_null( 
        arg_fp, "output_directory", &((*args)->output_directory) );

    fscanf_name_or_null( 
        arg_fp, "sam_output_fname", &((*args)->sam_output_fname) );


    fscanf_name_or_null( 
        arg_fp, "log_fname", &((*args)->log_fname) );
    // FILE* log_fp;

    fscanf_name_or_null( 
        arg_fp, "log_fname", &((*args)->log_fname) );

    fscanf( arg_fp, "min_match_penalty:\t%f\n", 
            &((*args)->min_match_penalty) );
    fscanf( arg_fp, "max_penalty_spread:\t%f\n", 
            &((*args)->max_penalty_spread) );
    fscanf( arg_fp, "min_num_hq_bps:\t%i\n", 
            &((*args)->min_num_hq_bps) );

    fscanf( arg_fp, "num_starting_locations:\t%i\n", 
            &((*args)->num_starting_locations) );

    fscanf( arg_fp, "indexed_seq_len:\t%i\n", 
            &((*args)->indexed_seq_len) );

    fscanf( arg_fp, "num_threads:\t%i\n", 
            &((*args)->num_threads) );

    fscanf( arg_fp, "search_type:\t%d\n",
            &((*args)->search_type) );

    fscanf( arg_fp, "input_file_type:\t%d\n", 
            &( (*args)->input_file_type) );
    fscanf( arg_fp, "assay_type:\t%d\n", 
            &( (*args)->assay_type) );
    
    fclose( arg_fp  );
    
    return;
}

void
read_config_file_from_disk( struct args_t** args  )
{
    read_config_file_fname_from_disk( "config.dat", args );
}

/*
// TEST CODE. I used this on a couple runs to ensure that 
// test_config was identical to the input file. 
int main( int argc, char** argv )
{
    struct args_t* args;
    read_config_file_fname_from_disk( argv[1], &args  );
    write_config_file_to_stream( stderr, args );
    FILE* fp = fopen( "test_config.dat", "w" );
    write_config_file_to_stream( fp, args );
    fclose( fp );

}
*/
