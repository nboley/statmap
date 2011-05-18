#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "error_correction.h"

void
init_error_data( 
    struct error_data_t** data,
    int max_read_length
)
{
    *data = calloc(  sizeof(struct error_data_t), 1 );
    
    (*data)->max_read_length = max_read_length;
    (*data)->position_mismatch_cnts = calloc( sizeof(int), max_read_length  );
    
    int j;
    for( j = 0; j < max_num_qual_scores; j++ )
    {
        /* add a little bias to the scores */
        (*data)->qual_score_cnts[j] = 0;
        (*data)->qual_score_mismatch_cnts[j] = 0;
    }
    
    return;
}

void 
free_error_data( struct error_data_t* data )
{
    free( data->position_mismatch_cnts );
    free( data );
    
    return;
}

void
update_error_data(
    struct error_data_t* data,
    char* genome_seq,
    char* read,
    char* error_str,
    int length
)
{
    data->num_unique_reads += 1;
    
    int i;
    for( i = 0; i < length; i++ )
    {
        if( toupper(read[i]) != toupper(genome_seq[i])  )
        {
            data->qual_score_mismatch_cnts[ (unsigned char) error_str[i] ] += 1;
            data->position_mismatch_cnts[ i ] += 1;
        }
        
        data->qual_score_cnts[ (unsigned char) error_str[i] ] += 1;
        data->num_unique_reads += 1;
    }
    
    return;
}

void
fprintf_error_data( FILE* stream, struct error_data_t* data )
{
    fprintf( stream, "Num Unique Reads:\t%i\n", data->num_unique_reads );

    fprintf( stream, "Loc Error Rates:\n" );

    int i;
    for( i = 0; i < data->max_read_length; i++ )
    {
        fprintf( stream, "%i\t%e\n", i+1, 
                 ((float)data->position_mismatch_cnts[ i ])/data->num_unique_reads );
    }
    
    fprintf( stream, "Qual Score Error Rates:\n" );
    for( i = 0; i < max_num_qual_scores; i++ )
    {
        fprintf( stream, "%i\t%e\n", i+1, 
                 ((float)data->qual_score_mismatch_cnts[ i ])/data->qual_score_cnts[ i ] );
    }
    
    fprintf( stream, "\n" );
    
    return;
}
