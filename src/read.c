#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <pthread.h>

#include "config.h"
#include "statmap.h"
#include "rawread.h"
#include "error_correction.h"
#include "quality.h"
#include "read.h"

void
init_read(
        struct read** r,
        char* readname
    )
{
    *r = malloc( sizeof(struct read) );

    /* Allocate memory for and copy the read name */
    (*r)->name = calloc( strlen(readname) + 1, sizeof(char) );
    memcpy( (*r)->name, readname, strlen(readname) );

    (*r)->r1 = NULL;
    (*r)->r2 = NULL;
}

void
init_subtemplate(
        struct subtemplate** st,
        char* char_seq,
        char* error_str,
        int length,
        int pos_in_template,
        enum READ_END end
    )
{
    *st = malloc( sizeof(struct subtemplate) );

    // allocate memory for char_seq and error_str
    (*st)->char_seq = malloc( length*sizeof(char) );
    (*st)->error_str = malloc( length*sizeof(char) );

    // copy strings
    memcpy( (*st)->char_seq, char_seq, length );
    memcpy( (*st)->error_str, error_str, length );

    (*st)->length = length;
    (*st)->pos_in_template = pos_in_template;
    (*st)->end = end;
}

void
free_subtemplate( struct subtemplate* st )
{
    // free the strings
    free( st->char_seq );
    free( st->error_str );

    free( st );
}

void
free_read( struct read* r )
{
    if( r == NULL ) return;

    free( r->name );

    // free the subtemplates
    if( r->r1 != NULL )
        free_subtemplate( r->r1 );
    if( r->r2 != NULL )
        free_subtemplate( r->r2 );

    free( r );
}

void
fprintf_subtemplate( FILE* fp, struct subtemplate* st )
{
    fprintf(fp, "Len: %d\tPos: %d\End: %u\n",
                st->length, st->pos_in_template, st->end);
    fprintf(fp, "%s\n", st->char_seq);
    fprintf(fp, "%s\n", st->error_str);
}

void
fprintf_subtemplate_to_fastq( FILE* fp, char* name, struct subtemplate* st )
{
    fprintf( fp, "@%s\n", name );
    fprintf( fp, "%.*s\n", st->length, st->char_seq );
    fprintf( fp, "+%s\n", name );
    fprintf( fp, "%.*s\n", st->length, st->error_str );
}

void
fprintf_read( FILE* fp, struct read* r )
{
    fprintf(fp, "==== %s\n\n", r->name);

    // print subtemplates
    fprintf_subtemplate( fp, r->r1 );
    fprintf_subtemplate( fp, r->r2 );

    fprintf(fp, "\n\n");
}

enum bool
filter_read(
        struct read* r,
        struct error_model_t* error_model
    )
{
    /***************************************************************
     * check to make sure this read is mappable
     * we consider a read 'mappable' if:
     * 1) There are enough hq bps
     *
     */

    /* Make sure the global option has been set 
       ( it's initialized to -1 ); */
    assert( min_num_hq_bps >= 0 );

    int num_hq_bps = 0;

    // loop over the subtemplates, counting hq basepairs
    int i;
    struct subtemplate* subtemplates[2] = { r->r1, r->r2 };
    for( i = 0; i < 2 && subtemplates[i] != NULL; i++ )
    {
        // get a pointer to the current subtemplate
        struct subtemplate* st = subtemplates[i];

        /* loop over each bp in the subtemplate */
        int bp;
        for( bp = 0; bp < st->length; bp++ )
        {
            /*
               compute the inverse probability of error (quality)
               NOTE when error_prb receieves identical bp's, it returns the
               inverse automatically
             */
            double error = error_prb( st->char_seq[bp], st->char_seq[bp], 
                                      st->error_str[bp], i, error_model );
            double qual = pow(10, error );
            
            /* count the number of hq basepairs */
            if( qual > 0.999 )
                num_hq_bps += 1;
        }
    }

    if ( num_hq_bps < min_num_hq_bps )
        return true;
        
    return false;
}

void
build_read_from_rawreads(
        struct read** r,
        struct rawread* r1,
        struct rawread* r2
    )
{
    assert( r1 != NULL );

    // Initialize the read with the name from the first read
    init_read( r, r1->name );

    /*** Add the rawreads as subtemplate(s) ***/

    // If we're working with single end raw reads
    if( r2 == NULL )
    {
        init_subtemplate(
                &((*r)->r1),
                r1->char_seq, r1->error_str,
                r1->length,
                0,
                r1->end
            );
    }
    // if we're working with paired end raw reads
    else
    {
        assert( r1 != NULL );
        assert( r2 != NULL );

        /* populate_rawread_from_fastq_file trims the "/" suffix, so paired
         * end reads should have the same read name */
        assert( !strcmp(r1->name, r2->name) );

        init_subtemplate(
                &((*r)->r1),
                r1->char_seq, r1->error_str,
                r1->length,
                0,
                r1->end
            );

        init_subtemplate(
                &((*r)->r2),
                r2->char_seq, r2->error_str,
                r2->length,
                -1,
                r2->end
            );
    }
}

int
get_next_read_from_rawread_db( 
        struct rawread_db_t* rdb,
        readkey_t* readkey,
        struct read** r,
        long max_readkey
    )
{
    pthread_spin_lock( rdb->lock );
    
    if( max_readkey >= 0
        && rdb->readkey >= (readkey_t) max_readkey )
    {
        pthread_spin_unlock( rdb->lock );
        *r = NULL;
        return EOF;
    }

    /* 
     * Store the return value - 
     * 0 indicates success, negative failure , EOF no more reads.
     */
    int rv = -10;

    /* If the reads are single ended */
    if( rdb->single_end_reads != NULL )
    {
        struct rawread* r1;

        rv = populate_rawread_from_fastq_file(
            rdb->single_end_reads, &r1, NORMAL );

        if( rv == EOF )
        {
            *r = NULL;
            pthread_spin_unlock( rdb->lock );
            return EOF;
        }

        build_read_from_rawreads( r, r1, NULL );
    }
    /* If the reads are paired */
    else {
        assert( rdb->paired_end_1_reads != NULL );
        assert( rdb->paired_end_2_reads != NULL );

        /* Load first paired end rawread */
        struct rawread *r1, *r2;
        rv = populate_rawread_from_fastq_file(
            rdb->paired_end_1_reads, &r1, FIRST );

        if( rv == EOF )
        {
            /* make sure the mate is empty as well */
            rv = populate_rawread_from_fastq_file(
                rdb->paired_end_2_reads, &r2, SECOND );
            assert( rv == EOF );

            assert( rawread_db_is_empty( rdb ) );

            *r = NULL;
            pthread_spin_unlock( rdb->lock );
            return EOF;
        }

        /* Load second paired end rawread */
        rv = populate_rawread_from_fastq_file(
            rdb->paired_end_2_reads, &r2, SECOND );

        /* if returned an EOF, then it should have returned for r1 */
        assert( rv == 0 );

        build_read_from_rawreads( r, r1, r2 );
    }

    /* Set the assay type on the read from the rawread db */
    (*r)->assay = rdb->assay;

    /* increment the read counter */
    *readkey = rdb->readkey;
    rdb->readkey += 1;

    pthread_spin_unlock( rdb->lock );
    return 0;
}
