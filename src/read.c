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

    (*r)->subtemplates = NULL;
    (*r)->num_subtemplates = 0;
}

void
free_read_subtemplate( struct read_subtemplate* st )
{
    if( st == NULL ) return;

    // free the dynamically allocated strings
    free( st->char_seq );
    free( st->error_str );

    /* Don't free the read subtemplate itself - it is allocated as part of a
     * contiguous array. */
}

void
free_read( struct read* r )
{
    if( r == NULL ) return;

    free( r->name );

    // free the read subtemplates
    int st_index;
    for( st_index = 0; st_index < r->num_subtemplates; st_index++ )
    {
        free_read_subtemplate( &(r->subtemplates[st_index]) );
    }
    free( r->subtemplates );

    free( r );
}

void
add_subtemplate_to_read(
        struct read* r,
        char* char_seq,
        char* error_str,
        int length,
        int pos_in_template,
        int num_reads_in_template
    )
{
    // Reallocate subtemplates array
    r->num_subtemplates += 1;
    r->subtemplates = realloc( r->subtemplates,
            sizeof(struct read_subtemplate) * r->num_subtemplates );

    // reference to new subtemplate
    struct read_subtemplate* rst = &(r->subtemplates[r->num_subtemplates-1]);

    // Allocate memory for strings
    rst->char_seq = malloc(sizeof(char)*length );
    rst->error_str = malloc(sizeof(char)*length );

    // Copy strings
    memcpy( rst->char_seq, char_seq, length );
    memcpy( rst->error_str, error_str, length );

    rst->length = length;

    /* initialize pos_in_template struct */
    rst->pos_in_template.pos = pos_in_template;
    rst->pos_in_template.number_of_reads_in_template = num_reads_in_template;
    rst->pos_in_template.is_full_fragment = false;
}

void
fprintf_read_subtemplate( FILE* fp, struct read_subtemplate* st )
{
    fprintf(fp, "Len: %d\tPos: %d\n",
                st->length, st->pos_in_template.pos );
    fprintf(fp, "%s\n", st->char_seq);
    fprintf(fp, "%s\n", st->error_str);
}

void
fprintf_read_subtemplate_to_fastq( FILE* fp, char* name, struct read_subtemplate* st )
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

    // print read subtemplates
    int st_index;
    for( st_index = 0; st_index < r->num_subtemplates; st_index++ )
    {
        fprintf_read_subtemplate( fp, &(r->subtemplates[st_index]) );
    }

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
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        // get a pointer to the current subtemplate
        struct read_subtemplate* rst = &(r->subtemplates[i]);

        /* loop over each bp in the subtemplate */
        int bp;
        for( bp = 0; bp < rst->length; bp++ )
        {
            /*
               compute the inverse probability of error (quality)
               NOTE when error_prb receieves identical bp's, it returns the
               inverse automatically
             */
            double error = error_prb( rst->char_seq[bp], rst->char_seq[bp], 
                                      rst->error_str[bp], i, error_model );
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
    int num_reads_in_template;

    // If we're working with single end raw reads
    if( r2 == NULL )
    {
        num_reads_in_template = 1;

        add_subtemplate_to_read(
                *r,
                r1->char_seq, r1->error_str,
                r1->length,
                1,
                num_reads_in_template
            );
    }
    // if we're working with paired end raw reads
    else
    {
        assert( r1 != NULL );
        assert( r2 != NULL );

        num_reads_in_template = 2;

        /* populate_rawread_from_fastq_file trims the "/" suffix, so paired
         * end reads should have the same read name */
        assert( !strcmp(r1->name, r2->name) );

        add_subtemplate_to_read(
                *r,
                r1->char_seq, r1->error_str,
                r1->length,
                1,
                num_reads_in_template
            );

        add_subtemplate_to_read(
                *r,
                r2->char_seq, r2->error_str,
                r2->length,
                -1,
                num_reads_in_template
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

        // free the rawread
        free_rawread( r1 );
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

        // free the rawreads
        free_rawread( r1 );
        free_rawread( r2 );
    }

    /* Set the assay type on the read from the rawread db */
    (*r)->assay = rdb->assay;

    /* increment the read counter */
    *readkey = rdb->readkey;
    rdb->readkey += 1;

    pthread_spin_unlock( rdb->lock );
    return 0;
}

/*** Indexable subtemplates ***/

void
init_indexable_subtemplate(
        struct indexable_subtemplate** ist
    )
{
    *ist = malloc( sizeof( struct indexable_subtemplate ) );

    // TODO - consider setting some of the below values as part of this
    // initialization function

    (*ist)->seq_length = 0;
    (*ist)->subseq_offset = 0;

    (*ist)->char_seq = NULL;
    (*ist)->error_str = NULL;

    (*ist)->origin = NULL;
}

void
free_indexable_subtemplate(
        struct indexable_subtemplate* ist
    )
{
    if( ist == NULL ) return;

    // NOTE - we don't free char_seq or error_str because they are pointers
    // into memory that was allocated as part of the parent read subtemplate
    free( ist );
}

void
init_indexable_subtemplates(
        struct indexable_subtemplates** ists
    )
{
    *ists = malloc( sizeof( struct indexable_subtemplates ) );

    (*ists)->container = NULL;
    (*ists)->length = 0;
}

void
free_indexable_subtemplates(
        struct indexable_subtemplates* ists
    )
{
    if( ists == NULL ) return;

    // free the indexable subtemplates
    int i;
    for( i = 0; i < ists->length; i++ )
    {
        free_indexable_subtemplate( &(ists->container[i]) );
    }

    free( ists );
}

// TODO - for now. May need to write a function build & add indexable
// subtemplates.
void
add_indexable_subtemplate_to_indexable_subtemplates(
        struct indexable_subtemplate* ist,
        struct indexable_subtemplates* ists
    )
{
    ists->length += 1;
    ists->container = realloc( ists->container,
            sizeof( struct indexable_subtemplate ) * ists->length );

    ists->container[ists->length-1] = *ist;
}