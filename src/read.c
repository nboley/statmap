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

    return;
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
    if( st->ists != NULL )
        free_indexable_subtemplates(st->ists);
    
    return;
}

void
free_read( struct read* r )
{
    if( r == NULL ) return;

    free( r->name );

    int i;
    for( i = 0; i < r->num_subtemplates; i++ )
    {
        free_read_subtemplate( r->subtemplates + i );
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
    
    rst->ists = NULL;
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
                POS_SINGLE_END,
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
                POS_PAIRED_END_1,
                num_reads_in_template
            );

        add_subtemplate_to_read(
                *r,
                r2->char_seq, r2->error_str,
                r2->length,
                POS_PAIRED_END_2,
                num_reads_in_template
            );
    }
}

int
get_next_read_from_rawread_db( 
        struct rawread_db_t* rdb,
        struct read** r,
        long max_readkey
    )
{
    lock_rawread_db( rdb );
    
    if( max_readkey >= 0
        && rdb->readkey >= (readkey_t) max_readkey )
    {
        unlock_rawread_db( rdb );
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
            unlock_rawread_db( rdb );
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
            if( rv != EOF ) {
                statmap_log(LOG_WARNING, "Read pairs appear out of sync. Skipping extra pair 2 reads.");
                while( rv != EOF ) {
                    rv = populate_rawread_from_fastq_file(
                        rdb->paired_end_2_reads, &r2, SECOND );
                }
            }

            assert( rv == EOF );
            assert( rawread_db_is_empty( rdb ) );

            *r = NULL;
            unlock_rawread_db( rdb );
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

    /***** Set the prior_read_information *****/
    /* TODO Set the max lengths from config.h, for now */
    (*r)->prior.max_ref_insert_length = REFERENCE_INSERT_LENGTH_MAX;
    (*r)->prior.max_fragment_length = FRAGMENT_LENGTH_MAX;
    /* Set the assay type on the read from the rawread db */
    (*r)->prior.assay = rdb->assay;

    /* increment the read counter */
    (*r)->read_id = rdb->readkey;
    rdb->readkey += 1;

    unlock_rawread_db( rdb );
    return 0;
}

/*** Indexable subtemplates ***/

void
init_indexable_subtemplate(
        struct indexable_subtemplate** ist,
        struct read_subtemplate* rst,

        int subseq_length,
        int subseq_offset
    )
{
    *ist = malloc( sizeof( struct indexable_subtemplate ) );

    (*ist)->subseq_length = subseq_length;
    (*ist)->subseq_offset = subseq_offset;

    (*ist)->char_seq = rst->char_seq + subseq_offset;

    (*ist)->fwd_penalty_array.length = subseq_length;
    (*ist)->fwd_penalty_array.array = rst->fwd_penalty_array->array + subseq_offset;

    /* The reverse penalty array is built from the reverse complemented read
     * sequenc\e. Since subseq_offset is the distance of the offset from the 5'
     * start of the read, subseq_offset will actually be at the end of the
     * penalty array. The offset we want then is the length of the subsequence
     * plus the length of the offset subtracted from the total read length. */
    (*ist)->rev_penalty_array.length = subseq_length;
    (*ist)->rev_penalty_array.array = rst->rev_penalty_array->array
        + ( rst->length - (subseq_offset + subseq_length));
}

void
free_indexable_subtemplate(
        struct indexable_subtemplate* ist
    )
{
    if( ist == NULL ) return;

    // NOTE - we don't free char_seq because it is a pointer into memory that
    // was allocated in the parent read subtemplate
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

    free( ists->container ); // contiguous array
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

void
cache_penalty_arrays_in_read_subtemplates(
    struct read* r, struct mapping_params* mapping_params )
{
    assert( r->num_subtemplates == mapping_params->num_penalty_arrays );
    int i;
    for( i = 0; i < r->num_subtemplates; i++ ) {
        r->subtemplates[i].fwd_penalty_array = mapping_params->fwd_penalty_arrays[i];
        r->subtemplates[i].rev_penalty_array = mapping_params->rev_penalty_arrays[i];
    }
    return;
}
