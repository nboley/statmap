#include "read.h"

void
init_read( struct read** r,
           char* readname,
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
free_read_subtemplate( struct read_subtemplate* rst )
{
    /* free the char*'s */
    free( rst->char_seq );
    free( rst->error_seq );

    /* don't free the read_subtemplate itself, for now - is being
     * allocated as part of a contiguous array */
}

void
free_read( struct read* r )
{
    free( r->name );

    /* free all of the subtemplates */
    if( r>subtemplates != NULL )
    {
        int i;
        for( i = 0; i < r->num_subtemplates; i++ )
            free_read_subtemplate( &(r->subtemplates[i]) );
        free( r->subtemplates );
    }

    free( r );
}

void
add_subtemplate_to_read( struct read* r,
        char* char_seq, char* error_seq,
        int len, int pos_in_template,
        enum READ_END end
    )
{
    r->num_subtemplates += 1;

    r->subtemplates = realloc( r->subtemplates,
            sizeof(struct read_subtemplate) * r->num_subtemplates );

    /* XXX dynamically allocate char*, or set equal to passed values?
     * I want to avoid needelessly copying these values*/
    r->subtemplates[num_subtemplates - 1].char_seq = char_seq;
    r->subtemplates[num_subtemplates - 1].error_seq = error_seq;
    r->subtemplates[num_subtemplates - 1].len = len;
    r->subtemplates[num_subtemplates - 1].pos_in_template = pos_in_template;
    r->subtemplates[num_subtemplates - 1].end = end;
}

void
fprintf_read_subtemplate( FILE* fp, struct read_subtemplate* rst )
{
    fprintf(fp, "Len: %d\tPos: %d\End: %u\n",
            rst->len, rst->pos_in_template, rst->end);
    fprintf(fp, "%s\n", rst->char_seq);
    fprintf(fp, "%s\n", rst->error_seq);
}

void
fprintf_read( FILE* fp, struct read* r )
{
    fprintf(fp, "==== %s\n\n", r->name);

    /* print out subtemplates */
    int i;
    for( i = 0; i < r->num_subtemplates; i++ )
        fprintf_rawread_subtemplate( &(r->subtemplates[i]) );

    fprintf(fp, "\n\n");
}

enum bool
filter_read( struct read* r,
             struct error_model_t* error_model )
{
    /* Might pass a NULL read (r2 in the single read case) 
       from find_candidate_mappings */
    if( r == NULL )
        return false; // Do not filter nonexistent read

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
    int subt, i;

    /* loop over each of the subtemplates in read */
    for( subt = 0; subt < r->num_subtemplates; subt++ )
    {
        /* get a pointer to the subtemplate of interest */
        struct read_subtemplate* rst = &(read->subtemplates[i]);

        /* loop over each bp in the subtemplate */
        for( i = 0; i < rst->length; i++ )
        {
            /*
               compute the inverse probability of error (quality)
               NOTE when error_prb receieves identical bp's, it returns the
               inverse automatically
             */
            double error = error_prb( rst->char_seq[i], rst->char_seq[i], 
                                      rst->error_str[i], i, error_model );
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
        struct rawread* r2,
    )
{
    assert( r1 != NULL );

    /* Initialize the read with the name from the first read */
    init_read( r, r1->name );

    /*** Add the rawreads as subtemplate(s) ***/

    /* If we're working with single end raw read */
    if( r2 == NULL )
    {
        add_subtemplate_to_read(
                r1->char_seq, r1->error_str,
                r1->length,
                0,
                r1->end 
            );
    }
    /* if we're working with paired end raw reads */
    else
    {
        assert( r1 != NULL );
        assert( r2 != NULL );

        /* populate_rawread_from_fastq_file trims the / suffix, so paired
         * end reads should have the same read name */
        assert( strcmp(r1->name, r2->name) );

        add_subtemplate_to_read(
                r1->char_seq, r1->error_str,
                r1->length,
                0,
                r1->end
            );

        add_subtemplate_to_read(
                r2->char_seq, r2->error_str,
                r2->length,
                0,
                r2->end
            );
    }
}

int
get_next_read_from_rawread_db( 
    struct rawread_db_t* rdb, readkey_t* readkey,
    struct read** r,
    long max_readkey )
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

        struct rawread* r1, r2;

        /* Load first paired end rawread */
        rv = populate_rawread_from_fastq_file(
            rdb->paired_end_1_reads, &r1, FIRST );

        if( rv == EOF )
        {
            /* make sure the mate is empty as well */
            assert( rawread_db_is_empty( rdb ) );
            *r == NULL;
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
    r->assay = rdb->assay;

    /* increment the read counter */
    rdb->readkey += 1;

    pthread_spin_unlock( rdb->lock );
    return 0;
}
