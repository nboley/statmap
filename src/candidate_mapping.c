/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "statmap.h"
#include "read.h"
#include "candidate_mapping.h"
#include "mapped_read.h"
#include "genome.h"
#include "diploid_map_data.h"
#include "util.h"

struct mapped_read_t;

/** 
    Some bastard code. I need this later on when I unpack the pseudo loc reads
    so I stick it here. 
***/

int
modify_mapped_read_location_for_index_probe_offset(  
    int read_location,
    const int chr,
    const enum STRAND strnd,
    const int subseq_offset,
    const int subseq_len,
    const int read_len,
    struct genome_data* genome
) 
{
    // Check for overflow error
    if( read_location < 0 ) {
        perror( "ERROR: The read locations was less than zero in modify_mapped_read_location_for_index_probe_offset. THIS SHOULD NEVER HAPPEN, PLEASE REPORT THIS BUG.\n" );
        return -1;
    }
    
    // If this is a pseudo chromosome, we need to do these checks later.
    if( chr == PSEUDO_LOC_CHR_INDEX ) {
        perror( "ERROR: Pseudo locs should NEVER be passed to modify_mapped_read_location_for_index_probe_offset. THIS SHOULD NEVER HAPPEN, PLEASE REPORT THIS BUG.\n" );
        return -1;
    }

    /* first deal with reads that map to the 5' genome */
    if( strnd == FWD )
    {
        /* if the mapping location of the probe is less than
           the length of the probe offset, then the actual 
           read is mapping before the start of the genome, which 
           is clearly impossible 
        */
        if( read_location < subseq_offset ) 
        {
            return -1;
        } 
        /* we shift the location to the beggining of the sequence, 
           rather than the subseq that we looked at in the index  */
        else {
            read_location -= subseq_offset;
        }
                
        /* if the end of the read extends past the end of the genome
           then this mapping location is impossible, so ignore it    */
        /* note that we just shifted the read start, so it's correct to
           add the full read length without substracting off the probe 
           offset. */
        if( read_location + read_len
            > (long) genome->chr_lens[chr]      )
        {
            return -1;
        }

    } else if( strnd == BKWD ) {
        /*
          This can be very confusing, so we need to draw it out:
                  
                  
          READ - 20 basepairs
          RRRR1RRRRRRRRRR2RRRR
          SUBSEQ - 12 BASEPAIRS w/ 4 BP offset
          SSSSSSSSSSSS
                  
          If the subsequence maps to the 3' genome, that means the reverse
          complement maps to the 5' genome.
                  
          GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
          2SSSSSSSSSS1
          L
          ( where L indicates the start position of the subsequence )
                  
          So the *start* of the read in the 3' genome is at position 
          L - 4 ( the subsequence offset ) + 16 ( the read length )
        */
                
        /* this moves the read start to the beginning of the read 
           <b>in the 3' genome</b>. */


        /** check not going past the end of the gneome */
                
        /* make sure that the genome is not too short, this case should
           be pretty rare but it is possible */
        if( (long) genome->chr_lens[chr] < read_len  )
        {
            return -1;
        }
                
        /* this will actually be the read end in the 5' genome,
           so we check to make sure that it won't make the read extend
           past the end of the genome */                
        if( read_location > 
            (long) genome->chr_lens[chr]
            - ( subseq_len + subseq_offset )
            ) {
            return -1;
        }
                
        read_location += ( subseq_len + subseq_offset );             
                
        /* now we subtract off the full read length, so that we have the 
           read *end* in the 5' genome. Which is what our coordinates are 
           based upon. We do it like this to prevent overflow errors. We
           first check to make sure we have enough room to subtract, and 
           then we do 
        */
        if( read_location < read_len )
        {
            return -1;
        } else {
            read_location -= read_len;
        }
    
    } else {
        perror("IMPOSSIBLE BRANCH:  WE SHOULD NEVER NOT KNOW A LOCATIONS STRAND - IGNORING IT BUT PLEASE REPORT THIS ERROR.");
        return -1;
    }
    
    return read_location;
}

/* 
 * Bastard code for diploid mapping - builds a maternal complement from a
 * paternal candidate mapping
 *
 * Modifies a paternal candidate mapping to return its maternal complement
 * Make sure to check the read location before using
 */
candidate_mapping
convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
        struct genome_data* genome,
        candidate_mapping cm
    )
{
    int paternal_chr_index = cm.chr;
    int paternal_loc = cm.start_bp;

    /* look up maternal chr_index */
    char* prefix = get_chr_prefix( genome->chr_names[paternal_chr_index] );
    int maternal_chr_index = find_diploid_chr_index(
            genome, prefix, MATERNAL
        );
    assert( maternal_chr_index > 0 ); // pseudo chr is not allowed
    free( prefix );

    /* look up associated diploid map data structure */
    int map_data_index = get_map_data_index_from_chr_index(
            genome, paternal_chr_index
        );
    assert( map_data_index >= 0 );

    /* get maternal start pos from diploid index */
    /* locations offset because diploid index is 1-indexed, but statmap is 0-indexed */
    int maternal_start = find_diploid_locations(
            &(genome->index->diploid_maps->maps[map_data_index]),
            paternal_loc + 1
        ) - 1;
    assert( maternal_start >= 0 );

    /* modify cm to be maternal complement of original paternal cm */
    cm.chr = maternal_chr_index;
    cm.start_bp = maternal_start;

    return cm;
}


/*********************************************************************************
 *
 * Candidate Mapping Code
 *
 *********************************************************************************/

void
init_candidate_mappings( candidate_mappings** mappings )
{
    *mappings = malloc( sizeof(candidate_mappings) );
    (*mappings)->mappings =
        malloc(CAND_MAPPING_RESULTS_GROWTH_FACTOR*sizeof(candidate_mapping));
    (*mappings)->length = 0;
    (*mappings)->allocated_length = CAND_MAPPING_RESULTS_GROWTH_FACTOR;
    return;
}

candidate_mapping
init_candidate_mapping_from_template(
        struct subtemplate* st,
        float max_penalty_spread
    )
{
    /****** initialize the mapped_location info that we know  ******/
    /* copy the candidate map location template */
    candidate_mapping cand_map;
    memset( &cand_map, 0, sizeof(cand_map) );

    /* Set the read length */
    cand_map.rd_len = st->length;

    /** Set the length of the subsequence. 
     * This is the length of the sequence that we go to the index for. If it
     * is not equal to read length, then we need to do a recheck.
     */
    /* TODO - allow for subsequences */        
    /*
    cand_map.subseq_len = indexed_seq_len;
    cand_map.subseq_offset = rp->subseq_offset;
    */

    /* if read length <= seq_length, then a recheck is unnecessary */
    if( max_penalty_spread > -0.1 ) {
        cand_map.recheck = RECHECK_PENALTY;
    } else {
        cand_map.recheck = VALID;
    }
    
    /* set which type of read this is */
    switch( st->end )
    {
    case 1:
        cand_map.rd_type = SINGLE_END;
        break;
    case 2:
        cand_map.rd_type = PAIRED_END_1;
        break;
    case 3:
        cand_map.rd_type = PAIRED_END_2;
        break;
    default:
        fprintf(stderr, "FATAL - unrecognized read end '%i'\n", st->end );
        exit( -1 );
    }

    return cand_map;
}

void
add_candidate_mapping( candidate_mappings* mappings,
                       candidate_mapping* mapping     )
{
    /* 
     * test to see if there is enough allocated memory in results
     * if there isn't then realloc
     */
    if( mappings->length == mappings->allocated_length )
    {
        mappings->allocated_length += CAND_MAPPING_RESULTS_GROWTH_FACTOR;
        mappings->mappings = realloc(
            mappings->mappings,
            mappings->allocated_length*sizeof(candidate_mapping)
        );
        
        if( mappings->mappings == NULL )
        {
            fprintf(stderr, "Failed realloc in add_candidate_mapping\n");
            exit(1);
        }
    }

    /* add the new results to the end of the results set */
    /* Note that this copies the mapping */
    (mappings->mappings)[mappings->length] = (*mapping);
    mappings->length++;
    
    return;
}

void
free_candidate_mappings( candidate_mappings* mappings )
{
    free( mappings->mappings );
    free( mappings );
}

void
print_candidate_mapping( candidate_mapping* mapping )
{
    printf("Recheck:      %u\n", mapping->recheck);
    printf("Chr:          %u\n", mapping->chr);
    printf("Start BP:     %u\n", mapping->start_bp);
    printf("Read_type:    %u\n", mapping->rd_type);
    printf("Read Len:     %u\n", mapping->rd_len);
    printf("Read Strand:  %u\n", mapping->rd_strnd);
    printf("Penalty:      %.2f\n", mapping->penalty);
    // printf("Subseq Off:   %u\n", mapping->subseq_offset);
    // printf("Subseq Len:   %u\n", mapping->subseq_len);
    printf("\n");
    return;
}

void
print_candidate_mappings( candidate_mappings* mappings )
{
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        print_candidate_mapping( mappings->mappings + i );
    }
}


 int
cmp_candidate_mappings( const candidate_mapping* m1, const candidate_mapping* m2 )
{

    /* 
     * Sort Order:
     *  Read Type
     *  Strand 
     *  Chromosome
     *  BP Position
     *
     * Returns: 
     * < 0 if m1 < m2 
     * = 0 if m1 = m2 
     * > 0 if m1 > m2 
     */ 

    /* first, sort by read type */
    if( m1->rd_type != m2->rd_type )
        return m1->rd_type - m2->rd_type;

    /* next, sort by strand */
    if( m1->rd_strnd != m2->rd_strnd )
        return m1->rd_strnd - m2->rd_strnd;

    /* next, sort by chromosome */
    if( m1->chr != m2->chr )
        return m1->chr - m2->chr;
    
    /* next, sort by bp position */
    if( m1->start_bp != m2->start_bp )
        return m1->start_bp - m2->start_bp;
    
    /* if they are *still* equal, then nothing else matters */
    return 0;
}

int sort_candidate_mappings( candidate_mappings* mappings )
{
    qsort( mappings->mappings, 
           mappings->length, 
           sizeof(candidate_mapping*),
           (int(*)(const void*, const void*))cmp_candidate_mappings
    );


    return 0;
}

/*
 * END Candidate Mapping
 *********************************************************************************/


/*** Generic Candidate Mappings DB Code *********************************************/

int
cmp_joining_queue_datum( const void* a, const void* b )
{
    /* a is of type joining_queue_datum** */
    joining_queue_datum** c = (joining_queue_datum**) a;
    joining_queue_datum** d = (joining_queue_datum**) b;

    /* Make NULL's compare largest */
    if( *c == NULL && *d == NULL )
        return 0;
    if( *d == NULL )
        return -1;
    if( *c == NULL )
        return 1;
    
    if( (*c)->read_id != (*d)->read_id )
        return (*c)->read_id - (*d)->read_id;
    
    return cmp_candidate_mappings( &((*c)->mapping), 
                                   &((*d)->mapping)  );

}

void
init_candidate_mappings_db( candidate_mappings_db* db,
                            char* candidate_mappings_prefix)
{
    int error;
    
    char* cmdir_name;
    if( candidate_mappings_prefix == NULL )
    {
        /* check if this branch is used - I dont think that it is */
        assert( 0 );

        char dir_name_template[100] = "cand_mappings_XXXXXX";
        cmdir_name =  mkdtemp(dir_name_template);
        if( cmdir_name == NULL )
        {
            fprintf( stderr, 
                     "FATAL       :  Could not create temp directory with template '%s'\n", 
                     dir_name_template );
            exit( -1 );
        }            
    } else {
        cmdir_name = candidate_mappings_prefix;
        int rv = mkdir( cmdir_name, S_IRWXU | S_IRWXG | S_IRWXO );
        if( rv != 0 )
        {
            fprintf( stderr,
                     "FATAL       :  Failed (errno %i) to create dir '%s'\n", 
                     rv, cmdir_name );
            exit( -1 );
        }
    }
    
    fprintf( stderr, 
             "NOTICE      :  Created Directory for Candidate Mappings: %s\n", 
             cmdir_name 
    );
    db->data_dir = malloc( sizeof(char)*(strlen( cmdir_name )+ 1) );
    memcpy( db->data_dir, cmdir_name, sizeof(char)*(strlen( cmdir_name )+ 1) );

    char tfname[1000];
    
    /* TODO - allow for buffered reads */
    db->buffer = NULL;

    db->full_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_1st_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_1st_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_2nd_fwd_rds = malloc(num_threads*sizeof(FILE*));
    db->full_2nd_bkwd_rds = malloc(num_threads*sizeof(FILE*));
    
    error = chdir( cmdir_name );
    if( error == -1 )
    {
        perror( "FATAL       : Cannot chdir into the candidate mappings directory " );
    }
    assert( error == 0 );
    
    int i = 0;
    for( i = 0; i < num_threads; i++ )
    {
        sprintf( tfname, "thread%i.full_fwd_reads.mapped",
                 i);
        db->full_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_bkwd_reads.mapped",
                 i);
        db->full_bkwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_1st_fwd_reads.mapped",
                 i);
        db->full_1st_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_1st_bkwd_reads.mapped",
                 i);
        db->full_1st_bkwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_2nd_fwd_reads.mapped",
                 i);
        db->full_2nd_fwd_rds[i] = open_check_error( tfname, "w+" );

        sprintf( tfname, "thread%i.full_2nd_bkwd_reads.mapped",
                 i);
        db->full_2nd_bkwd_rds[i] = open_check_error( tfname, "w+" );
    }

    /* TODO - get curr dir and change back to that */
    error = chdir( "../" );
    assert( error == 0 );
 
    return;
}

void
close_candidate_mappings_db( candidate_mappings_db* db )
{  
    int i;
    for( i = 0; i < num_threads; i++ )
    {
        fclose( db->full_fwd_rds[i] );
        fclose( db->full_bkwd_rds[i] );
        fclose( db->full_1st_fwd_rds[i] );
        fclose( db->full_1st_bkwd_rds[i] );
        fclose( db->full_2nd_fwd_rds[i] );
        fclose( db->full_2nd_bkwd_rds[i] );
    }

    free( db->full_fwd_rds );
    free( db->full_bkwd_rds );
    free( db->full_1st_fwd_rds );
    free( db->full_1st_bkwd_rds );
    free( db->full_2nd_fwd_rds );
    free( db->full_2nd_bkwd_rds );
    
    free( db->data_dir );
    
    return;
}


 FILE*
determine_fp_from_candidate_mapping( 
    candidate_mapping* mapping, 
    candidate_mappings_db* db,
    int thread_num )
{
    if( mapping->rd_type == SINGLE_END )
    {
        if( mapping->rd_strnd == FWD ) {
            return db->full_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_bkwd_rds[thread_num];
        }
    } else 
    if( mapping->rd_type == PAIRED_END_1 )
    {
        if( mapping->rd_strnd == FWD ) {
            return db->full_1st_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_1st_bkwd_rds[thread_num];
        }
    } else 
    {   assert( mapping->rd_type == PAIRED_END_2 );
            
        if( mapping->rd_strnd == FWD ) {
            return db->full_2nd_fwd_rds[thread_num];
        } else {
            assert( mapping->rd_strnd == BKWD );
            return db->full_2nd_bkwd_rds[thread_num];
        }
    }
}

void
add_candidate_mappings_to_db ( 
    candidate_mappings_db* db, 
    candidate_mappings* mappings,
    long read_id,
    int thread_num )
{
    int error;

    if( mappings->length == 0  )
        return;
    
    /* TODO - add proper DB support */
    /* Determine the correct file to write this read to  */        
    /* BUG - FIXME */
    FILE* fp =  
        determine_fp_from_candidate_mapping( 
            mappings->mappings, db, thread_num );
        
    /* write the key to the file */
    /* print this as a string so that we can use unix sort */
    error = fprintf( fp, "%lu", read_id );
    /* Add the tab so that we can use unix sort */
    error = fputc( '\t', fp );
    
    int i;
    for( i = 0; i < mappings->length; i++ )
    {
        /* Move to the correct candidate mapping */
        candidate_mapping* mapping = mappings->mappings + i;

        error = putc( '\0', fp );
        error = fwrite( mapping, sizeof(candidate_mapping), 1, fp );
    }
    error = putc( '\n', fp );

    return;
}

void
close_candidate_mappings_cursor( 
    candidate_mappings_db_cursor* cursor )
{
    int i;
    for( i = 0; i < cursor->queue_len; i++ )
    {
        /* free the queue_datum ) */
        if( (cursor->mappings_queue)[i] != NULL )
        {
            free( (cursor->mappings_queue)[i] );
        }
    }
    
    free( cursor->mappings_queue );

    /* free the cursor memory */
    free( cursor );
}

static  int
get_readname_from_stream( long* read_id, 
                          FILE* stream )
{
    int error;
    error = fscanf( stream, "%lu", read_id );
    if( 1 != error )
    {
        assert( feof( stream ) );
        return EOF;
    }
    
    assert( error == 1 );
    /* Get the tab */
    fgetc( stream );
    return 0;    
}

void
open_candidate_mappings_cursor (
    candidate_mappings_db* db,
    candidate_mappings_db_cursor** cursor)
{
    int i;
    int error;
    
    *cursor = malloc( sizeof( candidate_mappings_db_cursor ) );
    (*cursor)->db = db;
    (*cursor)->queue_len = 6*num_threads;
    (*cursor)->mappings_queue 
        = malloc(sizeof(joining_queue_datum*)*(*cursor)->queue_len);
    (*cursor)->curr_read_id = 0;
    
    FILE** fps = malloc(sizeof(FILE*)*6*num_threads);
    for( i = 0; i < num_threads; i++ )
    {
        fps[6*i + 0] = db->full_fwd_rds[i];
        fps[6*i + 1] = db->full_bkwd_rds[i];
        fps[6*i + 2] = db->full_1st_fwd_rds[i];
        fps[6*i + 3] = db->full_1st_bkwd_rds[i];
        fps[6*i + 4] = db->full_2nd_fwd_rds[i];
        fps[6*i + 5] = db->full_2nd_bkwd_rds[i];
    };

    /* initialize the queue */    
    for( i = 0; i < (*cursor)->queue_len; i++ )
    {
        /* Move the fp's to the beginning of the stream */ 
        rewind( fps[i] );

        /* initialize space for the queue */
        ((*cursor)->mappings_queue)[i] 
            = malloc( sizeof( joining_queue_datum ) );

        /* add the file to the stream */
        ((*cursor)->mappings_queue)[i]->stream = fps[i];
        
        /* Get the readnames from the stream while the reads 
           dont have mappings. That is, we keep grabbing 
           readnames and testing for a new line, which indicates
           that there are no new mappings.
         */
        do { 
            error = get_readname_from_stream( 
                &(((*cursor)->mappings_queue)[i]->read_id),
                fps[i]
            );            
            error = getc( fps[i] );
        /* I dont need to test for EOF here because EOF != \n, 
           so if I get EOF it will break anyways */
        } while( error == '\n' );
        
        /* if the file is at eof, set the queue to NULL */
        if( EOF == error )
        {
            free( ((*cursor)->mappings_queue)[i] );
            ((*cursor)->mappings_queue)[i] = NULL;
          
        } else {        
            /* read in the mapped location and add it into the queue */
            error = fread( &((((*cursor)->mappings_queue)[i])->mapping), 
                    sizeof(candidate_mapping), 1, fps[i] );
            
        }
    }

    /* sort the queue */
    qsort( (*cursor)->mappings_queue, 
           (*cursor)->queue_len, 
           sizeof( joining_queue_datum* ), 
           cmp_joining_queue_datum );

    free( fps );
    
    return;
}

int
get_next_candidate_mapping_from_cursor( 
    candidate_mappings_db_cursor* cursor, 
    candidate_mappings** mappings,
    long* read_id )
{
    /* initialize a place to store the reads */
    init_candidate_mappings( mappings );
    
    /* sort the queue */
    qsort( cursor->mappings_queue, 
           cursor->queue_len, 
           sizeof( joining_queue_datum* ), 
           cmp_joining_queue_datum );
    
    /* check to see if we are done */
    if( (cursor->mappings_queue)[0] == NULL )
    {
        free_candidate_mappings( *mappings );
        return CURSOR_EMPTY;
    }
    
    /* Set the read id */
    *read_id = cursor->curr_read_id;
    /* increment the read id, for the next mapping */
    cursor->curr_read_id += 1;
    
    /* now, keep adding data until the key changes */
    while( (cursor->mappings_queue)[0] != NULL 
           && *read_id == (cursor->mappings_queue)[0]->read_id )
    {
        /* Add the current mapped location into mapped locations */
        add_candidate_mapping( *mappings, 
                               &(((cursor->mappings_queue)[0])->mapping) );
        
        /* get the next char in the top file to check if this is a new read */
        /* If we are, this should be a new line */
        char rv = getc( ((cursor->mappings_queue)[0])->stream );
        
        /* If we have moved to the next read */
        if( rv == '\n' )
        {
            /* loop until we have a read with locations, or an EOF */
            /* This is to ensure that an active read is always in the queue */
            while( rv == '\n' )
            {
                rv = get_readname_from_stream( 
                    &((cursor->mappings_queue)[0]->read_id),
                    (cursor->mappings_queue)[0]->stream
                );
                
                rv = getc( ((cursor->mappings_queue)[0])->stream );
            }
        }
                
        if( rv == EOF )
        {
            assert( feof( ((cursor->mappings_queue)[0])->stream ) );
            free( (cursor->mappings_queue)[0] );
            (cursor->mappings_queue)[0] = NULL;
        } else {
            /* read in the next mapping */
            assert( rv == '\0' );
            
            rv = fread( &((cursor->mappings_queue)[0]->mapping),
                        sizeof(candidate_mapping), 1,
                        ((cursor->mappings_queue)[0])->stream
            );
        }      
        
        /* finally, resort the queue with the new item */
        qsort( cursor->mappings_queue, 
               cursor->queue_len, 
               sizeof( joining_queue_datum* ), 
               cmp_joining_queue_datum );
        
    }

    return 0;
}

void
join_all_candidate_mappings( candidate_mappings_db* cand_mappings_db,
                             struct mapped_reads_db* mpd_rds_db,
                             struct genome_data* genome )
{
    int error;

    /* Join all candidate mappings */

    /* get the cursor to iterate through the candidate mappings */    
    candidate_mappings_db_cursor* candidate_mappings_cursor;
    open_candidate_mappings_cursor(
        cand_mappings_db, &candidate_mappings_cursor );
    
    /* store the read key of the current read */
    long read_key;
    
    candidate_mappings* mappings;
    struct mapped_read_t* mpd_rd;

    /* Get the first read */
    error = get_next_candidate_mapping_from_cursor( 
        candidate_mappings_cursor, 
        &mappings,
        &read_key 
    );
    
    while( CURSOR_EMPTY != error ) 
    {
        if( read_key > 0 && read_key%4000000 == 0 )
            fprintf( stderr, "NOTICE       : Joined %li reads\n", read_key );
        
        build_mapped_read_from_candidate_mappings( 
            genome, mappings, &mpd_rd, read_key );

        /* if we were actual able to build a mapped read with > 0 prb */
        if( NULL != mpd_rd ) {
            add_read_to_mapped_reads_db( mpd_rds_db, mpd_rd );
            free_mapped_read( mpd_rd );
        }
        
        free_candidate_mappings( mappings );

        /* Get the reads */
        error = get_next_candidate_mapping_from_cursor( 
            candidate_mappings_cursor, 
            &mappings,
            &read_key 
        );
    }
    
    goto cleanup;
    
cleanup:
    /* close the cursor */
    close_candidate_mappings_cursor( candidate_mappings_cursor );
        
    return;
}

