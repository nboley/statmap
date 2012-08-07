#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include <sys/stat.h>
/* mmap() is defined in this header */
#include <sys/mman.h> 

#include "genome.h"
#include "fragment_length.h"
#include "mapped_read.h"
#include "candidate_mapping.h"
#include "find_candidate_mappings.h"
#include "pseudo_location.h"

// this is needed for the wiggle writing code
#include "iterative_mapping.h"
#include "trace.h"


/*************************************************************************
 *
 *  Mapped Reads
 * 
 *  Reads that have been joined, but unlike mapped reads proper have not
 *  had the 'extra' read inforamtion attached.
 *
 */

void
init_mapped_read( struct mapped_read_t** rd )
{
    (*rd) = malloc(sizeof(struct mapped_read_t));
    (*rd)->read_id = 0;
    (*rd)->num_mappings = 0;
    (*rd)->rdb = NULL;
    (*rd)->free_locations = true;
    (*rd)->locations = NULL;
    return;
}

void
free_mapped_read( struct mapped_read_t* rd )
{
    if( rd == NULL )
        return;
    
    if( true == rd->free_locations ) {
        free( rd->locations );
        rd->locations = NULL;
    }
    
    free( rd );
    return;
}

void
add_location_to_mapped_read( 
    struct mapped_read_t* rd, struct mapped_read_location* loc )
{
    assert( (rd->num_mappings + 1) > rd->num_mappings );
    
    /* Allocate new space for the location */
    rd->num_mappings += 1;
    
    rd->locations = realloc( 
        rd->locations, (rd->num_mappings)*sizeof(struct mapped_read_location)  );
    if( rd->locations == NULL )
    {
        fprintf(stderr, "FATAL       :  Memory allocation error ( size %zu )\n", (rd->num_mappings)*sizeof(struct mapped_read_location) );
        assert( false );
        exit( -1 );
    }

    rd->locations[rd->num_mappings - 1] = *loc; 
    
    return;
}


void
reset_read_cond_probs( struct cond_prbs_db_t* cond_prbs_db, struct mapped_read_t* rd )
{
    struct fragment_length_dist_t* fl_dist = NULL;
    if( rd->rdb != NULL )
        fl_dist = rd->rdb->fl_dist;
    
    if( 0 == rd->num_mappings )
        return;
    
    float *prbs = calloc( rd->num_mappings, sizeof(float)  );
    
    /* prevent divide by zero */
    double prb_sum = ML_PRB_MIN;
    MPD_RD_ID_T i;
    for( i = 0; i < rd->num_mappings; i++ )
    {
        struct mapped_read_location* loc = rd->locations + i;
        float cond_prob = get_seq_error_from_mapped_read_location( loc );
        if( get_flag_from_mapped_read_location( loc )&IS_PAIRED )
            cond_prob *= get_fl_prb( fl_dist, get_fl_from_mapped_read_location( loc ) );
        prbs[i] = cond_prob;
        prb_sum += cond_prob;
    }
    assert( rd->num_mappings == 0 || prb_sum > ML_PRB_MIN );

    for( i = 0; i < rd->num_mappings; i++ )
    {
        set_cond_prb( cond_prbs_db, rd->read_id, i, prbs[i]/prb_sum );
        
        assert( (prbs[i]/prb_sum < 1.001) && (prbs[i]/prb_sum) >= -0.001 );
    }
    
    free( prbs );
    
    return;
};

int 
write_mapped_read_to_file( struct mapped_read_t* read, FILE* of  )
{
    size_t num = 0;

    num = fwrite( &(read->read_id), sizeof( read->read_id ), 1, of );
    if( num != 1 )
        return -num;

    num = fwrite( &(read->num_mappings), sizeof( read->num_mappings ), 1, of );
    if( num != 1 )
        return -num;

    num = fwrite( read->locations, 
                  sizeof( struct mapped_read_location ), 
                  read->num_mappings, 
                  of );
    if( (int)num != read->num_mappings )
        return -num;
    
    return 0;
}

/* returns 1 for success, 0 for failure */
static inline int
convert_unpaired_candidate_mapping_into_mapped_read( 
    candidate_mapping* cm,
    struct mapped_read_location* loc    
)
{
    /* Ensure all of the flags are turned off */
    MRL_FLAG_TYPE flag = 0;

    /* deal with pseudo location */
    set_chr_in_mapped_read_location( loc, cm->chr );
    if( PSEUDO_LOC_CHR_INDEX == cm->chr )
    {
        flag |= FIRST_READ_IS_PSEUDO;
    }
 
    /* Add the location */
    if( BKWD == cm->rd_strnd )
        flag |= FIRST_READ_WAS_REV_COMPLEMENTED;
   
    set_start_and_stop_in_mapped_read_location (
        loc, cm->start_bp, cm->start_bp + cm->rd_len );
    
    set_seq_error_in_mapped_read_location( 
        loc, pow( 10, cm->penalty ) );
    assert( loc->seq_error >= 0.0 && loc->seq_error <= 1.0 );
    
    set_flag_in_mapped_read_location( loc, flag  );
    
    return 1;
}


/* returns 1 for success, 0 for failure */
static inline int
join_two_candidate_mappings( 
    candidate_mapping* first_read,
    candidate_mapping* second_read,
    struct mapped_read_location* loc    
)
{
    /* Ensure all of the flags are turned off */
    MRL_FLAG_TYPE flag = IS_PAIRED;
    
    /* Set the appropriate flags */
    if( first_read->start_bp < second_read->start_bp )
        flag |= FIRST_PAIR_IS_FIRST_IN_GENOME;
    
    if( FWD == first_read->rd_strnd )
        flag |= FIRST_READ_WAS_REV_COMPLEMENTED;
    
    if( PSEUDO_LOC_CHR_INDEX == first_read->chr )
        flag |= FIRST_READ_IS_PSEUDO;
    
    if( PSEUDO_LOC_CHR_INDEX == second_read->chr )
        flag |= SECOND_READ_IS_PSEUDO;
    
    /* Set the chr */
    set_chr_in_mapped_read_location( loc, first_read->chr );
    assert( first_read->chr == second_read->chr);
    
    if( first_read->start_bp < second_read->start_bp )
    {
        set_start_and_stop_in_mapped_read_location (
            loc,
            first_read->start_bp,
            second_read->start_bp + second_read->rd_len
        );
    } else {
        set_start_and_stop_in_mapped_read_location (
            loc,
            second_read->start_bp,
            first_read->start_bp + first_read->rd_len
        );
    }
    
    set_seq_error_in_mapped_read_location( 
        loc, pow( 10, first_read->penalty + second_read->penalty ) );
    
    set_flag_in_mapped_read_location( loc, flag );
    
    /* ignore reads with zero probability ( possible with FL dist ) */
    if( get_seq_error_from_mapped_read_location( loc ) > 2*ML_PRB_MIN ) {
        return 1;
    }
    
    return 0;
}

static inline double
mapped_read_from_candidate_mapping_arrays( 
    candidate_mapping* r1_array,
    int r1_array_len,
    candidate_mapping* r2_array,
    int r2_array_len,
    struct mapped_read_t* mpd_rd
)
{
    struct mapped_read_location loc;
    
    double prob_sum = 0;
    
    int i;
    for( i = 0; i < r1_array_len; i++ )
    {
        int j;
        for( j = 0; j < r2_array_len; j++ )
        {
            /* if the chrs mismatch, since these are sorted, we know that
             * there is no need to continue */
            if( r2_array[j].chr > r1_array[i].chr )
                break;
            
            if( r2_array[j].chr != r1_array[i].chr )
                continue;
            
            /* Determine which of the candidate mappings corresponds 
               with the first pair */
            candidate_mapping* first_read = NULL;
            candidate_mapping* second_read = NULL;
            if( PAIRED_END_1 == r1_array[i].rd_type ) {
                first_read = r1_array + i;
                second_read = r2_array + j;
            } else {
                assert( PAIRED_END_1 == r2_array[j].rd_type );
                first_read = r2_array + j;
                second_read = r1_array + i;
            }
            
            assert( first_read->chr == second_read->chr );
            
            int rv = join_two_candidate_mappings( first_read, second_read, &loc );
            /* if the location is valid ( non-zero probability ) */
            if( rv == 1 )
            {
                add_location_to_mapped_read( mpd_rd, &loc );
                prob_sum += get_seq_error_from_mapped_read_location( &loc );
            }
        }
    }
    
    return prob_sum;
}

double
add_paired_candidate_mappings_to_mapped_read(
    candidate_mapping* cm1,
    candidate_mapping* cm2,
    struct mapped_read_t* mpd_rd
)
{
    /* if the chrs mismatch, this match is impossible continue */
    double prob_sum = 0;
    if( cm1->chr != cm2->chr )
        return prob_sum;

    /*
     * Determine which of the candidate mappings corresponds 
     * with the first pair
     * since pair is a property of the read, the original candidate
     * mapping results still holds, even if it is a pseudo loc
     */
       
    candidate_mapping* first_read = NULL;
    candidate_mapping* second_read = NULL;

    if( PAIRED_END_1 == cm1->rd_type ) {
        first_read = cm1;
        second_read = cm2;
    } else {
        assert( PAIRED_END_1 == cm2->rd_type );
        first_read = cm2;
        second_read = cm1;
    }

    struct mapped_read_location loc;
    int rv = join_two_candidate_mappings( first_read, second_read, &loc );
    loc.position.start_pos -= cm1->subseq_offset; // ?

    /* if the location is valid ( non-zero probability ) */
    if( rv == 1 )
    {
        add_location_to_mapped_read( mpd_rd, &loc );
        prob_sum = get_seq_error_from_mapped_read_location( &loc );
    }

    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_paired_candidate_mappings( 
        struct mapped_read_t* mpd_rd,
        candidate_mappings* mappings
    )
{
    double prob_sum = 0;
    
    assert( mappings->length > 0 );
    assert( mappings->mappings[0].rd_type > SINGLE_END ); 

    int p2_start = -1;

    /* 
     * Find the start of the pair 2 mapped reads
     */
    int i;
    for( i=0; i < mappings->length; i++ )
    {
        if( mappings->mappings[i].rd_type == PAIRED_END_2 )
        {
            p2_start = i;
            break;
        }
    }
    
    /* If there were no paired end 2 reads, then we say no reads mapped */
    if( -1 == p2_start ) {
        /* make sure that the number of reads is 0 */
        assert( mpd_rd->num_mappings == 0);
        return 0;
    }

    prob_sum += mapped_read_from_candidate_mapping_arrays(
        mappings->mappings,
        p2_start,
        mappings->mappings + p2_start,
        mappings->length - p2_start,
        mpd_rd
    );

    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_unpaired_candidate_mappings( 
        struct mapped_read_t* mpd_rd,
        candidate_mappings* mappings
    )
{
    double prob_sum = 0;

    /* store local read location data */
    struct mapped_read_location loc;
    
    int i;    
    for( i = 0; i < mappings->length; i++ )
    {        
        /* we expect a read to be either paired end or not - never both */
        assert( mappings->mappings[i].rd_type == SINGLE_END );

        /* if the mapping hasnt been determined to be valid, ignore it */
        if( (mappings->mappings)[i].recheck != VALID )
            continue;

        int rv = 
            convert_unpaired_candidate_mapping_into_mapped_read( 
                mappings->mappings + i, &loc );
        
        /* if the conversion succeeded */
        if( 1 == rv )
        {
            assert( loc.seq_error >= 0.0 && loc.seq_error <= 1.0 );
            prob_sum += get_seq_error_from_mapped_read_location( &loc );
            add_location_to_mapped_read( mpd_rd, &loc );
        }
    }

    return prob_sum;
}

void
build_mapped_read_from_candidate_mappings( 
        struct mapped_read_t** mpd_rd,
        candidate_mappings* mappings, 
        MPD_RD_ID_T read_id
    )
{
    /* 
     * Building mapped reads has several components:
     * First, the normal reads
     * Second, the paired end reads
     *    Paired end read can match iff:
     *    1) They are from the same read
     *    2) They are on opposite strands
     *    3) They are on the same chromosome
     *
     *    Luckily, every passed candidate is from the same
     *    read, so (1) is established automatically. Second, 
     *    the candidates are sorted ( in cmp_candidate_mappings )
     *    by read_type, strand, chr, bp_position. So, the plan is
     *
     * 1) Get and print all of the single ended reads ( rd_type == SINGLE_END )
     * 2) Find the Start of the PAIRED_END_1 and PAIRED_END_2 reads
     * 3) Merge join them 
     * 
     */
    
    /* Initialize the packed mapped read */
    init_mapped_read( mpd_rd );
    (*mpd_rd)->read_id = read_id;

    if( mappings->length == 0 )
        return;
    
    /* store the sum of the marginal probabilities */
    double prob_sum;

    /* this assumes all reads are either paired or not */
    if( mappings->mappings[0].rd_type == SINGLE_END )
    {
        prob_sum = build_mapped_read_from_unpaired_candidate_mappings( 
            *mpd_rd, mappings );
    } else {
        prob_sum = build_mapped_read_from_paired_candidate_mappings( 
            *mpd_rd, mappings );
    }    
    
    /* If there are no proper mappings ( this can happen if, 
       for instance, only one pair maps ) free the read and return */
    if( 0 == prob_sum )
    {
        free_mapped_read( *mpd_rd );
        *mpd_rd = NULL;
        return;
    }
    
    return;
}


/*****************************************************************************
 *
 * Mapped Reads DB Code
 *
 ***************************************************************************/



static void
init_mapped_reads_db( 
    struct mapped_reads_db** rdb, char* fname, const char* mode )
{
    *rdb = malloc(sizeof(struct mapped_reads_db));
    (*rdb)->fp = fopen( fname, mode );
    if( (*rdb)->fp == NULL )
    {
        perror("FATAL       :  Could not open mapped reads file");
        assert( false );
        exit(-1);
    }
    
    /* number of mapped reads in the DB */
    (*rdb)->num_mapped_reads = 0;
    
    /* Initialize the mode to 0, this will be set by
       the mode specific init function */
    (*rdb)->mode = 0;

    /* use a mutex to eliminate the spurious helgrind warnings that we got
       when using a spinlock */
    pthread_mutexattr_t mta;
    pthread_mutexattr_init(&mta);
    (*rdb)->mutex = malloc( sizeof(pthread_mutex_t) );
    pthread_mutex_init( (*rdb)->mutex, &mta );
    
    /* mmapped data */
    (*rdb)->mmapped_data = NULL;    
    (*rdb)->mmapped_data_size = 0;

    /* index */
    (*rdb)->index = NULL;
 
    /* fl dist */
    (*rdb)->fl_dist = NULL;

    /* store the number of times that we have iterated through the read
       db in the iterative mapping code. This is very hacky, and probably 
       shouldnt be here... TODO - remove this */
    (*rdb)->num_succ_iterations = NULL;

    /* the current read location ( either the number of reads written in 'w' 
       mode, or the read we're on in 'r' mode */
    (*rdb)->current_read = 0;

    return;
}

void
open_mapped_reads_db_for_reading(
        struct mapped_reads_db** rdb,
        char* fname
    )
{
    init_mapped_reads_db( rdb, fname, "r+" );

    mmap_mapped_reads_db( *rdb );

    // read num_mapped_reads from start of file
    (*rdb)->num_mapped_reads = *((MPD_RD_ID_T*) (*rdb)->mmapped_data);

    index_mapped_reads_db( *rdb );

    (*rdb)->mode = 'r';
}

void
open_mapped_reads_db_for_writing(
        struct mapped_reads_db** rdb,
        char* fname
    )
{
    init_mapped_reads_db( rdb, fname, "w+" );

    /* write placeholder for size of mapped reads db, this 
       will be updated when we close the mapped read db*/
    MPD_RD_ID_T placeholder = 0;
    fwrite( &placeholder, sizeof(MPD_RD_ID_T), 1, (*rdb)->fp );
    
    (*rdb)->mode = 'w';
}

void
build_fl_dist_from_file( struct mapped_reads_db* rdb, FILE* fl_fp )
{
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    return;
}

void
build_fl_dist_from_filename( struct mapped_reads_db* rdb, char* filename )
{
    FILE* fl_fp = fopen( filename, "r" );
    if( fl_fp == NULL )
    {
        fprintf( stderr, "Failed to open fl_dist from filename %s\n", filename);
        exit(-1);
    }
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    fclose( fl_fp );
}

void
close_reading_specific_portions_of_mapped_reads_db( struct mapped_reads_db* rdb)
{
    munmap_mapped_reads_db( rdb );

    if( rdb->fl_dist != NULL ) {
        free_fl_dist( &(rdb->fl_dist) );
        rdb->fl_dist = NULL;
    }

    if( rdb->index != NULL ) {
        free( rdb->index );
        rdb->index = NULL;
    }

    return;
}

void
close_writing_specific_portions_of_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* update the number of reads that we have written since the 
       file was opened. */
    fseek( rdb->fp, 0, SEEK_SET );
    fwrite( &(rdb->num_mapped_reads), sizeof(MPD_RD_ID_T), 1, rdb->fp );
    
    fclose( rdb->fp );
    
    return;
}

void
close_mapped_reads_db( struct mapped_reads_db** rdb )
{
    if( NULL == *rdb )
        return;

    if( (*rdb)->mode == 'r' )
    {
        close_reading_specific_portions_of_mapped_reads_db( *rdb );
    } else if( (*rdb)->mode == 'w' ) {
        close_writing_specific_portions_of_mapped_reads_db( *rdb );
    } else {
        perror( "Unrecognized mode for open mapped read db." );
        assert( false );
        exit( -1 );
    }
    
    pthread_mutex_destroy( (*rdb)->mutex );
    free( (*rdb)->mutex );
    
    free( *rdb );
    *rdb = NULL;
    
    return;
}

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    struct mapped_read_t* rd)
{
    if ( rdb->mode == 'r' )
    {
        fprintf(stderr, "ERROR       :  Mapped Reads DB is read-only - cannot add read.\n");
        assert( false );
        exit( -1 );
    }

    int error;

    rd->rdb = rdb;
    
    pthread_mutex_lock( rdb->mutex );
    error = write_mapped_read_to_file( rd, rdb->fp );
    rdb->num_mapped_reads += 1;
    pthread_mutex_unlock( rdb->mutex );

    if( error < 0 )
    {
        fprintf(stderr, "FATAL       :  Error writing to packed mapped reads db.\n");
        assert( false );
        exit( -1 );
    }
    
    return;
}

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if( rdb->mode != 'r' )
    {
        fprintf(stderr, "FATAL       :  Can only rewind mapped reads db that is open for reading ( mode 'r' ).\n");
        assert( false );
        exit( -1 );
    }

    /* since the rdb is mmapped, we just need to reset the current read id */
    rdb->current_read = 0;
    
    return;
}

int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* rdb, 
    struct mapped_read_t** rd
)
{
    init_mapped_read( rd );
    (*rd)->rdb = rdb;

    /* Make sure the db is open for reading */
    if( rdb->mode != 'r' )
    {
        fprintf(stderr, "FATAL       :  Cannot get read from mapped reads db unless it is open for reading.\n" );
        assert( false );
        exit( -1 );
    }

    /** Get the next read **/
    pthread_mutex_lock( rdb->mutex );
    /* if we have read every read */
    if( rdb->current_read == rdb->num_mapped_reads )
    {
        pthread_mutex_unlock( rdb->mutex );
        free_mapped_read( *rd );
        *rd = NULL;
        return EOF;
    }
    
    MPD_RD_ID_T current_read_id = rdb->current_read;
    rdb->current_read += 1;
    pthread_mutex_unlock( rdb->mutex );

    assert( current_read_id < rdb->num_mapped_reads );
    assert( sizeof(char) == 1 );
    
    /* get a pointer to the current read */
    char* read_start = rdb->index[current_read_id].ptr;
    
    /* read a mapping into the struct */
    (*rd)->read_id = *((MPD_RD_ID_T*) read_start);

    read_start += sizeof(MPD_RD_ID_T);
    (*rd)->num_mappings = *((MPD_RD_ID_T*) read_start);

    read_start += sizeof(MPD_RD_ID_T);

    (*rd)->locations = (struct mapped_read_location*) read_start;
    (*rd)->free_locations = false;
        
    return 0;
}


void
reset_all_read_cond_probs( 
    struct mapped_reads_db* rdb, struct cond_prbs_db_t* cond_prbs_db )
                           
{
    rewind_mapped_reads_db( rdb );
    struct mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) ) 
    {
        reset_read_cond_probs( cond_prbs_db, r );
        free_mapped_read( r );
    }
    
    free_mapped_read( r );
}



/* use this for wiggles */
void
update_traces_from_read_densities( 
    struct mapped_reads_db* rdb,
    struct cond_prbs_db_t* cond_prbs_db,
    struct trace_t* traces
)
{    
    zero_traces( traces );

    struct mapped_read_t* r;

    while( EOF != get_next_read_from_mapped_reads_db( rdb, &r ) )     
    {
            /* Update the trace from this mapping */
        MPD_RD_ID_T j;
        double cond_prob_sum = 0;
        for( j = 0; j < r->num_mappings; j++ )
        {
            MRL_CHR_TYPE chr_index 
                = get_chr_from_mapped_read_location( r->locations + j );
            MRL_START_POS_TYPE start
                = get_start_from_mapped_read_location( r->locations + j );
            MRL_START_POS_TYPE stop
                = get_stop_from_mapped_read_location( r->locations + j );
            float cond_prob = get_cond_prb( cond_prbs_db, r->read_id, j );
            cond_prob_sum += cond_prob;
            
            assert( cond_prob >= -0.0001 );
            assert( stop >= start );            
            assert( chr_index < traces->num_chrs );
            assert( traces->chr_lengths[chr_index] >= stop );

            unsigned int k = 0;
            for( k = start; k < stop; k++ )
            {
                /* update the trace */
                traces->traces[0][chr_index][k] 
                    += (1.0/(stop-start))*cond_prob;
            }
        }
    }
    
    return;
}

void
mmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* get the file descriptor for the file we wish to mmap */
    int fdin = fileno( rdb->fp );
    
    /* make sure the entire file has been written to disk */
    fflush( rdb->fp );

    /* check that the file is not empty before trying to mmap */
    fseek(rdb->fp, 0L, SEEK_END);   // seek to end
    long fp_size = ftell(rdb->fp);  // tell() to get size
    if( fp_size == 0 )
    {
        fprintf( stderr,
                 "FATAL       :  Cannot mmap empty mapped reads db.\n" );
        exit( 1 );
    }
    rewind( rdb->fp ); // reset fp

    /* find the size of the opened file */
    struct stat buf;
    fstat(fdin, &buf);
    rdb->mmapped_data_size = buf.st_size;
    
    #ifdef MALLOC_READS_DB
    fseek( rdb->fp, 0, SEEK_SET );

    fprintf( stderr, 
             "NOTICE        : Allocating %zu bytes for the mapped reads db.", 
             buf.st_size );
    rdb->mmapped_data = malloc( buf.st_size );
    if( NULL == rdb->mmapped_data )
    {
        fprintf( stderr, 
                 "FATAL       : Failed to allocate %zu bytes for the mapped reads.\n",
                 (size_t) buf.st_size );
        exit( 1 );
    }
    
    fread( rdb->mmapped_data, buf.st_size, 1, rdb->fp );
           
    #else
    /* mmap the file */
    rdb->mmapped_data
        = mmap( NULL, rdb->mmapped_data_size,  
                PROT_READ|PROT_WRITE, 
		MAP_POPULATE|MAP_SHARED, fdin, (off_t) 0 );

    if( rdb->mmapped_data == (void*) -1 )
    {
        char* buffer;
        buffer = malloc( sizeof(char)*500 );
        sprintf(buffer, "Can not mmap the fdescriptor '%i'", fdin );
        perror( buffer );
        assert( false );
        exit( -1 );
    }
    #endif
    
    return;
}

void
munmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if( rdb->mmapped_data == NULL )
        return;

    #ifdef MALLOC_READS_DB
    free( rdb->mmapped_data );
    #else 
    int error = munmap( rdb->mmapped_data, rdb->mmapped_data_size );
    if( error != 0 )
    {
        perror( "Could not munmap mapped reads db" );
        assert( false );
        exit( -1 );
    }
    #endif

    rdb->mmapped_data = NULL;

    rdb->mmapped_data_size = 0;

    free( rdb->index );
    rdb->index = NULL;
    
    return;
}

int
cmp_mapped_reads_db_index_t(
        const struct mapped_reads_db_index_t* i1,
        const struct mapped_reads_db_index_t* i2
    )
{
    /* sort by read id */
    return i1->read_id - i2->read_id;
}

void
index_mapped_reads_db( struct mapped_reads_db* rdb )
{
    const int REALLOC_BLOCK_SIZE = 1000000;
    
    /* initialize dynamic array of mapped_reads_db_index_t */
    MPD_RD_ID_T num_allcd_reads = REALLOC_BLOCK_SIZE;
    rdb->index = malloc(sizeof(struct mapped_reads_db_index_t)*num_allcd_reads);

    /* Copy the reads data pointer (adding the offset from num_mapped_reads) */
    /* we use a char just to have a byte indexed memory block, meaning that we
       can use pointer ariuthmetic with sizeof */
    char* read_start = rdb->mmapped_data + sizeof(MPD_RD_ID_T);
    
    /* count mmapped reads to check they match the saved mapped reads count */
    MPD_RD_ID_T num_indexed_reads = 0;

    /* Loop through all of the reads */
    while( ( (size_t)read_start - (size_t)rdb->mmapped_data ) 
           < rdb->mmapped_data_size )
    {
        /* if the array is full */
        if( rdb->num_mapped_reads + 1 == num_allcd_reads )
        {
            num_allcd_reads += REALLOC_BLOCK_SIZE;
            rdb->index = realloc(rdb->index,
                    num_allcd_reads*sizeof(struct mapped_reads_db_index_t) );
        }
        assert( num_indexed_reads < num_allcd_reads );
        
        /* add the new index element */
        MPD_RD_ID_T read_id = *((MPD_RD_ID_T*) read_start);
        rdb->index[num_indexed_reads].read_id = read_id;
        rdb->index[num_indexed_reads].ptr = read_start;

        num_indexed_reads += 1;

        /* skip over the mapped read in the mmapped memory */

        /* skip the read ID */
        assert( 1 == sizeof(char) );
        read_start += sizeof(MPD_RD_ID_T);
        MPD_RD_ID_T num_mappings = *((MPD_RD_ID_T*) read_start);
        read_start += sizeof(MPD_RD_ID_T);
        /* skip the array of mapped locations */
        read_start += num_mappings*sizeof(struct mapped_read_location);
    }

    /* reclaim any wasted memory */
    rdb->index = realloc( rdb->index,
            num_indexed_reads*sizeof(struct mapped_reads_db_index_t) );

    /* sort the index by read id - this restores synchronization with the
     * reads database */
    qsort( rdb->index,
           num_indexed_reads,
           sizeof(struct mapped_reads_db_index_t),
           (int(*)(const void*, const void*))cmp_mapped_reads_db_index_t
    );

    /* make sure that we have indexed every read */
    if( num_indexed_reads != rdb->num_mapped_reads )
    {
        fprintf( stderr, 
                 "FATAL           :  The number of indexed reads (%i) is not equal to the number of reads in the mapped read db ( %i). This may indicate that the mapped read db is corrupt.", 
                 num_indexed_reads, rdb->num_mapped_reads );
    }
    
    return;
}

/* for debugging */
void
print_mapped_reads_db_index(
        struct mapped_reads_db* rdb
    )
{
    MPD_RD_ID_T i;
    for( i = 0; i < rdb->num_mapped_reads; i++ )
    {
        fprintf(stderr, "read_id: %u, ptr: %p\n",
                rdb->index[i].read_id, rdb->index[i].ptr );
    }
}

/*
 *  END Mapped DB Reads
 *
 **************************************************************************/


/*****************************************************************************
 *
 * Conditional Probs DB Code
 *
 ***************************************************************************/


void
init_cond_prbs_db_from_mpd_rdb( 
    struct cond_prbs_db_t** cond_prbs_db,
    struct mapped_reads_db* mpd_rdb
)
{
    *cond_prbs_db = malloc( sizeof( struct cond_prbs_db_t ) );
    
    /* reset the database position */
    rewind_mapped_reads_db( mpd_rdb );
    
    /* find the maximum readid */
    MPD_RD_ID_T max_rd_id = 0;
    struct mapped_read_t* mapped_rd;
    while( EOF != get_next_read_from_mapped_reads_db( mpd_rdb, &mapped_rd ) )
    {
        max_rd_id = MAX( max_rd_id, mapped_rd->read_id );
        free_mapped_read( mapped_rd );
    }
    (*cond_prbs_db)->max_rd_id = max_rd_id;
    
    free_mapped_read( mapped_rd );
    
    /* allocate space for the prb start pointers */
    (*cond_prbs_db)->cond_read_prbs = calloc( max_rd_id+1, sizeof(ML_PRB_TYPE*) );
    
    /* allocate space for the prbs */
    rewind_mapped_reads_db( mpd_rdb );
    while( EOF != get_next_read_from_mapped_reads_db( mpd_rdb, &mapped_rd ) )
    {
        (*cond_prbs_db)->cond_read_prbs[mapped_rd->read_id] 
            = calloc( mapped_rd->num_mappings, sizeof( ML_PRB_TYPE  )  );
        free_mapped_read( mapped_rd );
    }

    free_mapped_read( mapped_rd );

    return;
}

void
free_cond_prbs_db( struct cond_prbs_db_t* cond_prbs_db )
{
    MPD_RD_ID_T i;
    for( i = 0; i < cond_prbs_db->max_rd_id+1; i++ )
    {
        if( cond_prbs_db->cond_read_prbs[i] != NULL )
            free( cond_prbs_db->cond_read_prbs[i] );
    }
    
    free( cond_prbs_db->cond_read_prbs );
    free( cond_prbs_db );
}

/*
 *  END Conditional Prbs DB Code
 *
 **************************************************************************/

