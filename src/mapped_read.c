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
    long i;
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
    long num = 0;

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
    if( num != read->num_mappings )
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
add_pseudo_loc_to_mapped_read(
    struct genome_data* genome,
    candidate_mapping* cm,
    struct mapped_read_t* mpd_rd
)
{
    /* store local read location data */
    struct mapped_read_location loc;

    int rv = 
        convert_unpaired_candidate_mapping_into_mapped_read( 
            cm, &loc );

    /* if this is a pseudo location, we need to make the offset correction */
    int read_location = loc.position.start_pos;
    read_location = modify_mapped_read_location_for_index_probe_offset(  
        read_location,
        loc.chr,
        cm->rd_strnd,
        cm->subseq_offset,
        // subseq len? prbly a bug...
        cm->rd_len - cm->subseq_offset,
        cm->rd_len,
        genome
    );

    // If this location is impossible for some reason ( ie, less than 0 )
    // then skip adding it
    if( read_location < 0 ) {
        return 0;
    } else {
        loc.position.start_pos = read_location;
    }

    double seq_error = 0;
    /* if the conversion succeeded */
    if( 1 == rv )
    {
        assert( loc.seq_error >= 0.0 && loc.seq_error <= 1.0 );
        seq_error = get_seq_error_from_mapped_read_location( &loc );
        add_location_to_mapped_read( mpd_rd, &loc );
    }

    return seq_error;
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


/* CMA = Candidat eMapping Array */
static inline double
mapped_read_from_pseudo_CMA_and_pseudo_CMA( 
    candidate_mapping* pseudo_array_1,
    int pseudo_array_1_len,
    candidate_mapping* pseudo_array_2,
    int pseudo_array_2_len,
    struct genome_data* genome,
    struct pseudo_locations_t* ps_locs,
    struct mapped_read_t* mpd_rd
)
{
    double prob_sum = 0;
    
    int i, j, k, l;
    /* loop through each pseudo read */
    for( i = 0; i < pseudo_array_1_len; i++ )
    {
        int ps_loc_1_index = pseudo_array_1[i].start_bp;
        struct pseudo_location_t* ps_loc_1 = ps_locs->locs + ps_loc_1_index;

        /* loop through each real read */
        for( j = 0; j < pseudo_array_2_len; j++ )
        {
            int ps_loc_2_index = pseudo_array_2[j].start_bp;
            struct pseudo_location_t* ps_loc_2 = ps_locs->locs + ps_loc_2_index;

            /* loop through each location in the pseudo reads */
            for( k = 0; k < ps_loc_1->num; k++ )
            {
                GENOME_LOC_TYPE* gen_locs_1 = ps_loc_1->locs;
                
                pseudo_array_1[i].chr = gen_locs_1[k].chr;
                pseudo_array_1[i].start_bp = gen_locs_1[k].loc;
                
                /* loop through each location in the pseudo reads */
                for( l = 0; l < ps_loc_2->num; l++ )
                {
                    GENOME_LOC_TYPE* gen_locs_2 = ps_loc_2->locs;
                    
                    pseudo_array_2[j].chr = gen_locs_2[l].chr;
                    pseudo_array_2[j].start_bp = gen_locs_2[l].loc;

                    /* attempt to add both locs without modification */
                    prob_sum += add_paired_candidate_mappings_to_mapped_read(
                        &pseudo_array_1[i], &pseudo_array_2[j], mpd_rd
                    );

                    /* if the first ps_loc's current loc has both bit flags set,
                     * expand its maternal complement and attempt to add it to the mapped read
                     * with the original second loc */
                    if( gen_locs_1[k].is_paternal && gen_locs_1[k].is_maternal )
                    {
                        candidate_mapping maternal = 
                            convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                                genome, pseudo_array_1[i]
                            );

                        prob_sum += add_paired_candidate_mappings_to_mapped_read(
                            &maternal, &pseudo_array_2[j], mpd_rd
                        );
                    }

                    /* if the second ps_loc's current loc has both bit flags set,
                     * expand its maternal complement and add it to the mapped read
                     * with the original first loc */
                    if( gen_locs_2[l].is_paternal && gen_locs_2[l].is_maternal )
                    {
                        candidate_mapping maternal = 
                            convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                                genome, pseudo_array_2[j]
                            );

                        prob_sum += add_paired_candidate_mappings_to_mapped_read(
                            &maternal, &pseudo_array_1[i], mpd_rd
                        );
                    }

                    /* if both ps_loc's current locs have both bit flags set,
                     * expand both of their maternal complements and add them together */
                    if( gen_locs_1[k].is_paternal && gen_locs_1[k].is_maternal
                            &&
                        gen_locs_2[l].is_paternal && gen_locs_2[l].is_maternal )
                    {
                        candidate_mapping maternal_1 = 
                            convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                                genome, pseudo_array_1[i]
                            );
                        candidate_mapping maternal_2 = 
                            convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                                genome, pseudo_array_2[j]
                            );
                        prob_sum += add_paired_candidate_mappings_to_mapped_read(
                            &maternal_1, &maternal_2, mpd_rd
                        );
                    }

                }
            }
        }
    }
    
    return prob_sum;
}


/* CMA = Candidat eMapping Array */
static inline double
mapped_read_from_CMA_and_pseudo_CMA( 
    candidate_mapping* pseudo_array,
    int pseudo_array_len,
    candidate_mapping* r2_array,
    int r2_array_len,
    struct genome_data* genome,
    struct pseudo_locations_t* ps_locs,
    struct mapped_read_t* mpd_rd
)
{
    double prob_sum = 0;
    
    int i, j, k;
    /* loop through each pseudo read */
    for( i = 0; i < pseudo_array_len; i++ )
    {
        int ps_loc_index = pseudo_array[i].start_bp;
        struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;
            
        /* loop through each real read */
        for( j = 0; j < r2_array_len; j++ )
        {
            /* loop through each location in the pseudo reads */
            for( k = 0; k < ps_loc->num; k++ )
            {
                GENOME_LOC_TYPE* gen_locs = ps_loc->locs;
                
                pseudo_array[i].chr = gen_locs[k].chr;
                pseudo_array[i].start_bp = gen_locs[k].loc;

                /* add both mappings as they are */
                prob_sum += add_paired_candidate_mappings_to_mapped_read(
                    &pseudo_array[i], &r2_array[j], mpd_rd
                );

                /* if the pseudo location is shared sequence, add paired mappings
                 * with the maternal complement of the pseudo loc */
                if( gen_locs[k].is_paternal && gen_locs[k].is_maternal )
                {
                    candidate_mapping maternal =
                        convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                            genome, pseudo_array[i]
                        );

                    prob_sum += add_paired_candidate_mappings_to_mapped_read(
                        &maternal, &r2_array[j], mpd_rd
                    );
                }
            }
        }
    }
    
    return prob_sum;
}


/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_paired_candidate_mappings( 
    struct genome_data* genome,
    candidate_mappings* mappings,
    /* assume this has already been initialized */
    struct mapped_read_t* mpd_rd )
{
    double prob_sum = 0;
    
    int p1_start=-1, p1_stop=-1, p2_start=-1;

    assert( mappings->length > 0 );
    assert( mappings->mappings[0].rd_type > SINGLE_END ); 

    /* 
       Find relevant indexes, 
       1) the start of pseudo indexes is always 0
       2) p1_start - the start of non pseudo pairs, and the end of p1 pseudo
       3) p1_stop - the end of the first read pairs, and start of p2 pseudo
       4) p2_start - the start of pair 2, and end of pair 2 pseudo
       5) the end of pair two is mappings->length
    */
    int i;
    for( i=0; i < mappings->length; i++ )
    {
        if( p1_start == -1 && mappings->mappings[i].chr > 0 )
            p1_start = i;

        if( p1_stop == -1 && mappings->mappings[i].rd_type == PAIRED_END_2 )
            p1_stop = i;

        if( p2_start == -1
            && mappings->mappings[i].rd_type == PAIRED_END_2
            && mappings->mappings[i].chr > 0 )
            p2_start = i;
    }
    
    /* If we are done, return ( note that we dont have any 
       mapped paired end 2's so no paired ends mapped ) */
    if( -1 == p1_stop ) {
        /* make sure that the number of reads is 0 */
        assert( mpd_rd->num_mappings == 0);
        return 0;
    }

    /* if either of these weren't set in the previous loop, there were no
       non-pseudo reads for that set of read pairs. set start=stop to represent this. */
    if( p1_start == -1 )
        p1_start = p1_stop;
    if( p2_start == -1 )
        p2_start = mappings->length;

    /* join the pseudo location pairs */
    prob_sum += mapped_read_from_pseudo_CMA_and_pseudo_CMA(
        mappings->mappings,
        p1_start,
        mappings->mappings + p1_stop,
        p2_start - p1_stop,
        genome,
        genome->index->ps_locs,
        mpd_rd
    );
    
    /* join the first pseudo location locations */
    prob_sum += mapped_read_from_CMA_and_pseudo_CMA(
        mappings->mappings,
        p1_start,
        mappings->mappings + p2_start,
        mappings->length - p2_start,
        genome,
        genome->index->ps_locs,
        mpd_rd
    );
    
    /* join the non-pseudo location locations */
    prob_sum += mapped_read_from_candidate_mapping_arrays(
        mappings->mappings + p1_start,
        p1_stop - p1_start,
        mappings->mappings + p2_start,
        mappings->length - p2_start,
        mpd_rd
    );

    /* join the second pseudo location locations */
    prob_sum += mapped_read_from_CMA_and_pseudo_CMA(
        mappings->mappings + p1_stop,
        p2_start - p1_stop,
        mappings->mappings + p1_start,
        p1_stop - p1_start,
        genome,
        genome->index->ps_locs,
        mpd_rd
    );
    
    return prob_sum;
}

/* returns the sum of sequencing error probabilities - used for renormalization */
static inline double
build_mapped_read_from_unpaired_candidate_mappings( 
    struct genome_data* genome,
    candidate_mappings* mappings,
    /* assume this has already been initialized */
    struct mapped_read_t* mpd_rd )
{
    double prob_sum = 0;

    struct pseudo_locations_t* ps_locs = genome->index->ps_locs;

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

        /* deal with pseudo locations */
        if( EXPAND_UNPAIRED_PSEUDO_LOCATIONS
            && PSEUDO_LOC_CHR_INDEX == (mappings->mappings)[i].chr )
        {
            int ps_loc_index = mappings->mappings[i].start_bp;
            struct pseudo_location_t* ps_loc = ps_locs->locs + ps_loc_index;
            
            /* loop through each location in the pseudo reads */
            int k;
            for( k = 0; k < ps_loc->num; k++ )
            {
                GENOME_LOC_TYPE* gen_locs = ps_loc->locs;
                
                (mappings->mappings)[i].chr = gen_locs[k].chr;
                (mappings->mappings)[i].start_bp = gen_locs[k].loc;

                prob_sum += add_pseudo_loc_to_mapped_read( genome, mappings->mappings + i, mpd_rd );

                /*
                 * if both bit flags are set on a loc in ps_locs->locs,
                 * add a mapped read for the maternal complement
                 */
                if( gen_locs[k].is_paternal && gen_locs[k].is_maternal )
                {
                    candidate_mapping maternal =
                        convert_paternal_candidate_mapping_to_maternal_candidate_mapping(
                            genome, *(mappings->mappings + i)
                        );
                    prob_sum += add_pseudo_loc_to_mapped_read( genome, &maternal, mpd_rd );
                }
                
            }
        } else {
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
    }

    return prob_sum;
}

void
build_mapped_read_from_candidate_mappings( 
    struct genome_data* genome,
    candidate_mappings* mappings, 
    struct mapped_read_t** mpd_rd,
    long read_id )
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
    (*mpd_rd)->free_locations = true;
    (*mpd_rd)->read_id = read_id;
    (*mpd_rd)->rdb = NULL;

    if( mappings->length == 0 )
        return;
    
    /* store the sum of the marginal probabilities */
    double prob_sum;

    /* this assumes all reads are either paired or not */
    if( mappings->mappings[0].rd_type == SINGLE_END )
    {
        prob_sum = build_mapped_read_from_unpaired_candidate_mappings( 
            genome, mappings, *mpd_rd );
    } else {
        prob_sum = build_mapped_read_from_paired_candidate_mappings( 
            genome, mappings, *mpd_rd );
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
init_mapped_reads_db( struct mapped_reads_db** rdb, char* fname, const char* mode )
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

    /* Set by caller (reading or writing) */
    (*rdb)->mode = 0;

    (*rdb)->access_lock = malloc( sizeof(pthread_spinlock_t) );
    pthread_spin_init( (*rdb)->access_lock, PTHREAD_PROCESS_PRIVATE ); 

    /* mmapped data */
    (*rdb)->mmapped_data = NULL;    
    (*rdb)->mmapped_data_size = 0;

    /* index */
    (*rdb)->mmapped_reads_starts = NULL;
    (*rdb)->num_mmapped_reads = 0;
 
    /* fl dist */
    (*rdb)->fl_dist = NULL;

    (*rdb)->num_succ_iterations = NULL;

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
    (*rdb)->num_mapped_reads = *((int*) (*rdb)->mmapped_data);

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

    /* write placeholder for size of mapped reads db */
    int placeholder = 0;
    fwrite( &placeholder, sizeof(int), 1, (*rdb)->fp );
    /* updated on closing the mapped reads db (for writing) */

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
        fprintf( stderr, "Failed to open fl_dist from filename %s\n", filename );
        exit(-1);
    }
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    fclose( fl_fp );
}

void
close_reading_specific_portions_of_mapped_reads_db( struct mapped_reads_db** rdb )
{
    munmap_mapped_reads_db( *rdb );

    if( (*rdb)->fl_dist != NULL ) {
        free_fl_dist( &((*rdb)->fl_dist) );
        (*rdb)->fl_dist = NULL;
    }

    if( (*rdb)->mmapped_reads_starts != NULL ) {
        free( (*rdb)->mmapped_reads_starts );
        (*rdb)->mmapped_reads_starts = NULL;
    }

    return;
}

void
close_writing_specific_portions_of_mapped_reads_db( struct mapped_reads_db** rdb )
{
    /* if the db is open for writing, update the number of mapped reads */
    fseek( (*rdb)->fp, 0, SEEK_SET );
    fwrite( &((*rdb)->num_mapped_reads), sizeof(int), 1, (*rdb)->fp );

    fclose( (*rdb)->fp );
    
    return;
}

void
close_mapped_reads_db( struct mapped_reads_db** rdb )
{
    if( NULL == *rdb )
        return;

    if( (*rdb)->mode == 'r' )
    {
        close_reading_specific_portions_of_mapped_reads_db( rdb );
    } else if( (*rdb)->mode == 'w' ) {
        close_writing_specific_portions_of_mapped_reads_db( rdb );
    } else {
        assert( false );
    }

    pthread_spin_destroy( (*rdb)->access_lock );    
    free( (void*) (*rdb)->access_lock );
    
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
        /* TODO - be able to recover from this */
        exit( -1 );
    }

    int error;

    rd->rdb = rdb;
    
    pthread_spin_lock( rdb->access_lock );
    error = write_mapped_read_to_file( rd, rdb->fp );
    pthread_spin_unlock( rdb->access_lock );

    if( error < 0 )
    {
        fprintf(stderr, "FATAL       :  Error writing to packed mapped reads db.\n");
        exit( -1 );
    }

    rdb->num_mapped_reads += 1;
    
    return;
}

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* if rdb is mmapped */
    rdb->current_read = 0;
    
    return;
}

enum bool
mapped_reads_db_is_empty( struct mapped_reads_db* rdb )
{
    /* if the rdb was opened for reading (mmapped) */
    assert( rdb->mode == 'r' );
    assert( rdb->mmapped_data != NULL );

    if( rdb->current_read == rdb->num_mmapped_reads )
        return true;

    return false;
}

int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* const rdb, 
    struct mapped_read_t** rd )
{
    init_mapped_read( rd );
    (*rd)->rdb = rdb;

    /* Make sure the db is open for reading */
    if( rdb->mode != 'r' )
    {
        fprintf(stderr, "FATAL       :  Cannot get read from mapped reads db unless it is open for reading.\n" );
        exit( -1 );
    }

    /** Get the next read **/
    pthread_spin_lock( rdb->access_lock );
    /* if we have read every read */
    if( rdb->current_read == rdb->num_mmapped_reads )
    {
        pthread_spin_unlock( rdb->access_lock );
        free_mapped_read( *rd );
        *rd = NULL;
        return EOF;
    }
    
    unsigned int current_read_id = rdb->current_read;
    rdb->current_read += 1;
    pthread_spin_unlock( rdb->access_lock );

    assert( current_read_id < rdb->num_mmapped_reads );
    
    /* get a pointer to the current read */
    char* read_start = rdb->mmapped_reads_starts[current_read_id];
    
    /* read a mapping into the struct */
    (*rd)->read_id = *((MPD_RD_ID_T*) read_start);

    read_start += sizeof(MPD_RD_ID_T)/sizeof(char);
    (*rd)->num_mappings = *((MPD_RD_NUM_MAPPINGS_T*) read_start);

    read_start += sizeof(MPD_RD_NUM_MAPPINGS_T)/sizeof(char);

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
        unsigned int j;
        double cond_prob_sum = 0;
        for( j = 0; j < r->num_mappings; j++ )
        {
            int chr_index = get_chr_from_mapped_read_location( r->locations + j );
            unsigned int start = 
                get_start_from_mapped_read_location( r->locations + j );
            unsigned int stop = 
                get_stop_from_mapped_read_location( r->locations + j );
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
    
    /* update the read pointers */
    
    
    
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

    free( rdb->mmapped_reads_starts );
    rdb->mmapped_reads_starts = NULL;
    rdb->num_mmapped_reads = 0;
    
    return;
}

void
index_mapped_reads_db( struct mapped_reads_db* rdb )
{
    const int REALLOC_BLOCK_SIZE = 1000000;

    /* allocate space for the reads start */
    rdb->num_mmapped_reads = 0;

    unsigned long num_allcd_reads = REALLOC_BLOCK_SIZE;
    rdb->mmapped_reads_starts = malloc(sizeof(char*)*REALLOC_BLOCK_SIZE);

    /* Copy the reads data pointer (adding the offset from num_mapped_reads) */
    char* read_start = rdb->mmapped_data + sizeof(int);

    /* Loop through all of the reads */
    while( ((size_t)read_start - (size_t)rdb->mmapped_data) 
           < rdb->mmapped_data_size )
    {
        /* check to ensure the array is big enough */
        if( rdb->num_mmapped_reads + 1 == num_allcd_reads )
        {
            num_allcd_reads += REALLOC_BLOCK_SIZE;
            rdb->mmapped_reads_starts = realloc( 
                rdb->mmapped_reads_starts, num_allcd_reads*sizeof(char*) );
        }
        assert( rdb->num_mmapped_reads < num_allcd_reads );

        /* add the new read start */
        (rdb->mmapped_reads_starts)[rdb->num_mmapped_reads] = read_start;
        (rdb->num_mmapped_reads)++;

        /* read a mapping into the struct */
        /* skip the read ID */
        read_start += sizeof(MPD_RD_ID_T)/sizeof(char);
        MPD_RD_NUM_MAPPINGS_T num_mappings = *((MPD_RD_NUM_MAPPINGS_T*) read_start);
        read_start += sizeof(MPD_RD_NUM_MAPPINGS_T)/sizeof(char);
        /* skip the array of mapped locations */
        read_start += (num_mappings)*(sizeof(struct mapped_read_location)/sizeof(char));
    }

    /* reclaim any wasted memory */
    rdb->mmapped_reads_starts = realloc( rdb->mmapped_reads_starts, 
                                         rdb->num_mmapped_reads*sizeof(char*) );
    
    return;
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

