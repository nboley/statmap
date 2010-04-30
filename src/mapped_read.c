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
    (*rd)->locations = NULL;
    return;
}

void
free_mapped_read( struct mapped_read_t* rd )
{
    if( rd == NULL )
        return;
    free( rd->locations );
    rd->locations = NULL;
    free( rd );
    return;
}

void
add_location_to_mapped_read( 
    struct mapped_read_t* rd, struct mapped_read_location* loc )
{
    /* Allocate new space for the location */
    rd->num_mappings += 1;
    rd->locations = realloc( 
        rd->locations, rd->num_mappings*sizeof(struct mapped_read_location)  );
    if( rd->locations == NULL )
    {
        fprintf(stderr, "FATAL       :  Memory allocation error\n");
        assert( false );
        exit( -1 );
    }

    memcpy( rd->locations + rd->num_mappings - 1, 
            loc, sizeof(struct mapped_read_location) );
    
    return;
}


void
reset_read_cond_probs( struct mapped_read_t* rd  )
{
    /* prevent divide by zero */
    double prb_sum = FLT_MIN;
    int i;
    for( i = 0; i < rd->num_mappings; i++ )
    {
        struct mapped_read_location* loc = rd->locations + i;
        loc->cond_prob = loc->seq_error*loc->fl_prob;
        prb_sum += loc->cond_prob;
    }

    for( i = 0; i < rd->num_mappings; i++ )
    {
        rd->locations[i].cond_prob /= prb_sum;
    }

    return;
};

void
set_read_fl_probs( struct mapped_read_t* rd, 
                   struct fragment_length_dist_t* fl_dist  )
{
    int i;
    for( i = 0; i < rd->num_mappings; i++ )
    {
        struct mapped_read_location* loc = rd->locations + i;
        int fl = loc->stop_pos - loc->start_pos + 1;
        float fl_prb = get_fl_prb( fl_dist, fl );
        loc->fl_prob = fl_prb;
    }
    
    return;
};


void
re_weight_read_cond_prb_by_fl( struct mapped_read_t* r,
                               struct fragment_length_dist_t* fl_dist  )
{
    /* If the fl dist is NULL, assume equal weights for every possible fragment */
    if( NULL != fl_dist )
        return;

    int i;
    /* 
     * Store the total proability mass of paired end reads. 
     * This only matters if there is a mixture of paired and 
     * unpaired reads - ( this should be very rare, to the extent
     * that I'm raising an assertion earlier. But *this* code should
     * work. ) I use it so that single end reads dont have too much 
     * weight in the re-weighting.
     */
    double sum_paired_weight = 0;
    double fl_weights_sum = 0;
    /* For pass one, update the weights */
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* ignore single ended reads */
        if( !((r->locations[i].flag)&IS_PAIRED) )
            continue;
        
        sum_paired_weight += r->locations[i].cond_prob;
        /* we know this is safe because each thread has it's own read */
        double fl_prb = get_fl_prb(  
            fl_dist, r->locations[i].stop_pos - r->locations[i].start_pos + 1  );
        
        r->locations[i].cond_prob *= fl_prb;
        
        fl_weights_sum += r->locations[i].cond_prob;                
    }
    
    /* For pass two, renormalize the weights */
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* ignore single ended reads */
        if( !((r->locations[i].flag)&IS_PAIRED) )
            continue;
        
        r->locations[i].cond_prob /= (fl_weights_sum*sum_paired_weight);
        
        assert( r->locations[i].cond_prob >= -0.0001 );
    }            

};

void
fprintf_mapped_read( FILE* fp, struct mapped_read_t* r )
{
    int i;
    fprintf( fp, "Read ID: %lu\n", r->read_id);
    for( i = 0; i < r->num_mappings; i++ )
    {
        fprintf( fp, "\t%u\t%u\t%u\t%u\t%e\t%e\n",
                 (r->locations[i]).chr,  
                 (r->locations[i]).flag,  
                 (r->locations[i]).start_pos,
                 (r->locations[i]).stop_pos,
                 (r->locations[i]).seq_error,
                 (r->locations[i]).cond_prob
            );
    }
    fprintf( fp, "\n" );
}

int 
write_mapped_read_to_file( struct mapped_read_t* read, FILE* of  )
{
    int num = 0;

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


void
fprintf_nonpaired_mapped_read_as_sam( 
    FILE* sam_fp,
    struct mapped_read_location* pkd_rd,
    struct genome_data* genome,
    char* key, 
    char* seq,
    char* phred_qualities
)
{
    int rd_len = MAX( pkd_rd->start_pos, pkd_rd->stop_pos )
                  - MIN( pkd_rd->start_pos, pkd_rd->stop_pos );

    /* print the query name */
    fprintf( sam_fp, "%s\t", key );

    /* build and print the flag */
    unsigned short flag = 0;    
    if( pkd_rd->start_pos > pkd_rd->stop_pos  )
        flag |= BAM_FREVERSE;
    fprintf( sam_fp, "%hu\t", flag );
    
    /* print the refernce seq ( chr ) name */
    fprintf( sam_fp, "%s\t", genome->chr_names[pkd_rd->chr] );
    
    /* print the genome location */
    /* note that sam expects this in the 5' genome, but we 
       keep it relative to the read strand */
    fprintf( sam_fp, "%u\t", MIN( pkd_rd->start_pos, pkd_rd->stop_pos ) );
    
    /* print the phred quality score */
    /* this is the ( quoting samtools ) 
         'phred scaled posterior probability that the mapping position of this
          read is incorrect'
    */
    fprintf( sam_fp, "%u\t", (unsigned int) MIN( 254, (-10*log10( 1 - pkd_rd->cond_prob)) ) );
    
    /* print the cigar string ( we dont handle indels ) */
    fprintf( sam_fp, "%uM\t", rd_len );

    /* print the mate information - empty since this is not paired */
    fprintf( sam_fp, "0\t0\t*\t" );

    /* print the actual sequence */
    fprintf( sam_fp, "%.*s\t", rd_len, seq );

    /* print the quality string */
    fprintf( sam_fp, "%.*s\t", rd_len, phred_qualities );

    /* Print the optional fields */
    /* For single end reads, we use 
       MQ - Phred likelihood of the read, conditional on 
            the mapping locations being correct 
       
       I take this to mean the phred scaled probabilites of 
       observing the sequence, given that they came from the 
       location in the genome that this mapping refers to.
       
       I also use 
       XQ - which is just the probability of observing this
            sequence, given that they came from this location
       
       XP - the probability that the read came from this location,
            given that it came from somewhere on the genome
    */
    fprintf( sam_fp, "PQ:i:%u\t", (unsigned int) MIN( 254, (-10*log10(1 - pkd_rd->seq_error))) );
    fprintf( sam_fp, "XQ:i:%e\t", pkd_rd->seq_error);
    fprintf( sam_fp, "XP:i:%e\n", pkd_rd->cond_prob);

    return;
}

void
fprintf_paired_mapped_read_as_sam( 
    FILE* sam_fp,
    struct mapped_read_location* pkd_rd,
    struct genome_data* genome,
    char* r1_key, 
    char* r1_seq,
    char* r1_phred_qualities,
    int r1_len,

    char* r2_key, 
    char* r2_seq,
    char* r2_phred_qualities,
    int r2_len
)
{
    assert( pkd_rd->start_pos < pkd_rd->stop_pos );

    /* Print out the read 1 sam line */
    
    int r1_start, r2_start;
    if( (pkd_rd->flag&FIRST_PAIR_IS_FIRST_IN_GENOME) > 0 ) 
    {
        r1_start = pkd_rd->start_pos;
        r2_start = pkd_rd->stop_pos - r2_len;
    } else {
        r1_start = pkd_rd->stop_pos - r1_len;
        r2_start = pkd_rd->start_pos;
    }

    /* print the query name */
    fprintf( sam_fp, "%s\t", r1_key );

    /* build and print the flag */
    unsigned short flag = 0; 
    flag |= BAM_FPAIRED;
    flag |= BAM_FPROPER_PAIR;

    /* Since this is the first read */
    flag |= BAM_FREAD1;

    if( (flag&FIRST_READ_WAS_REV_COMPLEMENTED) > 0  )
        flag |= BAM_FREVERSE;
    fprintf( sam_fp, "%hu\t", flag );
    
    /* print the refernce seq ( chr ) name */
    fprintf( sam_fp, "%s\t", genome->chr_names[pkd_rd->chr] );
    
    /* print the genome location */
    /* note that sam expects this in the 5' genome, but we 
       keep it relative to the read strand */
    fprintf( sam_fp, "%u\t", r1_start );
    
    /* print the phred quality score */
    /* this is the ( quoting samtools ) 
         'phred scaled posterior probability that the mapping position of this
          read is incorrect'
    */
    fprintf( sam_fp, "%u\t", (unsigned int) 
             MIN( 254, (-10*log10(1 - pkd_rd->cond_prob))) );
    
    /* print the cigar string ( we dont handle indels ) */
    fprintf( sam_fp, "%uM\t", r1_len );

    /* print the mate reference name */
    fprintf( sam_fp, "%s\t", r2_key );

    /* print the mate start position */
    fprintf( sam_fp, "%u\t", r2_start );

    /* print the inferred insert size */
    /* This is stupid and screwed up, but quoting from the standard...
    
      If the two reads in a pair are mapped to the same reference, ISIZE 
      equals the difference between the coordinate of the 5ʼ-end of the 
      mate and of the 5ʼ-end of the current read; otherwise ISIZE equals 
      0 (by the “5ʼ-end” we mean the 5ʼ-end of the original read, so for 
      Illumina short-insert paired end reads this calculates the difference 
      in mapping coordinates of the outer edges of the original sequenced 
      fragment). ISIZE is negative if the mate is mapped to a smaller 
      coordinate than the current read.
    
      So I guess this means that I should be aiming for the fragment length, 
      subject of course to the correct sign conventions 
    */
    fprintf( sam_fp, "%i\t", r2_start - r1_start - r1_len + r1_len + r2_len );

    /* print the actual sequence */
    fprintf( sam_fp, "%.*s\t", r1_len, r1_seq );

    /* print the quality string */
    fprintf( sam_fp, "%.*s\t", r1_len, r1_phred_qualities );

    /* Print the optional fields */
    /* For paired end reads, we use 
       PQ - Phred likelihood of the read pair, conditional on both 
            the mapping locations being correct 
       
       I take this to mean the phred scaled probabilites of 
       observing the sequences, given that they came from the 
       location in the genome that this mapping refers to.
       
       I also use 
       XQ - which is just the probability of observing these 
            sequences, given that they came from this location
       
       XP - the probability that the read came from this location,
            given that it came from somewhere on the genome
    */
    fprintf( sam_fp, "PQ:i:%u\t", (unsigned int) 
             MIN( 254, (-10*log10(1-pkd_rd->seq_error))) );
    fprintf( sam_fp, "XQ:i:%e\t", pkd_rd->seq_error);
    fprintf( sam_fp, "XP:i:%e\n", pkd_rd->cond_prob);
   
    /*************************************************************/
    /*************************************************************/
    /*************************************************************/
    /** Print the read's second SAM line */
    
    if( (pkd_rd->flag&FIRST_PAIR_IS_FIRST_IN_GENOME) == 0 ) 
    {
        r1_start = pkd_rd->start_pos;
        r2_start = pkd_rd->stop_pos - r2_len;
    } else {
        r1_start = pkd_rd->stop_pos - r1_len;
        r2_start = pkd_rd->start_pos;
    }

    /* print the query name */
    fprintf( sam_fp, "%s\t", r2_key );

    /* build and print the flag */
    flag = 0; 
    flag |= BAM_FPAIRED;
    flag |= BAM_FPROPER_PAIR;

    /* Since this is the first read */
    flag |= BAM_FREAD1;

    if( (flag&FIRST_READ_WAS_REV_COMPLEMENTED) == 0  )
        flag |= BAM_FREVERSE;
    fprintf( sam_fp, "%hu\t", flag );
    
    /* print the refernce seq ( chr ) name */
    fprintf( sam_fp, "%s\t", genome->chr_names[pkd_rd->chr] );
    
    /* print the genome location */
    /* note that sam expects this in the 5' genome, but we 
       keep it relative to the read strand */
    fprintf( sam_fp, "%u\t", r1_start );
    
    /* print the phred quality score */
    /* this is the ( quoting samtools ) 
         'phred scaled posterior probability that the mapping position of this
          read is incorrect'
    */
    fprintf( sam_fp, "%u\t", (unsigned int) 
             MIN( 254, (-10*log10(1 - pkd_rd->cond_prob))) );
    
    /* print the cigar string ( we dont handle indels ) */
    fprintf( sam_fp, "%uM\t", r2_len );

    /* print the mate reference name */
    fprintf( sam_fp, "%s\t", r1_key );

    /* print the mate start position */
    fprintf( sam_fp, "%u\t", r2_start );

    /* print the inferred insert size */
    /* This is stupid and screwed up, but quoting fromt he standard... */
    /*
      If the two reads in a pair are mapped to the same reference, ISIZE 
      equals the difference between the coordinate of the 5ʼ-end of the 
      mate and of the 5ʼ-end of the current read; otherwise ISIZE equals 
      0 (by the “5ʼ-end” we mean the 5ʼ-end of the original read, so for 
      Illumina short-insert paired end reads this calculates the difference 
      in mapping coordinates of the outer edges of the original sequenced 
      fragment). ISIZE is negative if the mate is mapped to a smaller 
      coordinate than the current read.
    */
    /* So I guess this means that I should be aiming for the fragment length, 
       subject of course to the correct sign conventions 
    */
    fprintf( sam_fp, "%i\t", r2_start - r1_start + r2_len - r1_len - r2_len );

    /* print the actual sequence */
    fprintf( sam_fp, "%.*s\t", r2_len, r2_seq );

    /* print the quality string */
    fprintf( sam_fp, "%.*s\t", r2_len, r2_phred_qualities );

    /* Print the optional fields */
    /* For paired end reads, we use 
       PQ - Phred likelihood of the read pair, conditional on both 
            the mapping locations being correct 
       
       I take this to mean the phred scaled probabilites of 
       observing the sequences, given that they came from the 
       location in the genome that this mapping refers to.
      
       I also use 
       XQ - which is just the probability of observing these 
            sequences, given that they came from this location
       
       XP - the probability that the read came from this location,
            given that it came from somewhere on the genome
    */
    fprintf( sam_fp, "PQ:i:%u\t", (unsigned int) 
             MIN( 254, (-10*log10(pkd_rd->seq_error))) );
    fprintf( sam_fp, "XQ:i:%e\t", pkd_rd->seq_error);
    fprintf( sam_fp, "XP:i:%e\n", pkd_rd->cond_prob);

    return;
}


void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* mpd_rd,
    struct genome_data* genome,
    rawread* rr1,
    rawread* rr2
)
{
    int i = 0;
    for( i = 0; i < mpd_rd->num_mappings; i++ )
    {
        if( rr2 != NULL )
        {
            assert( rr1->end == FIRST  );

            fprintf_paired_mapped_read_as_sam( 
                sam_fp,
                mpd_rd->locations + i,
                genome,

                rr1->name,
                rr1->char_seq,
                rr1->error_str,
                rr1->length,

                rr2->name,
                rr2->char_seq,
                rr2->error_str,
                rr2->length
            );
        } else {
           fprintf_nonpaired_mapped_read_as_sam( 
                sam_fp,
                mpd_rd->locations + i,
                genome,

                rr1->name,
                rr1->char_seq,
                rr1->error_str
            );
        }
    }
    return;
}

void
build_mapped_read_from_candidate_mappings( 
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

    (*mpd_rd)->read_id = read_id;

    /* store the sum of the marginal probabilities */
    double prob_sum = 0;
    
    /* store local read location data */
    struct mapped_read_location loc;
    
    unsigned int i = 0, j = 0, p1_start = -1, p2_start = -1;
    /* Take care of all of the unknown read types ( this should never happen ) */
    while( i < mappings->length
           && mappings->mappings[i].rd_type == UNKNOWN )
    {
        assert( false );
        i++;
    }
    
    while( i < mappings->length
           && mappings->mappings[i].rd_type == SINGLE_END )
    {
        /* if the mapping hasnt been determined to be valid, ignore it */
        if( (mappings->mappings)[i].recheck != VALID )
        {
            i++;
            continue;
        }
        
        /* Add the location */
        /* Ensure all of the flags are turned off */
        loc.flag = 0;
        if( BKWD == (mappings->mappings)[i].rd_strnd )
            loc.flag |= FIRST_READ_WAS_REV_COMPLEMENTED;
        
        if( (mappings->mappings)[i].does_cover_snp )
        {
            loc.flag |= FIRST_READ_COVERS_SNP;
            loc.snps_bm_r1 = (mappings->mappings)[i].snp_bitfield;
        }

        loc.chr = (mappings->mappings)[i].chr;
        loc.start_pos = (mappings->mappings)[i].start_bp;;
        loc.stop_pos = loc.start_pos + (mappings->mappings)[i].rd_len;;
        loc.seq_error = pow( 10, (mappings->mappings)[i].penalty );
        prob_sum += loc.seq_error;
        loc.cond_prob = -1;

        add_location_to_mapped_read( *mpd_rd, &loc );
        
        i++;
    }

    /* If we are done, then return */
    /* ( This is an optimization for a common case, where the 
         reads are all single ) */
    if( i == mappings->length ) {
        goto renormalize_probabilities;
    }

    /* Skip pass the paired ends of indeterminate end */
    /* For now, I can not think of a reason for this to happen, so I forbid it */
    while( i < mappings->length
           && mappings->mappings[i].rd_type == PAIRED_END )
    {
        assert( false );
        i++;
    }

    /* If we are done, return */
    if( i == mappings->length ) {
        goto renormalize_probabilities;
    }
    
    /* Make sure that we are at the begining of the first paired ends */
    p1_start = i;
    while( i < mappings->length
           && mappings->mappings[i].rd_type == PAIRED_END_1 )
    {
        i++;
    }
    
    /* If we are done, return ( note that we dont have any 
       mapped paired end 2's so no paired ends mapped ) */
    if( i == mappings->length ) {
        goto renormalize_probabilities;
    }

    /* Make sure that we are at the beggining of the first paired ends */
    assert( mappings->mappings[i].rd_type == PAIRED_END_2 );
    p2_start = i;

    /* BUG - check the logic here - this actually may not be merge join */
    /* Merge join */
    i = p1_start;
    /* Loop every the inner relation */
    while( i < p2_start )
    {
        /* initialize j to be the start of the second half of the reads */
        j = p2_start;
        
        /* while the read matches are still possible */
        while( j < mappings->length 
               && mappings->mappings[i].chr 
                  == mappings->mappings[j].chr 
            )
        {
            /* Ensure all of the flags are turned off */
            loc.flag = IS_PAIRED;
            
            /* Determine which of the candidate mappings corresponds 
               with the first pair */
            candidate_mapping* first_read = NULL;
            candidate_mapping* second_read = NULL;
            if( PAIRED_END_1 == (mappings->mappings)[i].rd_type ) {
                first_read = mappings->mappings + i;
                second_read = mappings->mappings + j;
            } else {
                assert( PAIRED_END_1 == second_read->rd_type );
                first_read = mappings->mappings + j;
                second_read = mappings->mappings + i;
            }

            /* Set the appropriate flags */
            if( first_read->start_bp < second_read->start_bp )
                loc.flag |= FIRST_PAIR_IS_FIRST_IN_GENOME;
                        
            if( FWD == first_read->rd_strnd )
                loc.flag |= FIRST_READ_WAS_REV_COMPLEMENTED;

            if( first_read->does_cover_snp )
            {
                loc.flag |= FIRST_READ_COVERS_SNP;
                loc.snps_bm_r1 = first_read->snp_bitfield;
            }

            if( second_read->does_cover_snp )
            {
                loc.flag |= SECOND_READ_COVERS_SNP;
                loc.snps_bm_r2 = second_read->snp_bitfield;
            }            
            
            /* Set the chr */
            loc.chr = first_read->chr;
            assert( first_read->chr == second_read->chr);

            if( first_read->start_bp < second_read->start_bp )
            {
                loc.start_pos = first_read->start_bp;
                loc.stop_pos = second_read->start_bp + second_read->rd_len;
            } else {
                loc.start_pos = second_read->start_bp;
                loc.stop_pos = first_read->start_bp + first_read->rd_len;
            }
            
            loc.seq_error = pow( 10, first_read->penalty );
            loc.seq_error *= pow( 10, second_read->penalty );
            
            prob_sum += loc.seq_error;
            loc.cond_prob = -1;

            /* ignore reads with zero probability ( possible with FL dist ) */
            if( loc.seq_error > 2*FLT_MIN )
            {
                add_location_to_mapped_read( *mpd_rd, &loc );
            }

            j++;
        };
                
        i++;
    }

renormalize_probabilities:

    /* set the conditional probabilites - just renormalize */
    for( i = 0; i < (*mpd_rd)->num_mappings; i++ )
    {
        (*mpd_rd)->locations[i].cond_prob 
            = (*mpd_rd)->locations[i].seq_error/prob_sum;

    }

}


/*****************************************************************************
 *
 * Mapped Reads DB Code
 *
 ***************************************************************************/



void
init_mapped_reads_db( struct mapped_reads_db** rdb, char* fname )
{
    *rdb = malloc(sizeof(struct mapped_reads_db));
    (*rdb)->fp = fopen( fname, "w+" );
    if( (*rdb)->fp == NULL )
    {
        perror("FATAL       :  Could not open mapped reads file");
        exit(-1);
    }

    (*rdb)->locked = false;
    
    (*rdb)->mmapped_data = NULL;    
    (*rdb)->mmapped_data_size = 0;

    (*rdb)->mmapped_reads_starts = NULL;
    (*rdb)->num_mmapped_reads = 0;
 
    (*rdb)->fl_dist = NULL;
   
    return;
}

void
build_fl_dist_from_file( struct mapped_reads_db* rdb, FILE* fl_fp )
{
    init_fl_dist_from_file( &(rdb->fl_dist), fl_fp );
    return;
}

void
close_mapped_reads_db( struct mapped_reads_db* rdb )
{
    munmap_mapped_reads_db( rdb );
    fclose( rdb->fp );
    free( rdb );
    return;
}

void
add_read_to_mapped_reads_db( 
    struct mapped_reads_db* rdb,
    struct mapped_read_t* rd)
{
    if ( true == rdb->locked )
    {
        fprintf( stderr, "ERROR       :  Mapped Reads DBis locked - cannot add read.\n");
        /* TODO - be able to recover from this */
        exit( -1 );
    }

    int error;
    
    error = write_mapped_read_to_file( rd, rdb->fp );
    if( error < 0 )
    {
        fprintf( stderr, "FATAL       :  Error writing to packed mapped read db.\n");
        exit( -1 );
    }
    
    return;
}

void
rewind_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if ( true == rdb->locked )
    {
        fprintf( stderr, "ERROR       :  Mapped Reads DBis locked - cannot rewind.\n");
        /* TODO - be able to recover from this */
        exit( -1 );
    }

    assert( rdb->fp != NULL );
    /* make sure the fp is open */
    assert( ftell( rdb->fp ) >= 0 );
    fflush( rdb->fp );
    fseek( rdb->fp, 0L , SEEK_SET );
    // rewind( rdb->fp );
    return;
}

enum bool
mapped_reads_db_is_empty( struct mapped_reads_db* rdb )
{
    if( feof( rdb->fp ) )
        return true;

    return false;
}

int
get_next_read_from_mapped_reads_db( 
    struct mapped_reads_db* rdb, 
    struct mapped_read_t** rd )
{
    if ( true == rdb->locked )
    {
        fprintf( stderr, "ERROR       :  Mapped Reads DBis locked - cannot get next read.\n");
        /* TODO - be able to recover from this */
        exit( -1 );
    }

    size_t rv;

    init_mapped_read( rd );
    
    /* Read in the read id */
    rv = fread( 
        &((*rd)->read_id), sizeof((*rd)->read_id), 1, rdb->fp );
    if( 1 != rv )
    {
        assert( feof( rdb->fp ) );
        free_mapped_read( *rd );
        *rd = NULL;
        return EOF;
    }
    
    /* Read in the number of mappings */
    rv = fread(
        &((*rd)->num_mappings), sizeof((*rd)->num_mappings), 1, rdb->fp );
    if( 1 != rv )
    {
        fprintf( stderr, "FATAL       :  Unexpected end of file\n" );
        assert( false );
        exit( -1 );
    }

    /* read in the locations */
    (*rd)->locations = malloc( 
        (*rd)->num_mappings*sizeof(struct mapped_read_location) );
    
    rv = fread( (*rd)->locations, 
                sizeof(struct mapped_read_location), 
                (*rd)->num_mappings, 
                rdb->fp );
    
    if( (*rd)->num_mappings != rv )
    {
        fprintf( stderr, "FATAL       :  Unexpected end of file\n" );
        assert( false );
        exit( -1 );
    }

    return 0;
}

void
set_all_read_fl_probs( struct mapped_reads_db* rdb )
{
    long long i;
    for( i = 0; i < (long long) rdb->num_mmapped_reads; i++ )
    {
        char* read_start = rdb->mmapped_reads_starts[i];
        
        /* read a mapping into the struct */
        struct mapped_read_t r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (struct mapped_read_location*) read_start;

        set_read_fl_probs( &r, rdb->fl_dist );
    }
}

void
reset_all_read_cond_probs( struct mapped_reads_db* rdb )
{
    long long i;
    for( i = 0; i < (long long) rdb->num_mmapped_reads; i++ )
    {
        char* read_start = rdb->mmapped_reads_starts[i];
        
        /* read a mapping into the struct */
        struct mapped_read_t r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (struct mapped_read_location*) read_start;

        reset_read_cond_probs( &r );
    }
}


void
write_mapped_reads_to_sam( struct rawread_db_t* rdb,
                           struct mapped_reads_db* mappings_db,
                           struct genome_data* genome,
                           FILE* sam_ofp )
{
    int error;

    long readkey;
    
    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );
    
    rawread *rd1, *rd2;
    struct mapped_read_t* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    
    get_next_read_from_rawread_db( 
        rdb, &readkey, &rd1, &rd2 );

    while( rd1 != NULL ) 
    {         
        /* 
         * TODO - clean this up - this is kind of a weird place to put this. 
         *
         * If we couldnt map it anywhere,
         * print out the read to the non-mapping
         * fastq file. ( Theoretically, one could
         * rerun these with lower penalties, or 
         * re-align these to one another. )
         */
        /* We test for mapped read NULL in case the last read was unmappable */
        if( mapped_rd == NULL 
            || mapped_rd->num_mappings == 0 )
        {
            /* If this is a single end read */
            if( rd2 == NULL )
            {
                if( filter_rawread( rd1 ) )
                {
                    fprintf_rawread_to_fastq( 
                        rdb->unmappable_single_end_reads, rd1 );
                } else {
                    fprintf_rawread_to_fastq( 
                        rdb->non_mapping_single_end_reads, rd1 );                    
                }
            } else {
                if( filter_rawread( rd1 ) || filter_rawread( rd2 ) )
                {
                    fprintf_rawread_to_fastq( 
                        rdb->unmappable_paired_end_1_reads, rd1 );
                    
                    fprintf_rawread_to_fastq( 
                        rdb->unmappable_paired_end_2_reads, rd2 );
                } else {
                    fprintf_rawread_to_fastq( 
                        rdb->non_mapping_paired_end_1_reads, rd1 );
                    
                    fprintf_rawread_to_fastq( 
                        rdb->non_mapping_paired_end_2_reads, rd2 );
                }
            }
        /* otherwise, print it out to the sam file */
        } else {
            fprintf_mapped_read_to_sam( 
                sam_ofp, mapped_rd, genome, rd1, rd2 );
        }
        
        free_mapped_read( mapped_rd );
        
        /* Free the raw reads */
        free_rawread( rd1 );
        if( rd2 != NULL )
            free_rawread( rd2 );

        error = get_next_read_from_mapped_reads_db( 
            mappings_db, 
            &mapped_rd
        );

        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd1, &rd2 );
    }

    goto cleanup;

cleanup:
    free_mapped_read( mapped_rd );

    /* Free the raw reads */
    if( rd1 != NULL )
        free_rawread( rd1 );
    if( rd2 != NULL )
        free_rawread( rd2 );
        
    return;
}

/* use this for wiggles */
void
update_traces_from_read_densities( 
    struct mapped_reads_db* reads_db,
    struct trace_t* traces
)
{    
    int i;
    for( i = 0; i < traces->num_chrs; i++ )
    {
        memset( traces->traces[0] + i, 0, 
                sizeof(TRACE_TYPE)*(traces->trace_lengths[i]) );
    }
    
    /* Update the trace from the reads */
    for( i = 0; i < (long long) reads_db->num_mmapped_reads; i++ )
    {
        char* read_start = reads_db->mmapped_reads_starts[i];

        /* read a mapping into the struct */
        struct mapped_read_t r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (struct mapped_read_location*) read_start;
        
        /* Update the trace from this mapping */
        unsigned int j;
        double cond_prob_sum = 0;
        for( j = 0; j < r.num_mappings; j++ )
        {
            int chr_index = r.locations[j].chr;
            unsigned int start = r.locations[j].start_pos;
            unsigned int stop = r.locations[j].stop_pos;
            ML_PRB_TYPE cond_prob = r.locations[j].cond_prob;
            cond_prob_sum += cond_prob;
            
            assert( cond_prob >= -0.0001 );
            assert( stop >= start );            
            assert( chr_index < traces->num_chrs );
            assert( traces->trace_lengths[chr_index] >= stop );

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
write_mapped_reads_to_wiggle( struct mapped_reads_db* rdb, 
                              struct genome_data* genome,
                              FILE* wfp )
{
    const double filter_threshold = 1e-6;

    /* Print out the header */
    fprintf( wfp, "track type=wiggle_0 name=%s\n", "cage_wig_track" );

    /* build and update the chr traces */
    struct trace_t* traces;
    init_traces( genome, &traces, 1 );

    update_traces_from_mapped_reads( 
        rdb, traces,
        naive_update_trace_expectation_from_location
    );
    
    int i;
    unsigned int j;
    for( i = 0; i < traces->num_chrs; i++ )
    {
        fprintf( wfp, "variableStep chrom=%s\n", genome->chr_names[i] );
        
        for( j = 0; j < traces->trace_lengths[i]; j++ )
        {
            if( traces->traces[0][i][j] >= filter_threshold )
                fprintf( wfp, "%i\t%e\n", j+1, traces->traces[0][i][j] );
        }
    }
    
    close_traces( traces );
}




void
mmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    /* get the file descriptor for the file we wish to mmap */
    int fdin = fileno( rdb->fp );
    
    /* Lock the db to prevent writes */
    rdb->locked = true;
    
    /* make sure the entire file has been written to disk */
    fflush( rdb->fp );

    /* find the size of the opened file */
    struct stat buf;
    fstat(fdin, &buf);
    rdb->mmapped_data_size = buf.st_size;
        
    /* mmap the file */
    rdb->mmapped_data
        = mmap( NULL, rdb->mmapped_data_size,  
                PROT_READ|PROT_WRITE, 
		MAP_POPULATE|MAP_SHARED, fdin, (off_t) 0 );

    if( rdb->mmapped_data == (void*) -1 )
    {
        char* buffer;
        buffer = malloc( sizeof(char)*500 );
        sprintf(buffer, "Can not mmap the fdescriptor '%u'", fdin );
	perror( buffer );
        assert( false );
        exit( -1 );
    }

    return;
}

void
munmap_mapped_reads_db( struct mapped_reads_db* rdb )
{
    if( rdb->mmapped_data == NULL )
        return;

    int error = munmap( rdb->mmapped_data, rdb->mmapped_data_size );
    if( error != 0 )
    {
        perror( "Could not munmap mapped reads db" );
        assert( false );
        exit( -1 );
    }
    rdb->mmapped_data = NULL;

    rdb->locked = false;
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

    /* Copy the reads data pointer */
    char* read_start = rdb->mmapped_data;

    /* Loop therough all of the reads */
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
        read_start += sizeof(unsigned long)/sizeof(char);
        unsigned short num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        /* skip the array of mapped locations */
        read_start += (num_mappings)*(sizeof(struct mapped_read_location)/sizeof(char));
    }

    /* reclaim any wasted memory */
    rdb->mmapped_reads_starts = realloc( rdb->mmapped_reads_starts, 
                                         rdb->num_mmapped_reads*sizeof(char*) );
    
    return;
}

/*
 *  END Mapped Reads
 *
 **************************************************************************/

