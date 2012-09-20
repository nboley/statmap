#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "sam.h"

#include "genome.h"
#include "rawread.h"
#include "mapped_read.h"
#include "pseudo_location.h"
#include "dna_sequence.h"

void
fprintf_sam_headers_from_genome(
        FILE* sam_fp,
        struct genome_data* genome
    )
{
    /**
     * Write @SQ tag for each contig in the genome
     **/
    int chr;
    /* loop over contigs; chr=1 to skip pseudo chr */
    for( chr = 1; chr < genome->num_chrs; chr++ )
    {
        fprintf(sam_fp, "@SQ\tSN:%s\tLN:%d\n",
                genome->chr_names[chr], genome->chr_lens[chr]
            );
    }
}

void
fprintf_cigar_string(
        FILE* sam_fp,
        MRL_TRIM_TYPE trim_offset,
        int rd_len,
        MRL_FLAG_TYPE flag
    )
{
    /* If the read mapped to the reverse strand, the cigar string needs to
     * correspond to the reverse complemented sequence */
    if( 0 < (flag&FIRST_READ_WAS_REV_COMPLEMENTED) )
    {
        if( trim_offset > 0 )
        {
            fprintf( sam_fp, "%uM%uS\t", rd_len - trim_offset, trim_offset );
        } else {
            fprintf( sam_fp, "%uM\t", rd_len );
        }
    }
    else
    {
        /* if there is a trim offset, soft clip the initial trim_offset bases */
        if( trim_offset > 0 )
        {
            fprintf( sam_fp, "%uS%uM\t", trim_offset, rd_len - trim_offset );
        } else {
            fprintf( sam_fp, "%uM\t", rd_len );
        }
    }
}

void
fprintf_seq(
        FILE* sam_fp,
        char* seq,
        int rd_len,
        MRL_FLAG_TYPE flag
    )
{
    if( 0 < (flag&FIRST_READ_WAS_REV_COMPLEMENTED) )
    {
        /* print the reverse complemented sequence for reads that map to the
         * reverse strand */
        char* tmp_seq = calloc( rd_len + 1, sizeof(char) );

        rev_complement_read( seq, tmp_seq, rd_len );
        fprintf( sam_fp, "%.*s\t", rd_len, tmp_seq );

        free( tmp_seq );
    }
    else
    {
        fprintf( sam_fp, "%.*s\t", rd_len, seq );
    }
}

void
fprintf_nonpaired_mapped_read_as_sam( 
    FILE* sam_fp,
    struct mapped_read_location* loc,
    float cond_prob,
    struct genome_data* genome,
    char* key, 
    char* seq,
    char* phred_qualities
)
{
    const MRL_FLAG_TYPE loc_flag = get_flag_from_mapped_read_location( loc  );
    const int chr_index = get_chr_from_mapped_read_location( loc  );
    /* + 1 because Statmap's is 0-indexed, but SAM is 1-indexed */
    unsigned int start = get_start_from_mapped_read_location( loc  ) + 1;
    unsigned int stop = get_stop_from_mapped_read_location( loc  ) + 1;
    const float seq_error = get_seq_error_from_mapped_read_location( loc  );

    if( loc_flag&FIRST_READ_WAS_REV_COMPLEMENTED  ) {
        stop -= loc->rd1_trim_offset;
    } else{ 
        start += loc->rd1_trim_offset;
    }


    assert( cond_prob <= 1.0 );
    
    int rd_len = MAX( start, stop )
                  - MIN( start, stop );

    /* print the query name */
    fprintf( sam_fp, "%s\t", key );

    /* build and print the flag */
    unsigned short flag = 0;    
    if( 0 < (loc_flag&FIRST_READ_WAS_REV_COMPLEMENTED)  )
        flag |= BAM_FREVERSE;
    fprintf( sam_fp, "%hu\t", flag );
    
    /* print the refernce seq ( chr ) name */
    fprintf( sam_fp, "%s\t", genome->chr_names[chr_index] );
    
    /* print the genome location */
    /* note that sam expects this in the 5' genome, but we 
       keep it relative to the read strand */
    fprintf( sam_fp, "%u\t", MIN( start, stop ) );
    
    /* print the phred quality score */
    /* this is the ( quoting samtools ) 
         'phred scaled posterior probability that the mapping position of this
          read is incorrect'
    */
    fprintf( sam_fp, "%u\t", 
             (unsigned int) MIN( 254, (-10*log10( 1 - cond_prob)) ) );
    
    /* print the cigar string ( we dont handle indels ) */
    /* soft clip to account for assay specific corrections */
    fprintf_cigar_string( sam_fp, loc->rd1_trim_offset, rd_len, loc_flag );

    /* print the mate information - empty since this is not paired */
    fprintf( sam_fp, "*\t0\t0\t" );

    /* print the actual sequence */
    fprintf_seq( sam_fp, seq, rd_len, loc_flag );
    
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
    fprintf( sam_fp, "PQ:i:%u\t", (unsigned int) MIN(254, (-10*log10(1-seq_error))));
    fprintf( sam_fp, "XQ:f:%f\t", seq_error);
    fprintf( sam_fp, "XP:f:%f\n", cond_prob);

    return;
}

void
fprintf_paired_mapped_read_as_sam( 
    FILE* sam_fp,
    struct mapped_read_location* loc,
    const float cond_prob,

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
    const unsigned char mrl_flag = get_flag_from_mapped_read_location( loc  );
    const int chr_index = get_chr_from_mapped_read_location( loc  );
    /* + 1 because Statmap's is 0-indexed, but SAM is 1-indexed */
    const unsigned int start = get_start_from_mapped_read_location( loc  ) + 1;
    const unsigned int stop = get_stop_from_mapped_read_location( loc  ) + 1;
    const float seq_error = get_seq_error_from_mapped_read_location( loc  );


    if( cond_prob > 1.0 )
    {
        printf( "ERROR: %e\t%e\n", cond_prob, seq_error );
    }
    assert( cond_prob <= 1.0 );
    
    assert( start < stop );

    /* Print out the read 1 sam line */
    
    int r1_start, r2_start;
    if( (mrl_flag&FIRST_PAIR_IS_FIRST_IN_GENOME) > 0 ) 
    {
        r1_start = start;
        r2_start = stop - r2_len;
    } else {
        r1_start = stop - r1_len;
        r2_start = start;
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
    fprintf( sam_fp, "%s\t", genome->chr_names[chr_index] );
    
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
             MIN( 254, (-10*log10(1 - cond_prob))) );
    
    /* print the cigar string ( we dont handle indels ) */
    /* soft clip to account for assay specific corrections */
    fprintf_cigar_string( sam_fp, loc->rd1_trim_offset, r1_len, mrl_flag );

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
    fprintf_seq( sam_fp, r1_seq, r1_len, mrl_flag );

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
             MIN( 254, (-10*log10(1-seq_error))) );
    fprintf( sam_fp, "XQ:f:%f\t", seq_error);
    fprintf( sam_fp, "XP:f:%f\n", cond_prob);
   
    /*************************************************************/
    /*************************************************************/
    /*************************************************************/
    /** Print the read's second SAM line */
    
    if( (mrl_flag&FIRST_PAIR_IS_FIRST_IN_GENOME) == 0 ) 
    {
        r1_start = start;
        r2_start = stop - r2_len;
    } else {
        r1_start = stop - r1_len;
        r2_start = start;
    }

    /* print the query name */
    fprintf( sam_fp, "%s\t", r2_key );

    /* build and print the flag */
    flag = 0; 
    flag |= BAM_FPAIRED;
    flag |= BAM_FPROPER_PAIR;

    /* Since this is the first read */
    flag |= BAM_FREAD2;

    if( (flag&FIRST_READ_WAS_REV_COMPLEMENTED) == 0  )
        flag |= BAM_FREVERSE;
    fprintf( sam_fp, "%hu\t", flag );
    
    /* print the refernce seq ( chr ) name */
    fprintf( sam_fp, "%s\t", genome->chr_names[chr_index] );
    
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
             MIN( 254, (-10*log10(1 - cond_prob))) );
    
    /* print the cigar string ( we dont handle indels ) */
    fprintf( sam_fp, "%uM\t", r2_len );
    /* print the cigar string ( we dont handle indels ) */
    /* soft clip to account for assay specific corrections */
    fprintf_cigar_string( sam_fp, loc->rd2_trim_offset, r2_len, mrl_flag );

    /* print the mate reference name */
    fprintf( sam_fp, "%s\t", r1_key );

    /* print the mate start position */
    fprintf( sam_fp, "%u\t", r2_start );

    /* print the inferred insert size */
    /* This is stupid and screwed up, but quoting from the standard... */
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
    fprintf_seq( sam_fp, r2_seq, r2_len, mrl_flag );

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
             MIN( 254, (-10*log10(seq_error))) );
    fprintf( sam_fp, "XQ:f:%f\t", seq_error);
    fprintf( sam_fp, "XP:f:%f\n", cond_prob);

    return;
}


void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* mpd_rd,
    struct cond_prbs_db_t* cond_prbs_db,    
    struct genome_data* genome,
    struct rawread* rr1,
    struct rawread* rr2,
    enum bool expand_pseudo_locations
)
{
    assert( expand_pseudo_locations == false );

    size_t i = 0;
    for( i = 0; i < mpd_rd->num_mappings; i++ )
    {
        float cond_prob = get_cond_prb( cond_prbs_db, mpd_rd->read_id, i );
        
        /* if this is a paired end read */
        if( rr2 != NULL )
        {
            /* make sure the flag agrees */
            assert( rr1->end == FIRST  );

            fprintf_paired_mapped_read_as_sam( 
                sam_fp,
                mpd_rd->locations + i,
                cond_prob,

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
                cond_prob,
                
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
write_nonmapping_reads_to_fastq( 
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mappings_db
)
{
    int error;

    readkey_t readkey;
    
    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );
    
    struct rawread *rd1, *rd2;
    struct mapped_read_t* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    
    get_next_read_from_rawread_db( 
        rdb, &readkey, &rd1, &rd2, -1 );

    while( rd1 != NULL ) 
    {    
        /* if this read doesn't have an associated mapped reads */
        if( mapped_rd == NULL
            || ( mapped_rd != NULL
                   && (
                        readkey < mapped_rd->read_id
                        || 0 == mapped_rd->num_mappings  
                   )
                )
            )
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
        }
        
        /* if we need to get the next mapped read */
        /* if mapped_rd is null, we are out of mapped reads */
        if( mapped_rd != NULL
            && readkey == mapped_rd->read_id )
        {
            free_mapped_read( mapped_rd );
            error = get_next_read_from_mapped_reads_db( 
                mappings_db, 
                &mapped_rd
            );
        }

        /* Free the raw reads */
        free_rawread( rd1 );
        if( rd2 != NULL )
            free_rawread( rd2 );
        
        /* we always get the next raw read */
        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd1, &rd2, -1 );
    }

    goto cleanup;

cleanup:
    if( NULL != mapped_rd )
        free_mapped_read( mapped_rd );

    /* Free the raw reads */
    if( rd1 != NULL )
        free_rawread( rd1 );
    if( rd2 != NULL )
        free_rawread( rd2 );

    fflush( rdb->unmappable_single_end_reads );
    fflush( rdb->non_mapping_single_end_reads );
    fflush( rdb->unmappable_paired_end_1_reads );
    fflush( rdb->unmappable_paired_end_2_reads );
    fflush( rdb->non_mapping_paired_end_1_reads );
    fflush( rdb->non_mapping_paired_end_2_reads );
    
    return;
}


void
write_mapped_reads_to_sam( 
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mappings_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct genome_data* genome,
    enum bool reset_cond_read_prbs,
    /* whether or not to print out pseudo locations
       as real locations, or to print out each real loc
       that makes ups the pseudo location */
    enum bool expand_pseudo_locations,
    FILE* sam_ofp )
{
    assert( expand_pseudo_locations == false );

    /* first, write the SAM header out */
    fprintf_sam_headers_from_genome( sam_ofp, genome );

    int error;

    readkey_t readkey;
    
    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );
    
    struct rawread *rd1, *rd2;
    struct mapped_read_t* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );

    get_next_read_from_rawread_db( 
        rdb, &readkey, &rd1, &rd2, -1 );

    while( NULL != rd1 
           && NULL != mapped_rd
           && readkey < mapped_rd->read_id )
    {
        free_rawread( rd1 );
        rd1 = NULL;
        if( rd2 !=  NULL ) {
            free_rawread( rd2 );
            rd2 = NULL;
        }
        
        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd1, &rd2, -1 );
    }

    while( rd1 != NULL
           && mapped_rd != NULL ) 
    {    
        assert( readkey == mapped_rd->read_id );
        
        if( readkey > 0 && readkey%1000000 == 0 )
            fprintf( stderr, "NOTICE       : Written %u reads to sam\n", readkey );
        
        /* We test for mapped read NULL in case the last read was unmappable */
        if( mapped_rd != NULL 
            && mapped_rd->num_mappings > 0 )
        {
            /* sometimes we want the marginal distribution */
            if( reset_cond_read_prbs )
                reset_read_cond_probs( cond_prbs_db, mapped_rd );

            fprintf_mapped_read_to_sam( 
                sam_ofp, mapped_rd, cond_prbs_db, 
                genome, rd1, rd2, expand_pseudo_locations );
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
            rdb, &readkey, &rd1, &rd2, -1 );
        
        while( NULL != rd1 
               && NULL != mapped_rd
               && readkey < mapped_rd->read_id )
        {
            free_rawread( rd1 );
            rd1 = NULL;
            if( rd2 !=  NULL ) {
                free_rawread( rd2 );
                rd2 = NULL;
            }
            
            get_next_read_from_rawread_db( 
                rdb, &readkey, &rd1, &rd2, -1 );
        }
        
    }
    
    goto cleanup;

cleanup:
    if( NULL != mapped_rd )
        free_mapped_read( mapped_rd );

    /* Free the raw reads */
    if( rd1 != NULL )
        free_rawread( rd1 );
    if( rd2 != NULL )
        free_rawread( rd2 );
        
    return;
}
