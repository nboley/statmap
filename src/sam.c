#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "sam.h"

#include "genome.h"
#include "rawread.h"
#include "mapped_read.h"

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
    struct rawread* rr1,
    struct rawread* rr2
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
write_nonmapping_reads_to_fastq( 
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mappings_db
)
{
    int error;

    long readkey;
    
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
        rdb, &readkey, &rd1, &rd2 );

    while( rd1 != NULL ) 
    {         
        /* 
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


void
write_mapped_reads_to_sam( struct rawread_db_t* rdb,
                           struct mapped_reads_db* mappings_db,
                           struct genome_data* genome,
                           enum bool reset_cond_read_prbs,
                           FILE* sam_ofp )
{
    int error;

    long readkey;
    
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
        rdb, &readkey, &rd1, &rd2 );

    while( rd1 != NULL ) 
    {    
        if( readkey%100000 == 0 )
            fprintf( stderr, "NOTICE       : Written %li reads to sam\n", readkey );
        
        /* We test for mapped read NULL in case the last read was unmappable */
        if( mapped_rd != NULL 
            && mapped_rd->num_mappings > 0 )
        {
            if( reset_cond_read_prbs )
                reset_read_cond_probs( mapped_rd );

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
