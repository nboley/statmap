#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "sam.h"

#include "config.h"
#include "genome.h"
#include "read.h"
#include "mapped_read.h"
#include "pseudo_location.h"
#include "error_correction.h"

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
    const unsigned int start = get_start_from_mapped_read_location( loc  );
    const unsigned int stop = get_stop_from_mapped_read_location( loc  );
    const float seq_error = get_seq_error_from_mapped_read_location( loc  );

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
    fprintf( sam_fp, "PQ:i:%u\t", (unsigned int) MIN(254, (-10*log10(1-seq_error))));
    fprintf( sam_fp, "XQ:i:%e\t", seq_error);
    fprintf( sam_fp, "XP:i:%e\n", cond_prob);

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
    const unsigned int start = get_start_from_mapped_read_location( loc  );
    const unsigned int stop = get_stop_from_mapped_read_location( loc  );
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
             MIN( 254, (-10*log10(1-seq_error))) );
    fprintf( sam_fp, "XQ:i:%e\t", seq_error);
    fprintf( sam_fp, "XP:i:%e\n", cond_prob);
   
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
    flag |= BAM_FREAD1;

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
             MIN( 254, (-10*log10(seq_error))) );
    fprintf( sam_fp, "XQ:i:%e\t", seq_error);
    fprintf( sam_fp, "XP:i:%e\n", cond_prob);

    return;
}


void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    struct mapped_read_t* mpd_rd,
    struct cond_prbs_db_t* cond_prbs_db,    
    struct genome_data* genome,
    struct read* r,
    enum bool expand_pseudo_locations
)
{
    assert( expand_pseudo_locations == false );

    /* HACK - assumptions to get this to compile */
    assert( r->num_subtemplates == 1 || r->num_subtemplates == 2 );

    MPD_RD_ID_T i = 0;
    for( i = 0; i < mpd_rd->num_mappings; i++ )
    {
        float cond_prob = get_cond_prb( cond_prbs_db, mpd_rd->read_id, i );
        
        /* if this is a paired end read */
        if( r->num_subtemplates == 2 )
        {
            // reference to read subtemplates
            struct read_subtemplate* r1 = &(r->subtemplates[0]);
            struct read_subtemplate* r2 = &(r->subtemplates[1]);

            /* make sure the flag agrees */
            assert( r1->pos_in_template.pos == POS_PAIRED_END_1 );
            assert( r2->pos_in_template.pos == POS_PAIRED_END_2 );

            fprintf_paired_mapped_read_as_sam( 
                sam_fp,
                mpd_rd->locations + i,
                cond_prob,

                genome,

                r->name,
                r1->char_seq,
                r1->error_str,
                r1->length,

                r->name, // XXX - will this work? Gets rid of /0 and /1!
                r2->char_seq,
                r2->error_str,
                r2->length
            );
        } else {
            // reference to first read subtemplate
            struct read_subtemplate* r1 = &(r->subtemplates[0]);
            assert( r1->pos_in_template.pos == POS_SINGLE_END );

            fprintf_nonpaired_mapped_read_as_sam( 
                sam_fp,
                mpd_rd->locations + i,
                cond_prob,
                
                genome,
                
                r->name,
                r1->char_seq,
                r1->error_str
            );
        }
    }
    return;
}

/**
 * Writes a single (single or paired end) read out to the correct file pointer(s)
 */
void
write_nonmapping_read_to_fastq(
        struct read* r,
        FILE* single_end_reads_fp,
        FILE* paired_end_1_reads_fp,
        FILE* paired_end_2_reads_fp
    )
{
    /* HACK - assumptions to get this to compile */
    assert( r->num_subtemplates == 1 || r->num_subtemplates == 2 );

    /* If this is a single end read */
    if( r->num_subtemplates == 1 )
    {
        // reference to read subtemplate
        struct read_subtemplate* r1 = &(r->subtemplates[0]);

        // check position in template
        assert( r1->pos_in_template.pos == POS_SINGLE_END );

        fprintf_read_subtemplate_to_fastq( single_end_reads_fp, r->name, r1 );
    }
    /* if this is a paired end read */
    else {
        // references to read subtemplates
        struct read_subtemplate* r1 = &(r->subtemplates[0]);
        struct read_subtemplate* r2 = &(r->subtemplates[1]);

        // check position in template
        assert( r1->pos_in_template.pos == POS_PAIRED_END_1 );
        assert( r2->pos_in_template.pos == POS_PAIRED_END_2 );

        fprintf_read_subtemplate_to_fastq( paired_end_1_reads_fp, r->name, r1 );
        fprintf_read_subtemplate_to_fastq( paired_end_2_reads_fp, r->name, r2 );
    }
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
    
    struct read* rd;
    struct mapped_read_t* mapped_rd;
    
    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    assert( error == 0 );
    
    get_next_read_from_rawread_db( 
        rdb, &readkey, &rd, -1 );

    while( rd != NULL )
    {
        /* if this read doesn't have an associated mapped reads */
        if( mapped_rd == NULL
            || ( mapped_rd != NULL
                   &&
                 readkey < mapped_rd->read_id )
          )
        {
            /* then it was unmappable, and was never added to the mpd rd db */
            write_nonmapping_read_to_fastq( rd,
                    rdb->unmappable_single_end_reads,
                    rdb->unmappable_paired_end_1_reads,
                    rdb->unmappable_paired_end_2_reads
                );
        }
        /* if the read has an associated mapped reads */
        else if( mapped_rd != NULL &&
                 readkey == mapped_rd->read_id &&
                 mapped_rd->num_mappings == 0 )
        {
            /*
               If the read has an associate mapped read, but the mapped read
               has no mapping, it was declared mappable but did not map.

               Nonmapping reads are added to the mapped reads db; see
               find_candidate_mappings.c:1028
            */
            write_nonmapping_read_to_fastq( rd,
                    rdb->non_mapping_single_end_reads,
                    rdb->non_mapping_paired_end_1_reads,
                    rdb->non_mapping_paired_end_2_reads
                );
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

        /* Free the read */
        free_read( rd );
        
        /* we always get the next read */
        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd, -1 );

    }

    goto cleanup;

cleanup:
    if( NULL != mapped_rd )
        free_mapped_read( mapped_rd );

    /* Free the read */
    free_read( rd );

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

    int error;

    readkey_t readkey;

    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );

    struct read* rd;
    struct mapped_read_t* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    assert( error == 0 );

    get_next_read_from_rawread_db( 
        rdb, &readkey, &rd, -1 );

    while( NULL != rd 
           && NULL != mapped_rd
           && readkey < mapped_rd->read_id )
    {
        free_read( rd );
        
        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd, -1 );
    }

    while( rd != NULL
           && mapped_rd != NULL ) 
    {    
        if( readkey != mapped_rd->read_id )
        {
            fprintf( stderr, "readkey: %d\n", readkey );
            fprintf( stderr, "mapped_rd->read_id: %d\n", mapped_rd->read_id );
        }
        assert( readkey == mapped_rd->read_id );
        
        if( readkey > 0 && readkey%1000000 == 0 )
            fprintf( stderr, "NOTICE       : Written %u reads to sam\n", readkey );
        
        /* We test for mapped read NULL in case the last read was unmappable */
        if( mapped_rd != NULL 
            && mapped_rd->num_mappings > 0 )
        {
            /* sometimes we want the marginal distribution */
            if( reset_cond_read_prbs )
                reset_read_cond_probs( cond_prbs_db, mapped_rd, mappings_db );

            fprintf_mapped_read_to_sam( 
                sam_ofp, mapped_rd, cond_prbs_db, 
                genome, rd, expand_pseudo_locations );
        }
        
        free_mapped_read( mapped_rd );
        
        /* Free the read */
        free_read( rd );

        error = get_next_read_from_mapped_reads_db( 
            mappings_db, 
            &mapped_rd
        );

        get_next_read_from_rawread_db( 
            rdb, &readkey, &rd, -1 );
        
        while( NULL != rd 
               && NULL != mapped_rd
               && readkey < mapped_rd->read_id )
        {
            free_read( rd );
            
            get_next_read_from_rawread_db( 
                rdb, &readkey, &rd, -1 );
        }
        
    }
    
    goto cleanup;

cleanup:
    if( NULL != mapped_rd )
        free_mapped_read( mapped_rd );

    /* Free the read */
    if( rd != NULL )
        free_read( rd );
        
    return;
}
