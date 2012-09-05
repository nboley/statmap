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

/*************************************************************************
 *
 *  Write mapped reads to SAM
 * 
 */

enum bool
last_sublocation_in_mapped_read_location( mapped_read_sublocation* sub )
{
    if( !sub->next_subread_is_ungapped &&
        !sub->next_subread_is_gapped )
    {
        /* then there is no "next" subread, and this is the last sublocation in
         * the mapped_read_location */
        return true;
    }

    return false;
}

enum bool
end_of_sublocations_group( mapped_read_sublocation* sub_loc )
{
    /* Sublocations are grouped by consecutive gapped sublocations.
     *
     * If a sublocation is either
     * 1) The last sublocation in a group of gapped sublocations
     * 2) The last sublocation in a mapped_read_location
     *
     * it is at the end of a group of sublocations.
     */
    if( sub_loc->next_subread_is_ungapped ||
            last_sublocation_in_mapped_read_location( sub_loc ))
        return true;

    return false;
}

char*
build_cigar_string_for_sublocations_group(
        mapped_read_sublocation* group_start )
{
    /* Note - this function allocates dynamic memory. It is the responsibility
     * of the caller to free it. */

    /* temporary buffer for building the cigar string */
    char buf[512];
    char* buf_pos = buf;

    /* ptr to iterate through the sublocations */
    char* ptr = (char*) group_start;

    /* Stop position of previous sublocation in the group. Initialized to -1 for
     * the first sublocation, since it has no prev_stop */
    int prev_stop = -1;

    /* loop over the sublocations */
    while(true)
    {
        mapped_read_sublocation* current_sublocation
            = (mapped_read_sublocation*) ptr;

        if( prev_stop > 0 )
        {
            /* If prev_stop is a real position, represent the intron gap
             * in the CIGAR string with N */
            assert( current_sublocation->start_pos > prev_stop );
            buf_pos += sprintf( buf_pos, "%iN", 
                    current_sublocation->start_pos - prev_stop );
        }

        /* Write the region of aligned sequence with M */
        buf_pos += sprintf( buf_pos, "%iM", current_sublocation->length );

        /* Save prev_stop for the next iteration */
        prev_stop = current_sublocation->start_pos + current_sublocation->length;

        /* If the next sublocation is ungapped, or we're at the end of the
         * sublocations, we are done building the cigar string */
        if( end_of_sublocations_group( current_sublocation ) )
            break;

        /* move ptr to the next sublocation */
        ptr += sizeof( mapped_read_sublocation );
    }

    /* append terminating null byte */
    *buf_pos = '\0';

    /* dynamically allocate memory to store and return the final string */
    char* cigar_string = calloc( buf_pos - buf, sizeof(char) );
    strcpy( cigar_string, buf );

    return cigar_string;
}

int
get_start_for_sublocations_group(
        mapped_read_sublocation* group )
{
    /* simply return the start of the first sublocation in the group */
    return group->start_pos;
}

int
get_stop_for_sublocations_group(
        mapped_read_sublocation* group )
{
    /* loop through to the last sublocation, and return its
     * start + length */

    char* ptr = (char*) group;

    while( true )
    {
        if( end_of_sublocations_group( (mapped_read_sublocation*) ptr ) )
            break;

        ptr += sizeof( mapped_read_sublocation );
    }

    mapped_read_sublocation* final_sublocation
        = (mapped_read_sublocation*) ptr;

    return final_sublocation->start_pos + final_sublocation->length;
}

unsigned short
build_flag_for_sublocations_group(
        mapped_read_sublocation* group,
        struct read* read,
        int subtemplate_index )
{
    unsigned short flag = 0;

    /* HACK - hardcode for single and paired end for now */
    // TODO generalize this
    
    assert( read->num_subtemplates == 1 ||
            read->num_subtemplates == 2 );

    /* if the first subread was reverse complemented, set this flag */
    if( group->rev_comp ) {
        flag |= BAM_FREVERSE;
    }

    if( read->num_subtemplates == 2 )
    {
        assert( subtemplate_index == 0 ||
                subtemplate_index == 1 );

        flag |= BAM_FPAIRED;
        flag |= BAM_FPROPER_PAIR;

        if( subtemplate_index == 0 )
        {
            flag |= BAM_FREAD1;
        } else { // subtemplate_index == 1
            flag |= BAM_FREAD2;
        }
    }

    return flag;
}

char*
build_mate_info_for_sublocations_group(
        mapped_read_sublocation* group,
        struct read* read,
        int subtemplate_index,
        int pnext )
{
    /* HACK - hardcode for single and paired end for now */
    // TODO generalize this
    
    assert( read->num_subtemplates == 1 ||
            read->num_subtemplates == 2 );

    char buf[512];
    char* buf_pos = buf;

    if( read->num_subtemplates == 1 )
    {
        /* empty since this is not paired */
        buf_pos += sprintf( buf_pos, "*\t0\t0\t" );
    } else { // read->num_subtemplates == 2

        /* append the mate reference name */
        /* note - I think we can use "=" here (equivalent, shorter) */
        /* this is what the original code was doing */
        buf_pos += sprintf( buf_pos, "%s\t", read->name );

        /* append the mate start position */
        buf_pos += sprintf( buf_pos, "%u\t", pnext );

        /* append the inferred insert size */
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
        // TODO compute the inferred insert size properly
        buf_pos += sprintf( buf_pos, "%i\t", 0 );
    }

    /* append terminating null byte */
    *buf_pos = '\0';

    char* mate_info = calloc( buf_pos - buf, sizeof(char) );
    strcpy( mate_info, buf );

    return mate_info;
}

void
fprintf_sam_line_from_sublocations_group(
        FILE* sam_fp,
        mapped_read_location_prologue* location_prologue,
        mapped_read_sublocation* group,
        struct read* read,
        int subtemplate_index,
        struct genome_data* genome,
        float cond_prob )
{
    struct read_subtemplate* subtemplate
        = read->subtemplates + subtemplate_index;

    const int chr_index = location_prologue->chr;
    const int start = get_start_for_sublocations_group( group );
    const int stop = get_stop_for_sublocations_group( group );
    const float seq_error = location_prologue->seq_error;

    int rd_len = MAX( start, stop )
               - MIN( start, stop );

    /* print the query name */
    fprintf( sam_fp, "%s\t", read->name );

    /* build and print the flag */
    unsigned short flag = build_flag_for_sublocations_group(
            group, read, subtemplate_index );
    fprintf( sam_fp, "%hu\t", flag );

    /* print the chromosome name */
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
             (unsigned int) MIN( 254, (-10*log10( 1- cond_prob )) ) );

    /* print the cigar string */
    char* cigar_string = build_cigar_string_for_sublocations_group( group );
    fprintf( sam_fp, "%s\t", cigar_string );
    free( cigar_string );

    /* print the mate information */
    // TODO figure out pnext (not 0)
    char* mate_info = build_mate_info_for_sublocations_group(
            group, read, subtemplate_index, 0 );
    fprintf( sam_fp, "%s", mate_info );
    free( mate_info );

    /* print the actual sequence */
    fprintf( sam_fp, "%.*s\t", rd_len, subtemplate->char_seq );

    /* print the quality string */
    fprintf( sam_fp, "%.*s\t", rd_len, subtemplate->error_str );

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
    fprintf( sam_fp, "PQ:i:%u\t",
             (unsigned int) MIN(254, (-10*log10(1-seq_error))));
    fprintf( sam_fp, "XQ:i:%e\t", seq_error);
    fprintf( sam_fp, "XP:i:%e\n", cond_prob);

    return;
}

void
fprintf_mapped_read_to_sam( 
    FILE* sam_fp,
    mapped_read_t* mpd_rd,
    mapped_read_index* mpd_rd_index,
    struct cond_prbs_db_t* cond_prbs_db,    
    struct genome_data* genome,
    struct read* r )
{
    /* HACK - assumptions to get this to compile */
    assert( r->num_subtemplates == 1 || r->num_subtemplates == 2 );

    MPD_RD_ID_T i = 0; // TODO why MPD_RD_ID_T? it is confusing
    for( i = 0; i < mpd_rd_index->num_mappings; i++ )
    {
        mapped_read_location* mapping = mpd_rd_index->mappings[i];

        float cond_prob =
            get_cond_prb( cond_prbs_db, mpd_rd_index->read_id, i );

        char* ptr = (char*) mapping;

        mapped_read_location_prologue* prologue
            = (mapped_read_location_prologue*) ptr;

        /* move pointer to the start of the first mapped_read_sublocation */
        ptr += sizeof( mapped_read_location_prologue );
        enum bool are_more_sublocations = true;

        /* Find groups of mapped_read_sublocations to pass to
         * fprintf_sam_line. A group of sublocations is simply defined by
         * a pointer to the start. It ends when next_subread_is_ungapped,
         * indicating a gap in the template. */
        int current_read_subtemplate_index = 0;
        mapped_read_sublocation* current_subgroup
            = (mapped_read_sublocation*) ptr;

        while( are_more_sublocations )
        {
            assert( current_read_subtemplate_index < r->num_subtemplates );

            mapped_read_sublocation* current_sublocation
                = (mapped_read_sublocation*) ptr;

            if( end_of_sublocations_group( current_sublocation ) )
            {
                /* print out a SAM line for this group of sublocations */
                fprintf_sam_line_from_sublocations_group(
                        sam_fp,
                        prologue,
                        current_subgroup,
                        r,
                        current_read_subtemplate_index,
                        genome,
                        cond_prob
                    );

                /* if that was the last sublocation in this
                 * mapped_read_location, break */
                if( last_sublocation_in_mapped_read_location(
                            current_sublocation ) )
                {
                    are_more_sublocations = false;
                }

                /* update the current_subgroup pointer to point to the next
                 * sublocation, which is the start of the next subgroup */
                current_subgroup = (mapped_read_sublocation*)
                    (ptr + sizeof( mapped_read_sublocation ));
                current_read_subtemplate_index++;
            }

            ptr += sizeof( mapped_read_sublocation );
        }
    }

    return;
}

void
write_mapped_reads_to_sam( 
    struct rawread_db_t* rdb,
    struct mapped_reads_db* mappings_db,
    struct cond_prbs_db_t* cond_prbs_db,
    struct genome_data* genome,
    enum bool reset_cond_read_prbs,
    FILE* sam_ofp )
{
    int error;

    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );

    struct read* rd;
    mapped_read_t* mapped_rd;

    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    assert( error == 0 );

    get_next_read_from_rawread_db( 
        rdb, &rd, -1 );

    /* iterate through rawreads until we get to the rawread that matches the
     * first mapped read. This sets up the following while loop */
    while( NULL != rd 
           && NULL != mapped_rd
           && rd->read_id < get_read_id_from_mapped_read( mapped_rd ) )
    {
        free_read( rd );
        
        get_next_read_from_rawread_db( 
            rdb, &rd, -1 );
    }
    /* since the rawread db has every read, we should be guaranteed of having a 
       match at this point */
    
    
    while( rd != NULL
           && mapped_rd != NULL ) 
    {
        mapped_read_index* mapped_rd_index;
        init_mapped_read_index( &mapped_rd_index, mapped_rd );

        if( rd->read_id != mapped_rd_index->read_id )
        {
            fprintf( stderr, "read_id: %d\n", rd->read_id );
            fprintf( stderr, "mapped_rd->read_id: %d\n",
                     mapped_rd_index->read_id );
        }
        assert( rd->read_id == mapped_rd_index->read_id );
        
        if( rd->read_id > 0 && rd->read_id%1000000 == 0 )
            fprintf( stderr, "NOTICE       : Written %u reads to sam\n", rd->read_id );
        
        /* We test for mapped read NULL in case the last read was unmappable */
        if( mapped_rd != NULL 
            && mapped_rd_index->num_mappings > 0 )
        {
            /* sometimes we want the marginal distribution */
            if( reset_cond_read_prbs )
                reset_read_cond_probs( cond_prbs_db, mapped_rd, mappings_db );

            fprintf_mapped_read_to_sam( 
                sam_ofp, mapped_rd, mapped_rd_index,
                cond_prbs_db, genome, rd );
        }
        
        free_mapped_read_index( mapped_rd_index );
        
        /* Free the read */
        free_read( rd );

        error = get_next_read_from_mapped_reads_db( 
            mappings_db, 
            &mapped_rd
        );

        get_next_read_from_rawread_db( 
            rdb, &rd, -1 );
        
        while( NULL != rd 
               && NULL != mapped_rd
               && rd->read_id < get_read_id_from_mapped_read( mapped_rd ) )
        {
            free_read( rd );
            
            get_next_read_from_rawread_db( 
                rdb, &rd, -1 );
        }
        
    }
    
    goto cleanup;

cleanup:

    /* Free the read */
    if( rd != NULL )
        free_read( rd );
        
    return;
}

/*************************************************************************
 *
 *  Write nonmapping/unmappable reads to FASTQ
 * 
 */

/*
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
    
    /* Join all candidate mappings */
    /* get the cursor to iterate through the reads */
    rewind_rawread_db( rdb );
    rewind_mapped_reads_db( mappings_db );
    
    struct read* rd;
    mapped_read_t* mapped_rd;
    
    error = get_next_read_from_mapped_reads_db( 
        mappings_db, 
        &mapped_rd
    );
    assert( error == 0 );
    
    get_next_read_from_rawread_db( 
        rdb,  &rd, -1 );

    while( rd != NULL )
    {
        mapped_read_index* mapped_rd_index = NULL;

        /* if the mapped reads db is empty and there are still rawreads, these
         * rawreads were nonmapping */
        if( mapped_rd == NULL )
        {
            write_nonmapping_read_to_fastq( rd,
                    rdb->unmappable_single_end_reads,
                    rdb->unmappable_paired_end_1_reads,
                    rdb->unmappable_paired_end_2_reads );
        } else {
            init_mapped_read_index( &mapped_rd_index, mapped_rd );

            /* if this rawread doesn't have an associated mapped reads */
            if( rd->read_id < mapped_rd_index->read_id )
            {
                /* then it was unmappable, and was never added to the mpd rd db */
                write_nonmapping_read_to_fastq( rd,
                        rdb->unmappable_single_end_reads,
                        rdb->unmappable_paired_end_1_reads,
                        rdb->unmappable_paired_end_2_reads );
            }
            /* if the read has an associated mapped reads */
            else if( rd->read_id == mapped_rd_index->read_id &&
                     mapped_rd_index->num_mappings == 0 )
            {
                /*
                   If the read has an associated mapped read, but the mapped read
                   has no mappings, it was declared mappable but did not map.

                   Nonmapping reads are added to the mapped reads db; see
                   find_candidate_mappings.c:1028
                */
                write_nonmapping_read_to_fastq( rd,
                        rdb->non_mapping_single_end_reads,
                        rdb->non_mapping_paired_end_1_reads,
                        rdb->non_mapping_paired_end_2_reads );
            }
        }

        /* if we need to get the next mapped read */
        /* if mapped_rd is null, we are out of mapped reads */
        if( mapped_rd != NULL
            && rd->read_id == mapped_rd_index->read_id )
        {
            error = get_next_read_from_mapped_reads_db( 
                mappings_db, 
                &mapped_rd
            );
        }

        free_mapped_read_index( mapped_rd_index );

        /* Free the read */
        free_read( rd );
        
        /* we always get the next read */
        get_next_read_from_rawread_db( 
            rdb, &rd, -1 );

    }

    goto cleanup;

cleanup:
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


