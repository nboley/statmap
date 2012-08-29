/* Copyright (c) 2009-2012, Nathan Boley */

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
#include "genome.h"
#include "diploid_map_data.h"
#include "util.h"

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
        struct read_subtemplate* rst,
        float max_penalty_spread
    )
{
    /****** initialize the mapped_location info that we know  ******/
    /* copy the candidate map location template */
    candidate_mapping cand_map;
    memset( &cand_map, 0, sizeof(cand_map) );

    /* Set the read length */
    cand_map.rd_len = rst->length;

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
    assert( rst->pos_in_template.number_of_reads_in_template == 1 ||
            rst->pos_in_template.number_of_reads_in_template == 2 );
    /* The number of reads in the tempate tells us whether this read
     * subtemplate is from a single or paired end read */
    if( rst->pos_in_template.number_of_reads_in_template == 1 )
    {
        assert( rst->pos_in_template.pos == POS_SINGLE_END );
        cand_map.rd_type = SINGLE_END;
    }
    else if ( rst->pos_in_template.number_of_reads_in_template == 2 )
    {
        /* The position in the template tells us which end of the paired end
         * read this subtemplate represents */

        assert( rst->pos_in_template.pos == POS_PAIRED_END_1 ||
                rst->pos_in_template.pos == POS_PAIRED_END_2 );

        if( rst->pos_in_template.pos == POS_PAIRED_END_1 )
        {
            cand_map.rd_type = PAIRED_END_1;
        } else if ( rst->pos_in_template.pos == POS_PAIRED_END_2 )
        {
            cand_map.rd_type = PAIRED_END_2;
        }
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

void
sort_candidate_mappings( candidate_mappings* mappings )
{
    qsort( mappings->mappings, 
           mappings->length, 
           sizeof(candidate_mapping),
           (int(*)(const void*, const void*))cmp_candidate_mappings
    );
}

/* Append the candidate mappings from src onto dest */
void
append_candidate_mappings(
        candidate_mappings* dest,
        candidate_mappings* src
    )
{
    int i;
    for( i = 0; i < src->length; i++ )
    {
        add_candidate_mapping( dest, &(src->mappings[i]) );
    }
}

/*
 * END Candidate Mapping
 *********************************************************************************/
