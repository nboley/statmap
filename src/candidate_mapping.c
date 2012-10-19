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
    
    // Pseudo locations should already have been expanded
    if( chr == PSEUDO_LOC_CHR_INDEX ) {
        perror( "ERROR: Pseudo locs should NEVER be passed to modify_mapped_read_location_for_index_probe_offset. THIS SHOULD NEVER HAPPEN, PLEASE REPORT THIS BUG.\n" );
        return -1;
    }

    /* first deal with reads that map to the 5' genome */
    if( strnd == FWD )
    {
        /* if the mapping location of the probe is less than the length of the
           probe offset, then the actual read is mapping before the start of
           the genome, which is clearly impossible. We add the softclip length
           here under the assumption that any soft clipped basepairs may not
           have come from the genome (may have been part of a primer sequence,
           etc.)
        */
        if( read_location + softclip_len < subseq_offset ) 
        {
            return -1;
        } 
        /* we shift the location to the beggining of the sequence, 
           rather than the subseq that we looked at in the index  */
        else {
            read_location -= (subseq_offset - softclip_len);
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
            - ( subseq_len + subseq_offset - softclip_len )
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

int
get_length_from_cigar_string( candidate_mapping* mapping )
{
    int fragment_length = 0;

    int i;
    for( i = 0; i < mapping->cigar_len; i++ )
    {
        fragment_length += mapping->cigar[i].len;
    }

    return fragment_length;
}

/*********************************************************************************
 *
 * Candidate Mappings Code
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

void
copy_candidate_mappings(
        candidate_mappings* dst,
        candidate_mappings* src
    )
{
    /* assumes dst is initialized */
    assert( dst != NULL );

    dst->length = src->length;
    dst->allocated_length = src->allocated_length;

    /* Allocate new memory equivalent to the amount allocated in the original
     * struct. We need to realloc since dst has been initialized and already
     * has memory allocated for mappings. */
    dst->mappings = realloc( dst->mappings,
            dst->allocated_length * sizeof(candidate_mapping) );
    assert( dst->mappings != NULL );

    /* Copy the candidate mappings */
    memcpy( dst->mappings, src->mappings, src->length * sizeof(candidate_mapping) );

    return;
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
print_candidate_mapping_cigar_string( candidate_mapping* mapping )
{
    printf("Cigar:        ");
    int i;
    for( i = 0; i < mapping->cigar_len; i++ )
    {
        printf("%c%i", mapping->cigar[i].op, mapping->cigar[i].len );
    }
    printf("\n");
}

void
print_candidate_mapping( candidate_mapping* mapping )
{
    printf("Chr:          %u\n", mapping->chr);
    printf("Start BP:     %u\n", mapping->start_bp);
    printf("Read Len:     %u\n", mapping->mapped_length);
    printf("Read Strand:  %u\n", mapping->rd_strnd);
    printf("Penalty:      %.2f\n", mapping->penalty);

    print_candidate_mapping_cigar_string( mapping );

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

    /* first, sort by subtemplate index (single or paired end) */
    if( m1->pos_in_template != m2->pos_in_template )
        return m1->pos_in_template - m2->pos_in_template;

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
