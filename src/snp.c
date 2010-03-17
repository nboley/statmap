#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "genome.h"
#include "quality.h"
#include "index_genome.h"

typedef struct {
    char alt_bp;
    float alt_frac;
    GENOME_LOC_TYPE loc;    
} snp_t;

/********************************************************************************
 * Methods for dealing with the snp type
 *
 */

int 
cmp_snp( const void* void_snp1, const void* void_snp2 )
{
    const snp_t* snp1 = void_snp1;
    const snp_t* snp2 = void_snp2;

    int cmp = snp1->loc.chr - snp2->loc.chr;;

    if( cmp == 0 )
        cmp = snp1->loc.loc - snp2->loc.loc;

    return cmp;
}

int 
cmp_snp_to_location( snp_t* snp, 
                     int chr,
                     unsigned int start )
{
    int cmp = snp->loc.chr - chr;

    if( cmp == 0 )
        cmp = snp->loc.loc - start;

    return cmp;
}


void
fprintf_snp( FILE* fp, snp_t* snp )
{
    fprintf( fp, "Chr: %u\tLoc: %u\tAlt Allele: %c\tProb: %e\n",
             snp->loc.chr, snp->loc.loc, snp->alt_bp, snp->alt_frac );
    return;
}


/********************************************************************************
 * Methods for dealing with the snp database ( for now just a sorted array )
 *
 */

struct snp_db_t {
    int num_snps;
    snp_t* snps;    
};


void
init_snp_db( struct snp_db_t** snp_db  )
{
    *snp_db = malloc( sizeof( struct snp_db_t ) );
    if( *snp_db == NULL )
    {
        fprintf( stderr, "FATAL:     memory allocation failed in init_snp_db." );
        assert( false );
        exit( -1 );
    }

    (*snp_db)->num_snps = 0;
    (*snp_db)->snps = NULL;
    
    return;
}

void
free_snp_db( struct snp_db_t* snp_db )
{
    free( snp_db->snps );
    free( snp_db );
}

void
sort_snp_db( struct snp_db_t* snp_db )
{
    qsort( snp_db->snps, snp_db->num_snps, sizeof(snp_t), cmp_snp );
    return;
}

void
find_snps_in_snp_db( struct snp_db_t* snp_db, int chr, int start, 
                     int stop, /* NOT including stop */ 
                     int* num_snps, int** snp_indexes )
{
    *num_snps = 0;
    *snp_indexes = malloc( (stop-start)*sizeof(int) );

    /* we use a binary search to find the insert location for start */
    /* then, we scan along until we have passed stop */
    
    /* 
     * find the correct insert location - we use a binary sort to keep 
     * the items in order - we dont use the stdlib because it doesnt find
     * the insert location
     */
    int low = 0;
    int high = snp_db->num_snps;
    snp_t* snps = snp_db->snps;

    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        int comparison = cmp_snp_to_location( snps + mid, chr, start );
        if( comparison < 0 ) {
            low = mid + 1; 
        } else {
            /* can't be high = mid-1: here A[mid] >= value, */
            /* so high can't be < mid if A[mid] == value    */
            high = mid; 
        }
    }

    assert( snps[low].loc.chr == chr );
    assert( snps[low].loc.loc >= start );
    /* TODO - optimize the start */
    while( low < snp_db->num_snps
           && snps[low].loc.chr == chr
           && snps[low].loc.loc < stop )
    {
        (*snp_indexes)[*num_snps] = low;
        *num_snps += 1;
        low += 1;
    }
}

void
fprintf_snp_db( FILE* fp, struct snp_db_t* snp_db )
{
    int i;
    for( i = 0; i < snp_db->num_snps; i++ )
    {
        fprintf_snp( fp, snp_db->snps + i );
    }
    return;
}

void
add_snp_to_snp_db( struct snp_db_t* snp_db, snp_t* snp )
{
    snp_db->snps = realloc(snp_db->snps, sizeof(snp_t)*(snp_db->num_snps+1));
    if( snp_db->snps == NULL )
    {
        fprintf( stderr, "FATAL:     memory allocation failed in add_snp_to_snp_db." );
        assert( false );
        exit( -1 );
    }

    snp_db->snps[snp_db->num_snps] = *snp;
    snp_db->num_snps += 1;
    
    return;
}

/*
 * Parse a snp cov file and add all of the snps.
 *
 * chr pos snpid A C G T N refallele allele1 allele2 state1 state2 count1 count2 assay cell group
 * chr1 981801 tmp:1:981801 0 11 0 2 0 C C T H H 11 2 Pol2HudsonAlpha GM12878 BU
 * chr1 1479937 tmp:1:1479937 7 0 0 0 0 A A G H H 7 0 DNase GM12878 BU
 * chr1 1600387 tmp:1:1600387 0 1 2 10 0 T C T H H 1 10 CTCF GM12878 BU
 *
 */

void
build_snp_db_from_snpcov_file( FILE* fp, genome_data* genome )
{
    int line_num = 0;
    
    char* fgets_ret = NULL;
    int fscanf_ret = 0;
 
    /* initialize the snp db */
    init_snp_db( &(genome->snp_db) );
   
    /* make sure that we are at the begining of the file */
    fseek( fp, 0, SEEK_SET );
    
    /* store temporary lines */
    char buffer[1000];
    memset( buffer, 0, 1000*sizeof(char) );

    /* check for a header */
    fgets_ret = fgets( buffer, 999, fp );
    /* If the line is NULL or doesnt contain a new line, exit */
    if( NULL == fgets_ret || NULL == strrchr( buffer, '\n' ) )
    {
        fprintf(stderr, "FATAL:     Could not read a full line from the snp file.\n");
        perror( "Error is");
        assert( false );
        exit( -1 );
    }

    /* Search for snpid. If dont we find it, assume there is no header */
    if( NULL == strstr( buffer, "snpid" ) )
    {
        fprintf(stderr, "NOTICE      :  No snpcov file header.\n" );
        fseek( fp, 0, SEEK_SET );
    } else {
        line_num += 1;
    }

    char chr[255];
    unsigned int pos;
    char refallele;
    char altallele;
    char allele1;
    char allele2;
    
    int cnt1;
    int cnt2;
    float alt_frac_est;
    
    snp_t snp;
    
    /* Assume that we are at the first worthwhile line */
    while( !feof(fp) )
    {
        line_num += 1;
        
        fscanf_ret = fscanf( 
            fp, 
            "%s\t%u\t%*s\t%*i\t%*i\t%*i\t%*i\t%*i\t%c\t%c\t%c\t%*c\t%*c\t%i\t%i\t",
            chr, &pos, &refallele, &allele1, &allele2, &cnt1, &cnt2
        );
        
        fgets_ret = fgets( buffer, 255, fp );
        
        /* If we dont read everything that we expect */
        if( NULL == fgets_ret )
        {
            /* If we are at an end of file, then we are done reading */
            if( feof( fp ) ) {
                break;
            } 
        }

        /* FIXME - figure out what this means */
        /* 
         * If allele1 == allele2, we think this means that 
         * both the mother and father are the same ( ie, its just an
         * error in the ref genome ) so skip it.
         */
        if( allele1 == allele2 )
            continue;
 
        /* determine the alternate allele */
        if( refallele == allele1 ) {
            altallele = allele2;
            alt_frac_est = ( (float) cnt2 ) / ( cnt1 + cnt2  );
        } else {
            if( altallele != allele2 ) {
                fprintf(stderr, "WARNING     : ref allele is not the same as allele1 *or* allele 2. Skipping line %i\n", line_num);
                continue;
            }
            
            altallele = allele1;
            alt_frac_est = ( (float) cnt1 ) / ( cnt1 + cnt2  );
        }
        
        /*** add this allele to the snp count */
        /* zero out the data */
        memset( &snp, 0, sizeof(snp_t) );

        /** first, edit the genome location */
        /* first, find the chr index */
        int chr_index = find_chr_index( genome, chr );
        /* If we couldnt find this chr name */
        if( -1 == chr_index )
        {
            fprintf(stderr, "WARNING     : couldn't find chr '%s'. Skipping line %i\n", chr, line_num);
            continue;
        } else {
            assert( chr_index >= 0 && chr_index <= CHR_NUM_MAX );
            snp.loc.chr = chr_index;
        }
        /* add the location */
        assert( pos <= LOCATION_MAX );
        snp.loc.loc = pos;
        /* note that the other fields are zero from the memset */

        /* set the alt allele */
        snp.alt_bp = altallele;
        snp.alt_frac = alt_frac_est;

        add_snp_to_snp_db( genome->snp_db, &snp );
    }

    /* sort the snp db */
    sort_snp_db( genome->snp_db );

    fprintf( stderr, "NOTICE      :  Loaded %i/%i snps\n", 
             genome->snp_db->num_snps, line_num-1 );
    
    return;
}

void
index_snp_sequences( struct snp_db_t* snp_db, genome_data* genome )
{
    /* initialize the constant loc elements */
    GENOME_LOC_TYPE loc;
    loc.snp_coverage = 0;

    /* Not a junction read */
    loc.read_type = 0;

    int seq_len = genome->index->seq_length;
    char* tmp_seq = malloc(seq_len*sizeof(char));

    int i;
    unsigned int j;
    /*
     * Iterate through each snp in the list. This snp 
     * is *always* set - other snps that overlap the same
     * ref sequence are handled inside the loop, and may or
     * may not be set.
     *
     */
    for( i = 0; i < snp_db->num_snps; i++ )
    {
        snp_t* snp = snp_db->snps + i;
        
        for( j = MAX(0, snp->loc.loc - seq_len + 1); 
             j <= MIN( snp->loc.loc, genome->chr_lens[snp->loc.chr]-seq_len );
             j++
            )
        {
            memcpy( tmp_seq, genome->chrs[snp->loc.chr] + j, 
                    sizeof(char)*seq_len );

            /* find out how many snp's are upstream of this one */
            /* that is, in the interval that contains the i'th snp
               and is less than  snp->loc.loc + j + seq_len 
               OBVIOUSLY FIX ME 
            */

            /* add in the alternate allele */
            tmp_seq[ snp->loc.loc - j  ] = snp->alt_bp;

            LETTER_TYPE *translation;
            translate_seq( tmp_seq, seq_len, &(translation));
            
            /* if we cant add this sequence ( probably an N ) continue */
            if( translation == NULL ) {
                continue;
            }
            
            loc.chr = snp->loc.chr;
            loc.loc = j;
            
            /* Add the sequence into the tree */
            add_sequence(genome->index, translation, seq_len, loc);
            
            free( translation );
        }

    }

    free( tmp_seq );

    return;
}

