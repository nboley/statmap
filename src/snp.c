#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "snp.h"
#include "mapped_read.h"
#include "quality.h"
#include "index_genome.h"

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
    fprintf( fp, "Chr: %u\tLoc: %u\tAlt Allele: %c\tProb: %e\tRef Cnt: %.2f\tAlt Cnt: %.2f\n",
             snp->loc.chr, snp->loc.loc, snp->alt_bp, snp->alt_frac, snp->ref_cnt, snp->alt_cnt );
    return;
}


/********************************************************************************
 * Methods for dealing with the snp database ( for now just a sorted array )
 *
 */

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

/*
 * Return the snps that cover basepairs in [start, stop)
 *
 */
void
find_snps_in_snp_db( struct snp_db_t* snp_db, 
                     int chr, int start, 
                     int stop, /* NOT including stop */ 
                     int* num_snps, snp_t** snps )
{
    *num_snps = 0;
    
    /* we use a binary search to find the insert location for start */
    /* then, we scan along until we have passed stop */
    
    /* 
     * find the correct insert location - we use a binary sort to keep 
     * the items in order - we dont use the stdlib because it doesnt find
     * the insert location
     */
    int low = 0;
    int high = snp_db->num_snps;
    snp_t* all_snps = snp_db->snps;

    while (low < high) {
        int mid = low + ((high - low) / 2) ;
        int comparison = cmp_snp_to_location( all_snps + mid, chr, start );
        if( comparison < 0 ) {
            low = mid + 1; 
        } else {
            /* can't be high = mid-1: here A[mid] >= value, */
            /* so high can't be < mid if A[mid] == value    */
            high = mid; 
        }
    }

    assert( all_snps[low].loc.chr == chr );
    assert( all_snps[low].loc.loc >= start );
    *snps = all_snps + low;
    // fprintf(stderr, "Low: %i\t Start: %i\n", low, start);
    /* TODO - optimize the start */
    while( low < snp_db->num_snps
           && all_snps[low].loc.chr == chr
           && all_snps[low].loc.loc < stop )
    {
        *num_snps += 1;
        low += 1;
    }

    return;
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
build_snp_db_from_snpcov_file( FILE* fp, struct genome_data* genome )
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
    char altallele = '\0';
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
        
        /* set the new items */
        snp.alt_cnt = 0;
        snp.ref_cnt = 0;

        add_snp_to_snp_db( genome->snp_db, &snp );
    }

    /* sort the snp db */
    sort_snp_db( genome->snp_db );

    fprintf( stderr, "NOTICE      :  Loaded %i/%i snps\n", 
             genome->snp_db->num_snps, line_num-1 );
    
    return;
}

void
update_snp_estimates_from_candidate_mappings( 
    struct mapped_reads_db* rdb, struct genome_data* genome )
{
    /* zero out the snp counts */
    int snp_index;
    for( snp_index = 0; snp_index < genome->snp_db->num_snps; snp_index++ )
    {
        (genome->snp_db->snps)[snp_index].alt_cnt = 0;
        (genome->snp_db->snps)[snp_index].ref_cnt = 0;
    }

    /* update the counts */
    struct mapped_read_t* rd;    
    rewind_mapped_reads_db( rdb );
    while( EOF != get_next_read_from_mapped_reads_db( rdb, &rd ) )
    {
        int i;
        for( i = 0; i < rd->num_mappings; i++ )
        {
            struct mapped_read_location* loc = rd->locations + i;
            if( ( get_flag_from_mapped_read_location( loc ) 
                  & FIRST_READ_COVERS_SNP ) > 0 )
            {
                int num_snps;
                snp_t* snps;
                
                find_snps_in_snp_db( 
                    genome->snp_db, 
                    get_chr_from_mapped_read_location( rd->locations + i ),
                    get_start_from_mapped_read_location( rd->locations + i ),
                    get_stop_from_mapped_read_location( rd->locations + i ),
                    &num_snps,
                    &snps
                );

                /* TODO optimize this to use a bitshift at each loop */
                for( snp_index = 0; snp_index < num_snps; snp_index++ )
                {                        
                    /* If the correct bit is set */
                    /* We cap the min penalties to enable the fast math optimizations */
                    if( ( (loc->snps_bm_r1)&(1<<snp_index) ) > 0 ) 
                    {
                        snps[snp_index].alt_cnt += 
                            get_cond_prob_from_mapped_read_location(
                                rd->locations + i);
                    } else {
                        snps[snp_index].ref_cnt += 
                            get_cond_prob_from_mapped_read_location(
                                rd->locations + i);
                    }
                }                    
            }
        }        
    }
}


void
write_snps_to_file( FILE* ofp, struct genome_data* genome )
{
    /* Write the header */
    fprintf( ofp, "chr_name\tposition\t \tRef_Allelle\tAllele1\tAllele2\t \tCount1\tCount2\n" );

    /* Loop through each snp in the db, and print them */
    int snp_index;
    for( snp_index = 0; snp_index < genome->snp_db->num_snps; snp_index++ )
    {
        /* Store a local ref - this should be optimized out */
        snp_t* snp;
        snp = genome->snp_db->snps + snp_index;
        
        /* determine the reference basepair */
        char refallele = toupper( genome->chrs[ snp->loc.chr ][ snp->loc.loc ] );
        
        fprintf( 
            ofp, 
            "%s\t%u\t \t%c\t%c\t%c\t \t%.0f\t%.0f\n",
            genome->chr_names[ snp->loc.chr ], snp->loc.loc, 
            refallele, refallele, snp->alt_bp, 
            snp->ref_cnt, snp->alt_cnt
        );
    }
    
}

/* returns size written */
size_t
write_snp_db_to_binary_file( struct snp_db_t* snp_db, FILE* ofp )
{
    int rv;
    
    size_t size_written = 0;
    
    /* if the snp db is empty, there are just 0 snps */
    if( snp_db == NULL )
    {
        int zero = 0;
        
        rv = fwrite( &zero, sizeof(int), 1, ofp );
        assert( rv == 1 );
        size_written += sizeof(int);
        
        rv = fwrite( &zero, sizeof(size_t), 1, ofp );
        assert( rv == 1 );
        size_written += sizeof(size_t);
        
        return size_written;
    }
    
    /* write the number of snps */
    rv = fwrite( &(snp_db->num_snps), sizeof(int), 1, ofp ); 
    assert( rv == 1 );
    size_written += sizeof(int);
    
    /* this is redundant - but it provides a crude architecture check */
    /* write the size of the snp array */
    size_t size = (snp_db->num_snps)*sizeof(snp_t);
    rv = fwrite( &size, sizeof(size_t), 1, ofp );     
    assert( rv == 1 );
    size_written += sizeof(size_t);

    /* write the snp array */
    rv = fwrite( snp_db->snps, sizeof(snp_t), snp_db->num_snps, ofp ); 
    assert( rv == snp_db->num_snps );
    size_written += sizeof(snp_t)*snp_db->num_snps;
    
    return size_written;
}

void
load_snp_db_from_mmapped_data( struct genome_data* genome, char* data )
{
    int num_snps = -1;
    size_t size = 0;
    
    num_snps = *( (int*) data );
    data += sizeof( int );
    
    size = *( (size_t*) data );
    data += sizeof( size_t  );
    
    if( num_snps == 0 )
    {
        // BUG BUG
        // assert( 0 == size );
        genome->snp_db = NULL;
        return;
    } 

    init_snp_db( &(genome->snp_db) );
    
    genome->snp_db->num_snps = num_snps;
    
    if( size != num_snps*sizeof( snp_t )  )
    {
        fprintf( stderr, "FATAL       :  Error loading snp db from file - there appears to be a size mismatch.\n" );
        exit( 1 );
    }

    genome->snp_db->snps = ( snp_t* ) data;
    
    return;
}



