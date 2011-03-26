/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include "math.h"

#include "statmap.h"
#include "find_candidate_mappings.h"
#include "quality.h"
#include "index_genome.h"
#include "snp.h"
#include "mapped_location.h"

int
find_optimal_subseq_offset( 
    struct rawread* r,
    /* store the desired subsequences length */
    int subseq_len
) {
    assert( subseq_len <= r->length );
    
    /* XXX for now, we just use the first subseq_len characters,
       so the offset is always 0 */
    return 0;
};

void
search_index( struct index_t* index, 
              
              float min_match_penalty,
              float max_penalty_spread,
              mapped_locations** results,

              struct rawread* r,
              float* bp_mut_rates,

              float* lookuptable_position,
              float* inverse_lookuptable_position,
              float* reverse_lookuptable_position,
              float* reverse_inverse_lookuptable_position
    )
{
    /**** Prepare the read for the index search */
    /* 
       first, we need to find the subseq of the read that we want to 
       probe the index for. This is controlled by the index->seq_length.
    */
    int subseq_offset = find_optimal_subseq_offset( r, index->seq_length );
    int subseq_length = index->seq_length;
    
    /* Store a copy of the read */
    /* This read has N's replaced with A's, and might be RC'd */
    char* sub_read = calloc(subseq_length + 1, sizeof(char));
    assert( sub_read != NULL );
    /* note that the NULL ending is pre-set from the calloc */
    memcpy( sub_read, r->char_seq + subseq_offset, sizeof(char)*(subseq_length) );
    
    /** Deal with the read on the fwd strand */
    /* Store the translated sequences here */
    LETTER_TYPE *fwd_seq;
    fwd_seq = translate_seq( sub_read, subseq_length, &fwd_seq );
    /* If we couldnt translate it */
    if( fwd_seq == NULL )
    {
        // fprintf(stderr, "Could Not Translate: %s\n", r->char_seq);
        return;
    }
    assert( fwd_seq != NULL );
    
    /** Deal with the read on the opposite strand */
    LETTER_TYPE *bkwd_seq;
    char* tmp_read = calloc(subseq_length + 1, sizeof(char));
    rev_complement_read( sub_read, tmp_read, subseq_length );
    bkwd_seq = translate_seq( tmp_read, subseq_length, &bkwd_seq );
    // BUG why did I have this here?
    //replace_ns_inplace( tmp_read, r->length );
    assert( bkwd_seq != NULL );

    init_mapped_locations( results );
    (*results)->subseq_len = subseq_length;
    (*results)->subseq_offset = subseq_offset;
    
    /* map the full read */
    find_matches_from_root( index, 
                            
                            min_match_penalty,
                            max_penalty_spread,
                            *results,

                            /* length of the reads */
                            subseq_length,
                            
                            /* the fwd stranded sequence */
                            fwd_seq, 
                            lookuptable_position,
                            inverse_lookuptable_position,
                            
                            /* the bkwd stranded sequence */
                            bkwd_seq, 
                            reverse_lookuptable_position,
                            reverse_inverse_lookuptable_position,
                            
                            bp_mut_rates
        );
    
    /* Free the allocated memory */
    free( fwd_seq );
    free( bkwd_seq );

    free( sub_read );
    free( tmp_read );

    return;
};


static inline void 
recheck_location( struct genome_data* genome, 
                  struct rawread* r, candidate_mapping* loc,
                  const float* const lookuptable_position,
                  const float* const inverse_lookuptable_position,
                  const float* const reverse_lookuptable_position,
                  const float* const reverse_inverse_lookuptable_position,
                  const float* const bp_mut_rates
)
{
    const float* correct_lookuptable_position = NULL;
    const float* correct_inverse_lookuptable_position = NULL;
    
    /* find a pointer to the sequence at this genomic location */
    
    char* genome_seq = find_seq_ptr( 
        genome, 
        loc->chr, 
        loc->start_bp,
        r->length
    );
                    
    /* if the genome_seq pointer is null, the sequence isn't valid
       for some reason ( ie, it falls off the end of the chromosome). 
       In such cases, mark the location as invalid and continue */
    if( NULL == genome_seq )
    {
        loc->recheck = INVALID;
        return;
    }
                    
    char* mut_genome_seq = NULL;
    mut_genome_seq = malloc(sizeof(char)*(r->length+1));
    assert( mut_genome_seq != NULL ); 
                    
    if( BKWD == loc->rd_strnd )
    {
        rev_complement_read( genome_seq, mut_genome_seq, r->length );
        correct_lookuptable_position = reverse_lookuptable_position;
        correct_inverse_lookuptable_position = reverse_inverse_lookuptable_position;
    } else {
        memcpy( mut_genome_seq, genome_seq, sizeof(char)*r->length );
        mut_genome_seq[r->length] = '\0';
        correct_lookuptable_position = lookuptable_position;
        correct_inverse_lookuptable_position = inverse_lookuptable_position;
    }
    
    float rechecked_penalty = recheck_penalty( 
        mut_genome_seq, 
        // char* observed,
        r->char_seq,
        // const int seq_length,
        r->length,
        correct_lookuptable_position,
        correct_inverse_lookuptable_position,
        bp_mut_rates
    );

    /*
     * DEBUG
    printf( "%.2f\t%.2f\t%i\n", loc->penalty, rechecked_penalty, (loc->rd_strnd == BKWD) );
    printf( "%.*s\t%.*s\t%.*s\n", 
            r->length, r->char_seq, 
            r->length, genome_seq, 
            r->length, mut_genome_seq );
     */
    
    loc->penalty = rechecked_penalty;
    
    free( mut_genome_seq );
}


/* TODO - revisit the read length vs seq length distinction */
/* 
 * I use a struct for the parameters so that I can initialzie threads
 * with this function. 
 *
 */
void*
find_candidate_mappings( void* params )    
{
    /* 
     * recreate the struct parameters for readability
     * this should be optimized out 
     *
     */

    struct single_map_thread_data* td = params;

    int thread_id = td->thread_id;

    struct genome_data* genome = td->genome;
    
    unsigned int* mapped_cnt = td->mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex = td->mapped_cnt_mutex;

    struct rawread_db_t* rdb = td->rdb;
    
    pthread_mutex_t* mappings_db_mutex = td->mappings_db_mutex;
    candidate_mappings_db* mappings_db = td->mappings_db;
    /* The minimum absolute penalty
       that a valid read can have */
    float min_match_penalty = td->min_match_penalty;
    /* The minimum difference between the lowest penalty read and a
       valid read, set <= -1 to disable */
    float max_penalty_spread = td->max_penalty_spread;

    /* END parameter 'recreation' */

    assert( genome->index != NULL );

    clock_t start;
    start = clock();
    
    /* build the basepair mutation rate lookup table */
    float* bp_mut_rates;
    determine_bp_mut_rates( &bp_mut_rates );

    /* The current read of interest */
    long readkey;
    struct rawread *r1, *r2;
    /* 
     * While there are still mappable reads in the read DB. All locking is done
     * in the get next read functions 
     */
    while( EOF != 
           get_next_mappable_read_from_rawread_db( 
               rdb, &readkey, &r1, &r2 )  
         ) 
    {                
        /* We dont lock mapped_cnt because it's read only and we dont 
           really care if it's wrong 
         */
        if( readkey % 100000 == 0 && readkey > 0 )
        {
            fprintf(stderr, "DEBUG       :  Mapped %li reads, %i successfully\n", 
                    readkey, *mapped_cnt);
        }
        
        /* consider both read pairs */
        int j = 0;
        struct rawread* reads[2] = { r1, r2 };
        for( j = 0; j < 2 && reads[j] != NULL; j++ )
        {
            struct rawread* r = reads[j];

            /* Build the quality lookup tables */
            float* lookuptable_position = malloc(sizeof(float)*r->length);
            float* inverse_lookuptable_position = malloc(sizeof(float)*r->length);
            float* reverse_lookuptable_position = malloc(sizeof(float)*r->length);
            float* reverse_inverse_lookuptable_position = malloc(sizeof(float)*r->length);
            build_lookup_table_from_rawread(
                r, 
                lookuptable_position, inverse_lookuptable_position,
                reverse_lookuptable_position, reverse_inverse_lookuptable_position
            );

            /**** go to the index for mapping locations */
            mapped_locations *results = NULL;

            search_index( genome->index, 

                          min_match_penalty,
                          max_penalty_spread,

                          &results,

                          r,
                          bp_mut_rates,

                          lookuptable_position,
                          inverse_lookuptable_position,
                          reverse_lookuptable_position,
                          reverse_inverse_lookuptable_position
                );

            int subseq_len = results->subseq_len;
            int subseq_offset = results->subseq_offset;
            
            /****** Prepare the template candidate_mapping objects ***********/
            candidate_mapping template_candidate_mapping 
                = init_candidate_mapping_from_template( 
                    r, max_penalty_spread 
            );
            assert( template_candidate_mapping.rd_type != 0 );

            /**** TODO - get rid of this requirement */
            /* Assert that subseqs are false */
            assert( subseq_len <= r->length );
            
            /***** COPY information from the index lookup into the result set
             * build and populate an array of candidate_mapping's. 
             */        
            candidate_mappings* mappings;
            init_candidate_mappings( &mappings );

            /* 
             * We keep track of the max observed penalty so that we can filter
             * out penalties that are too low. The index will never return 
             * results that are strictly below the minimum penalty, but it may 
             * return  results below the relative penalty. ( Read the indexing
             *  header for details )
             */
            float max_penalty = min_match_penalty;

            unsigned int i;
            for( i = 0; i < results->length; i++ )
            {
                /* 
                 * I'm scribbling on the base relation, but it doesnt 
                 * matter because the add will copy it and then I will
                 * overwrite what I scribbled on anyways. 
                 */
                /* hopefully this will be optimized out */
                mapped_location* result;
                result = results->locations + i;

                /* set the strand */
                if( result->strnd == FWD )
                {
                    template_candidate_mapping.rd_strnd = FWD;
                } else {
                    assert( result->strnd == BKWD );
                    template_candidate_mapping.rd_strnd = BKWD;
                }

                /* set the chr */
                template_candidate_mapping.chr = (result->location).chr;

                /* set the location */
                int read_location = (result->location).loc;
                if( (result->location).chr != PSEUDO_LOC_CHR_INDEX ) {
                    read_location -= subseq_offset;
                    
                    if( result->strnd == BKWD ) {
                        read_location -= ( r->length - subseq_len );
                    }
                } 
                
                template_candidate_mapping.start_bp = read_location;
                
                /* set the penalty */
                template_candidate_mapping.penalty = result->penalty;

                /* if necessary, update the max observed penalty */
                if( result->penalty > max_penalty )
                    max_penalty = result->penalty;

                if( result->location.covers_snp == 1 )
                {
                    template_candidate_mapping.recheck = COVERS_SNP;
                    template_candidate_mapping.snp_bitfield = 
                        result->location.snp_coverage;
                    template_candidate_mapping.does_cover_snp = true;
                } else {
                    template_candidate_mapping.does_cover_snp = false;
                    assert( 0 == result->location.snp_coverage );
                    template_candidate_mapping.snp_bitfield = 0;
                }

                add_candidate_mapping( mappings, &template_candidate_mapping );
            }

            // print_mapped_locations( results );
            // print_candidate_mappings( mappings );
            // printf( "==========================================\n\n" );

            free_mapped_locations( results );

            /****** Do the recheck ******/
            /* 
             * Currently, everything should be set *except* the gene strand. 
             * This is because we dont know what it is. Therefore, we will 
             * add these to the db with that bit unset, and then during the 
             * merging stage add to the penalty and say that it was equally 
             * likely to have come from either gene strand. This corresponds
             * with us being equally certain that the read is from either gene.
             *
             * Also, we need to check that we dont have any low quality reads.
             * ( This is possible if the path search went awry, and we found 
             *   a low quality read before a high quality read )
             */


            /* first, if necessary, recheck the locations */
            int j;
            for( j = 0; j < mappings->length; j++ )
            {
                /*  if the location doesnt need to be rechecked, then continue */
                if( true ) // OR probe length is less than index length
                {
                    recheck_location( genome, r, mappings->mappings + j,
                                      lookuptable_position,
                                      inverse_lookuptable_position,
                                      reverse_lookuptable_position,
                                      reverse_inverse_lookuptable_position,
                                      bp_mut_rates
                        );

                    /* we may need to update the max penalty */
                    if( (mappings->mappings + j)->penalty > max_penalty ) {
                        max_penalty = (mappings->mappings + j)->penalty;
                    }
                }
            }

            for( j = 0; j < mappings->length; j++ )
            {
                /* this should be optimized out */
                candidate_mapping* loc = mappings->mappings + j;

                /* 
                 * We always need to do this because of the way the search queue
                 * handles poor branches. If our brnach prediction fails we could
                 * add a low probability read, and then go back and find a very 
                 * good read which would invalidate the previous. Then, the read
                 * may not belong, but it isnt removed in the index searching method.
                 */
                /* I set the safe bit to 1e-6, which is correct for a float */
                                
                /* we check max_penalty_spread > -0.00001 to make sure that the
                   likelihood ratio threshold is actually active */
                if(  max_penalty_spread > -0.00001  ) 
                {
                    if(   loc->penalty < ( max_penalty - max_penalty_spread ) )
                    {
                        loc->recheck = INVALID;
                    } 
                }
                
                /* make sure that the penalty isn't too low */
                if( loc->penalty < min_match_penalty  )
                {
                    loc->recheck = INVALID;
                } 

                /* if it's passed all of the tests, it must be valid */
                if( loc->recheck != INVALID )
                    loc->recheck = VALID;                
            }

            /* increment the number of reads that mapped, if any pass the rechecks */
            for( j = 0; j < mappings->length; j++ ) {
                if( (mappings->mappings + j)->recheck == VALID )
                {
                    pthread_mutex_lock( mapped_cnt_mutex );
                    *mapped_cnt += 1;
                    pthread_mutex_unlock( mapped_cnt_mutex );
                    break;
                }
            }
            
            /****** add the results to the database ******/
            /* BUG - assert the readnames are identical */
            /* BUG - why is this mutex necessary? */
            pthread_mutex_lock( mappings_db_mutex );
            assert( thread_id < num_threads );
            /* note that we add to the DB even if there are 0 that map,
               we do this so that we can join with the rawreads easier */
            add_candidate_mappings_to_db( 
                mappings_db, mappings, readkey, thread_id );
            pthread_mutex_unlock( mappings_db_mutex );

            free_candidate_mappings( mappings );

            /* cleanup the unmarshalled reads */
            free_rawread( r );

            /* free the mutation lookup tables */
            free( lookuptable_position );
            free( inverse_lookuptable_position );
            free( reverse_lookuptable_position );
            free( reverse_inverse_lookuptable_position );
        }
    }

    /* cleanup the bp mutation rates */
    free( bp_mut_rates );
    
    return NULL;
}

void
find_all_candidate_mappings( struct genome_data* genome,
                             FILE* log_fp,
                             struct rawread_db_t* rdb,

                             candidate_mappings_db* mappings_db,
                             float min_match_penalty,
                             float max_penalty_spread,
                             float max_seq_length

    )
{
    clock_t start = clock();

    /* 
     * init mutexes to guard access to the log_fp, the fastq 'db', and
     * the mappings db. They will be stored in the same structure as the 
     * passed data
     */

    /* initialize the necessary mutex's */
    pthread_mutex_t log_fp_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mapped_cnt_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t mappings_db_mutex = PTHREAD_MUTEX_INITIALIZER;

    /* put the search arguments into a structure */
    struct single_map_thread_data td_template;
    td_template.genome = genome;
    
    td_template.log_fp = log_fp;
    td_template.log_fp_mutex = &log_fp_mutex;

    unsigned int mapped_cnt = 0;
    td_template.mapped_cnt = &mapped_cnt;
    td_template.mapped_cnt_mutex = &mapped_cnt_mutex;

    td_template.rdb = rdb;
    
    td_template.mappings_db = mappings_db;
    td_template.mappings_db_mutex = &mappings_db_mutex;
    
    td_template.min_match_penalty = min_match_penalty;
    td_template.max_penalty_spread = max_penalty_spread;
    td_template.max_subseq_len = max_seq_length;

    /* initialize the threads */
    
    long t;
    int rc;
    void* status;
    pthread_t thread[num_threads];

    pthread_attr_t attrs[num_threads];
    
    struct single_map_thread_data tds[num_threads];
    
    for( t = 0; t < num_threads; t++ )
    {  
        memcpy( tds+t,  &td_template, sizeof(td_template) );
        tds[t].thread_id = t;
        
        pthread_attr_init(attrs + t);
        pthread_attr_setdetachstate(attrs + t, PTHREAD_CREATE_JOINABLE);
        
        rc = pthread_create( thread + t, 
                             attrs + t, 
                             find_candidate_mappings, 
                             (void *)(tds + t) 
            ); 
        
        if (rc) {
            fprintf(stderr, 
                    "ERROR; return code from pthread_create() is %d\n", 
                    rc
            );
            exit(-1);
        }
    }
    
    /* Free attribute and wait for the other threads */    
    size_t num_reads2 = 0;
    for(t=0; t < num_threads; t++) {
        rc = pthread_join(thread[t], &status);
        pthread_attr_destroy(attrs+t);
        if (rc) {
            fprintf( stderr, 
                     "ERROR; return code from pthread_join() is %d\n", 
                     rc
            );
            exit(-1);
        }
        num_reads2 += (size_t) status;
    }
    
    /* Find all of the candidate mappings */    
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Mapped (%i/%li) Partial Reads in %.2lf seconds ( %e/thread-hour )\n",
            mapped_cnt, rdb->readkey, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)mapped_cnt)*CLOCKS_PER_SEC*3600)/(stop-start)
        );
}
