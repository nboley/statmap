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
    FILE* log_fp = td->log_fp;
    pthread_mutex_t* log_fp_mutex = td->log_fp_mutex;
    
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
    int max_subseq_len = td->max_subseq_len;

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

            /* If we are logging, print the read */
            if( log_fp != NULL ) {
                pthread_mutex_lock( log_fp_mutex );
                fprintf(log_fp, "%li:  %s\t%.*s\t%.*s", 
                        readkey, r->name, 
                        r->length, r->char_seq, 
                        r->length, r->error_str );
                fflush( log_fp );
                pthread_mutex_unlock( log_fp_mutex );
            }
                
            /**** go to the index for mapping locations */
            mapped_locations *results;
            init_mapped_locations( &results );

            search_index( genome->index, 

                          min_match_penalty,
                          max_penalty_spread,

                          results,

                          r,
                          bp_mut_rates
                );

            /* increment the number of reads that mapped */
            /* 
             *  FIXME - this isnt locked with the read key, 
             *  so the mapped cnt and cnt may be out of sync
             *  due to the threads. That is, when we print out
             *  status reports they are typically off by a read 
             *  or 2.
             *  
             *  FIXME - Also, what is the cost of these locks? I wonder if,
             *  since the counts are off anyways, it may be better 
             *  to let them accumulate and then update them out 
             *  every, for instance, 100 reads.
             *  
             */
            pthread_mutex_lock( mapped_cnt_mutex );
            *mapped_cnt += ( results->length > 0 ) ? true : false;
            pthread_mutex_unlock( mapped_cnt_mutex );

            if( log_fp != NULL ) {
                pthread_mutex_lock( log_fp_mutex );
                fprintf(log_fp, "\t%i\t%e\n", 
                        (int) results->length, 
                        readkey/(((double)(clock()-start))/CLOCKS_PER_SEC) ); 

                fflush( log_fp );
                pthread_mutex_unlock( log_fp_mutex );
            }

            /****** Prepare the template candidate_mapping objects ***********/
            candidate_mapping template_candidate_mapping 
                = init_candidate_mapping_from_template( 
                    r, max_subseq_len, max_penalty_spread 
            );
            assert( template_candidate_mapping.rd_type != 0 );

            /**** TODO - get rid of this requirement */
            /* Assert that subseqs are false */
            assert( template_candidate_mapping.subseq_len == r->length );

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
                template_candidate_mapping.start_bp = (result->location).loc;
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
            int j;
            for( j = 0; j < mappings->length; j++ )
            {
                /* recheck location */
                
                
                /* recheck penalty */
                /* 
                 * We always need to do this because of the way the search queue
                 * handles poor branches. If our brnach prediction fails we could
                 * add a low probability read, and then go back and find a very 
                 * good read which would invalidate the previous. Then, the read
                 * may not belong, but it isnt removed in the index searching method.
                 */
                /* I set the safe bit to 1e-6, which is correct for a float */
                if(  max_penalty_spread > -0.00001 
                     &&
                     ( (mappings->mappings + j)->penalty 
                       < ( max_penalty - max_penalty_spread ) ) )
                {
                    (mappings->mappings + j)->recheck = INVALID;
                } else {
                    (mappings->mappings + j)->recheck = VALID;
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
