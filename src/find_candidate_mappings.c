/* Copyright (c) 2009-2010, Nathan Boley */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>

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

    genome_data* genome = td->genome;
    FILE* log_fp = td->log_fp;
    pthread_mutex_t* log_fp_mutex = td->log_fp_mutex;
    
    unsigned int* mapped_cnt = td->mapped_cnt;
    pthread_mutex_t* mapped_cnt_mutex = td->mapped_cnt_mutex;

    rawread_db_t* rdb = td->rdb;
    unsigned int* read_key = td->read_key;
    pthread_mutex_t* ifp_mutex = td->ifp_mutex;
    
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

    int error;
    clock_t start;
    start = clock();
    
    /* build the basepair mutation rate lookup table */
    float* bp_mut_rates;
    determine_bp_mut_rates( &bp_mut_rates );
    
    /* to make this thread safe, we take the lock inside of 
       the while loop and break if the raw read db is empty  */
    while( 1  ) 
    {        
        /****** get the read that we want to map ******/
        rawread *r1, *r2;
        
        pthread_mutex_lock( ifp_mutex );
        (*read_key)++;
        unsigned int curr_read_key = *read_key;
        
        /* check to see if ifp is empty. if it is, break out of 
           the while loop */
        /* RDB */
        if( rawread_db_is_empty( rdb ) )
        {
            (*read_key)--;
            pthread_mutex_unlock( ifp_mutex );
            break;
        }

        /* get the next read in the file */
        error = 
            get_next_read_from_rawread_db( rdb, &r1, &r2 );
        
        /* release the lock */
        pthread_mutex_unlock( ifp_mutex );
        
        if( curr_read_key % 10000 == 0 && curr_read_key > 0 )
        {
            pthread_mutex_lock( mapped_cnt_mutex );
            fprintf(stderr, "DEBUG       :  Mapped %i reads, %i successfully\n", 
                    curr_read_key, *mapped_cnt);
            pthread_mutex_unlock( mapped_cnt_mutex );
        }

        /* Test to see if the read is mappable. If it is not, continue */
        /* 
         * Note that both pairs of a paired end read have to be mappable
         * for the read to be mappable.
         */
        /* If this is a paired end read */
        if( r2 != NULL )
        {
            if ( filter_rawread( r1 ) == true
                 || filter_rawread( r2 ) == true )
            {
                free_rawread( r1 );
                free_rawread( r2 );
                continue;
            }
        } 
        /* Otherwise, this is a single ended read */
        else {
            if( filter_rawread( r1 ) == true  )
            {   
                fprintf_rawread_to_fastq( 
                    rdb->unmappable_single_end_reads, r1 );
                free_rawread( r1 );
                continue;
            }
        }

        /* consider both read pairs */
        int j = 0;
        rawread* reads[2] = { r1, r2 };
        for( j = 0; j < 2 && reads[j] != NULL; j++ )
        {
            rawread* r = reads[j];

            /* If we are logging, print the read */
            if( log_fp != NULL ) {
                pthread_mutex_lock( log_fp_mutex );
                fprintf(log_fp, "%i:  %s\t%.*s\t%.*s", 
                        curr_read_key, r->name, 
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
             *  BUG - this isnt locked with the read key, 
             *  so the mapped cnt and cnt may be out of sync
             *  due to the threads. That is, when we print out
             *  status reports they are typically off by a read 
             *  or 2.
             */
            pthread_mutex_lock( mapped_cnt_mutex );
            *mapped_cnt += ( results->length > 0 ) ? true : false;
            pthread_mutex_unlock( mapped_cnt_mutex );

            if( log_fp != NULL ) {
                pthread_mutex_lock( log_fp_mutex );
                fprintf(log_fp, "\t%i\t%e\n", 
                        (int) results->length, 
                        curr_read_key/(((double)(clock()-start))/CLOCKS_PER_SEC) ); 

                fflush( log_fp );
                pthread_mutex_unlock( log_fp_mutex );
            }

            /****** Prepare the template candidate_mapping objects ***********************************/
            candidate_mapping template_candidate_mapping 
                = init_candidate_mapping_from_template( 
                    r, max_subseq_len, max_penalty_spread 
            );

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
             * out penalties that are too low. The index will never return results
             * that are strictly below the minimum penalty, but it may return 
             * results below the relative penalty. ( Read the indexing header for
             * details )
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

                if( result->location.snp_coverage != 0 )
                {
                    template_candidate_mapping.recheck = COVERS_SNP;
                }

                add_candidate_mapping( mappings, &template_candidate_mapping );
            }

            // print_mapped_locations( results );
            // print_candidate_mappings( mappings );
            // printf( "==========================================\n\n" );

            free_mapped_locations( results );

            /****** Do the recheck ******/
            /* 
             * Currently, everything should be set *except* the gene strand. This
             * is because we dont know what it is. Therefore, we will add these to the db
             * with that bit unset, and then during the merging stage add to the penalty
             * and say that it was equally likely to have come from either gene strand.
             * This corresponds with us being equally certain that the read is 
             * from either gene.
             *
             * Also, we need to check that we dont have any low quality reads. ( This is 
             *   possible if the path search went awry, and we found a low quality read 
             *   before a high quality read )
             */
            for( i = 0; i < mappings->length; i++ )
            {
                /* recheck location */
                /* FIXME - TODO */
                if( (mappings->mappings + i)->recheck == COVERS_SNP )
                {
                    // find all snps that cover this read
                    // read start = (mappings->mappings + i)->start_bp
                    // read_length = (mappings->mappings + i)->rd_len
                    int num_snps;
                    int* snps;

                    find_snps_in_snp_db( 
                        genome->snp_db, 
                        (mappings->mappings + i)->chr,
                        (mappings->mappings + i)->start_bp,
                        (mappings->mappings + i)->start_bp
                        + (mappings->mappings + i)->rd_len,
                        &num_snps,
                        &snps
                    );

                    printf("FOUND %i SNP READ(s)\n", num_snps);
                    int m = 0;
                    for( m = 0; m < num_snps; m++ )
                    {
                        printf( "Snp index: %i %i\n", snps[m], num_snps );
                    }

                    free( snps );
                }

                /* BUG - rounding error potentially */
                /* recheck penalty */
                if(  max_penalty_spread > -0.01 
                     &&
                     ( (mappings->mappings + i)->penalty 
                       < ( max_penalty - max_penalty_spread ) ) )
                {
                    (mappings->mappings + i)->recheck = INVALID;
                } else {
                    (mappings->mappings + i)->recheck = VALID;
                }
            }

            /****** add the results to the database ******/
            /* BUG - assert the readnames are identical */
            pthread_mutex_lock( mappings_db_mutex );
            assert( thread_id < NUM_THREADS );
            add_candidate_mappings_to_db( 
                mappings_db, mappings, thread_id, r->name, curr_read_key );
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
find_all_candidate_mappings( genome_data* genome,
                             FILE* log_fp,
                             rawread_db_t* rdb,

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
    /* 
     * This locks the input file(s). It also locks read_key - the global read cnt variable.
     * We need read_key in order to rejoin the reads after the candidate mappings are
     * found. This lets us maintain the read order despite the fact that 
     * the reads are mapped in separate threads.
     */
    pthread_mutex_t ifp_mutex = PTHREAD_MUTEX_INITIALIZER;
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

    unsigned int read_key = 0;
    td_template.read_key = &read_key;
    td_template.ifp_mutex = &ifp_mutex;
    
    td_template.mappings_db = mappings_db;
    td_template.mappings_db_mutex = &mappings_db_mutex;
    
    td_template.min_match_penalty = min_match_penalty;
    td_template.max_penalty_spread = max_penalty_spread;
    td_template.max_subseq_len = max_seq_length;

    /* initialize the threads */
    
    int num_threads = NUM_THREADS;
    long t;
    int rc;
    void* status;
    pthread_t thread[num_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    struct single_map_thread_data* tds;
    tds = malloc(num_threads*sizeof(td_template));
    for( t = 0; t < num_threads; t++ )
    {  
        memcpy( tds +t,  &td_template, sizeof(td_template) );
        tds[t].thread_id = t;
        rc = pthread_create( &(thread[t]), 
                             &attr, 
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
    pthread_attr_destroy(&attr);
    for(t=0; t<NUM_THREADS; t++) {
        rc = pthread_join(thread[t], &status);
        if (rc) {
            fprintf( stderr, 
                     "ERROR; return code from pthread_join() is %d\n", 
                     rc
            );
            exit(-1);
        }
        num_reads2 += (size_t) status;
    }
    
    /* Free the thread parameter data structures */
    free( tds );
    
    /* Find all of the candidate mappings */    
    clock_t stop = clock();
    fprintf(stderr, "PERFORMANCE :  Mapped (%i/%i) Partial Reads in %.2lf seconds ( %e/thread-hour )\n",
            mapped_cnt, read_key, 
            ((float)(stop-start))/CLOCKS_PER_SEC,
            (((float)mapped_cnt)*CLOCKS_PER_SEC*3600)/(stop-start)
        );

}
