#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <fcntl.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <errno.h>
#include <float.h>

#include "statmap.h"
#include "iterative_mapping.h"
#include "mapped_location.h"
#include "snp.h"

static TRACE_TYPE 
min( TRACE_TYPE a, TRACE_TYPE b )
{
    if( a > b )
        return b;
    return a;
}

static TRACE_TYPE 
max( TRACE_TYPE a, TRACE_TYPE b )
{
    if( a > b )
        return a;
    return b;
}

void
naive_update_trace_expectation_from_location( 
    const traces_t* const traces, 
    const mapped_read_location* const loc )
{
    int chr_index = loc->chr;
    unsigned int start = loc->start_pos;
    unsigned int stop = loc->stop_pos;
    ML_PRB_TYPE cond_prob = loc->cond_prob;

    unsigned int i;
    for( i = start; i <= stop; i++ )
        traces->traces[0][chr_index][i] += cond_prob;

    return;
}

void
update_traces_from_mapped_reads( 
    struct mapped_reads_db* reads_db,
    traces_t* traces,
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc)
)
{    
    /* zero traces */
    zero_traces( traces );
    
    /* Update the trace from the reads */
    /* 
     * FIXME - openmp ( for some reason ) throws a warning on an unsigned iteration 
     * variable. I dont understand why, but I am a bit nervous because the spec used to
     * be that iteration variable *had* to be unsigned so, I am wasting a ton of 
     * register space and upcasting all of the unsigned longs to long longs. When 
     * I move to open mp 3.0, I want to remove this 
     */

    long long i;
    const int chunk = 10000;
    #pragma omp parallel for schedule(dynamic, chunk) num_threads( num_threads )
    for( i = 0; i < (long long) reads_db->num_mmapped_reads; i++ )
    {
        char* read_start = reads_db->mmapped_reads_starts[i];

        /* read a mapping into the struct */
        mapped_read r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (mapped_read_location*) read_start;
        
        /* Update the trace from this mapping */
        unsigned int j;
        for( j = 0; j < r.num_mappings; j++ ) {
            update_trace_expectation_from_location( traces, r.locations + j );
        }
    }
    
    return;
}

double
update_mapped_reads_from_trace( 
    struct mapped_reads_db* reads_db,
    traces_t* traces,
    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
    )
{    
    /* store the total accumulated error */
    double abs_error = 0;

    /* 
     * FIXME - openmp ( for some reason ) throws a warning on an unsigned iteration 
     * variable. I dont understand why, but I am a bit nervous because the spec used to
     * be that iteration variable *had* to be unsigned so, I am wasting a ton of 
     * register space and upcasting all of the unsigned longs to long longs. When 
     * I move to open mp 3.0, I want to remove this 
     */
    long long k;
    const int chunk = 10000;
    #pragma omp parallel for schedule(dynamic, chunk) reduction(+:abs_error) num_threads( num_threads )
    for( k = 0; k < (long long) reads_db->num_mmapped_reads; k++ )
    {
        char* read_start = reads_db->mmapped_reads_starts[k];

        /* read a mapping into the struct */
        mapped_read r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (mapped_read_location*) read_start;

        abs_error += update_mapped_read_prbs( traces, &r );        
    }
    

    return abs_error;
}

int
update_mapping(
    struct mapped_reads_db* rdb,
    traces_t* starting_trace,
    int max_num_iterations,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),

    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
    )
{
    clock_t start, stop;
    
    double abs_error = 0;
    int num_iterations = 0;
    for( num_iterations = 0; 
         num_iterations < max_num_iterations; 
         num_iterations++ )
    {
        start = clock();

        update_traces_from_mapped_reads( 
            rdb, starting_trace, 
            update_trace_expectation_from_location
        );

        abs_error = update_mapped_reads_from_trace(
            rdb, starting_trace, 
            update_mapped_read_prbs
        );
        
        stop = clock( );
        fprintf( stderr, "Iter %i: Error: %e \tUpdated trace in %.2f sec\tTrace Sum: %e\n", 
                 num_iterations, abs_error, 
                 ((double)stop-(double)start)/CLOCKS_PER_SEC,
                 sum_traces( starting_trace )/rdb->num_mmapped_reads
            );
        
        if( abs_error < 1e-2 )
            break;
    }
    
    return 0;
}

void
build_random_starting_trace( 
    traces_t* traces, 
    struct mapped_reads_db* rdb,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),

    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
    )
{
    /* zero traces */
    zero_traces( traces );
    
    /* Update the trace from the reads */
    /* 
     * FIXME - openmp ( for some reason ) throws a warning on an unsigned iteration 
     * variable. I dont understand why, but I am a bit nervous because the spec used to
     * be that iteration variable *had* to be unsigned so, I am wasting a ton of 
     * register space and upcasting all of the unsigned longs to long longs. When 
     * I move to open mp 3.0, I want to remove this 
     */

    long long i;
    for( i = 0; i < (long long) rdb->num_mmapped_reads; i++ )
    {
        char* read_start = rdb->mmapped_reads_starts[i];
        
        /* read a mapping into the struct */
        mapped_read r;
        r.read_id = *((unsigned long*) read_start);
        read_start += sizeof(unsigned long)/sizeof(char);
        r.num_mappings = *((unsigned short*) read_start);
        read_start += sizeof(unsigned short)/sizeof(char);
        r.locations = (mapped_read_location*) read_start;

        /* reset the read cond prbs under a uniform prior */
        reset_read_cond_probs( &r  );
        
        /* update the read conditional probabilities from the trace */
        update_mapped_read_prbs( traces, &r );        

        /* Choose a random location, and use this */
        unsigned int j;
        float random_num = (float)rand()/(float)RAND_MAX;
        float cum_dist = 0;
        for( j = 0; j < r.num_mappings; j++ ) 
        {
            cum_dist += r.locations[j].cond_prob;
            if( random_num <= cum_dist )
            {
                float tmp_cond_prb = r.locations[j].cond_prob;
                r.locations[j].cond_prob = 1.0;
                update_trace_expectation_from_location( traces, r.locations + j );
                r.locations[j].cond_prob = tmp_cond_prb;
                break;
            }
        }
        
        /* deal with potential rounding errors */
        /* If this happens, no read would have been added */
        if( random_num > cum_dist )
        {
                float tmp_cond_prb = r.locations[r.num_mappings-1].cond_prob;
                r.locations[r.num_mappings-1].cond_prob = 1.0;
                update_trace_expectation_from_location( 
                    traces, r.locations + r.num_mappings - 1 );
                r.locations[r.num_mappings-1].cond_prob = tmp_cond_prb;
        }
    }
    
    return;
}


/*
 * END Chr Trace Code
 *
 *****************************************************************************/


/******************************************************************************
 *
 * High-Level Functions 
 *
 *****************************************************************************/

int
sample_random_traces( 
    struct mapped_reads_db* rdb, 
    struct genome_data* genome,
    int trace_dim,
    int num_samples,
    int max_num_iterations,
    
    void (* const update_trace_expectation_from_location)(
        const traces_t* const traces, 
        const mapped_read_location* const loc),
    
    double (* const update_mapped_read_prbs)( const traces_t* const traces, 
                                              const mapped_read* const r  )
                          
)
{
    /* Seed the random number generator */
    srand ( time(NULL) );

    /* initialize the max and min traces */
    traces_t *max_trace, *min_trace;
    init_traces( genome, &min_trace, trace_dim );
    init_traces( genome, &max_trace, trace_dim );
            
    /* build bootstrap samples. Then take the max and min. */
    int i;
    for( i = 0; i < num_samples; i++ )
    {
        traces_t* sample_trace;
        init_traces( genome, &sample_trace, trace_dim );

        build_random_starting_trace( 
            sample_trace, rdb, 
            update_trace_expectation_from_location,
            update_mapped_read_prbs
        );
        
        if( SAVE_STARTING_SAMPLES )
        {
            char buffer[100];
            sprintf( buffer, "%ssample%i.wig", STARTING_SAMPLES_PATH, i+1 );
            
            int j;
            for( j = 0; j < trace_dim; j++ )
            {
                char buffer2[100];
                sprintf( buffer2, "start_sample_%i_track_%i", i+1, j+1 );

                write_wiggle_from_traces( 
                    sample_trace, j, genome->chr_names, 
                    buffer2, "start_sample", 1e-2 );
            }
        }

        /* update the mapping */
        update_mapping( 
            rdb, sample_trace, max_num_iterations,
            update_trace_expectation_from_location,
            update_mapped_read_prbs
        );
        
        if( SAVE_SAMPLES )
        {
            char buffer[100];
            sprintf( buffer, "%ssample%i.wig", RELAXED_SAMPLES_PATH, i+1 );

            int j;
            for( j = 0; j < trace_dim; j++ )
            {
                char buffer2[100];
                sprintf( buffer2, "relaxed_sample_%i_track_%i", i+1, j+1 );

                write_wiggle_from_traces( 
                    sample_trace, j, genome->chr_names, 
                    buffer2, "relaxed_sample", 1e-2 );
            }
        }

        aggregate_over_traces( max_trace, sample_trace, max );

        aggregate_over_traces( min_trace, sample_trace, min );

        close_traces( sample_trace );
        
        printf( "Sample %i\n", i+1 );
    }

    int j;
    for( j = 0; j < trace_dim; j++ )
    {
        char buffer[100];
        sprintf( buffer, "max_trace_track_%i", j+1 );
        write_wiggle_from_traces( max_trace, 0, genome->chr_names, 
                                  "max_trace.wig", buffer, 1e-2 );

        sprintf( buffer, "min_trace_track_%i", j+1 );
        write_wiggle_from_traces( min_trace, 0, genome->chr_names, 
                                  "min_trace.wig", buffer, 1e-2 );
        
    }
    
    close_traces( max_trace );
    close_traces( min_trace );
    
    return 0;
}


/*
 *
 * END High-Level Functions 
 *
 ******************************************************************************/


/*****************************************************************************
 * 
 * ChIP Seq specific functions 
 *
 *****************************************************************************/

void 
update_chipseq_trace_expectation_from_location(
    const traces_t* const traces, 
    const mapped_read_location* const loc )
{
    int chr_index = loc->chr;
    unsigned char flag = loc->flag;
    unsigned int start = loc->start_pos;
    unsigned int stop = loc->stop_pos;
    ML_PRB_TYPE cond_prob = loc->cond_prob;
    
    assert( cond_prob >= -0.0001 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_traces );
    assert( traces->trace_lengths[chr_index] >= stop );
    
    /* update the trace */
    /* If the reads are paired */
    unsigned int k = 0;
    if( (flag&IS_PAIRED) != 0 )
    {
        for( k = start; k < stop; k++ )
        {
            assert( chr_index < traces->num_chrs );
            assert( k < traces->trace_lengths[chr_index] );
            
            #pragma omp atomic
            traces->traces[0][chr_index][k] 
                += (1.0/(stop-start))*cond_prob;
        }
    } 
    /* If the read is *not* paired */
    else {
        /* FIXME - hope that we actually have a fragment length */
        /* FIXME get the real fragment length */
        /* FIXME - cleanup the stop condition */
        unsigned int frag_len = 400;
        for( k = 0; k < frag_len; k++ )
        {
            if( start + k < traces->trace_lengths[chr_index] )
            {
                #pragma omp atomic
                traces->traces[0][chr_index][start + k] 
                    += k*cond_prob/( frag_len*frag_len  );
            }
            
            if( stop >= k )
            {
                #pragma omp atomic
                traces->traces[0][chr_index][stop - k] 
                    += k*cond_prob/( frag_len*frag_len  );
            }
        }
    }

    return;
}

double 
update_chipseq_mapped_read_prbs( const traces_t* const traces, 
                                 const mapped_read* const r  )
{
    double abs_error = 0;
    
    /* allocate space to store the temporary values */
    ML_PRB_TYPE* new_prbs = malloc( sizeof(double)*(r->num_mappings) );
    
    /* Update the reads from the trace */
    double density_sum = 0;
    unsigned int i;
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* calculate the mean density */
        double window_density = 0;
        
        int chr_index = r->locations[i].chr;
        unsigned char flag = r->locations[i].flag;
        unsigned int start = r->locations[i].start_pos;
        unsigned int stop = r->locations[i].stop_pos;
        
        /* If the reads are paired */
        unsigned int j = 0;
        if( (flag&IS_PAIRED) > 1 )
        {
            for( j = start; j <= stop; j++ )
            {
                assert( j > 0 && j < traces->trace_lengths[chr_index] );
                
                window_density += traces->traces[0][chr_index][j];
            }
        } 
        /* If the read is *not* paired */
        else {
            /* FIXME - hope that we actually have a fragment length */
            /* FIXME get the real fragment length */
            /* FIXME - cleanup the stop condition */
            unsigned int frag_len = 400;
            for( j = start; 
                 j < MIN(traces->trace_lengths[chr_index], start + frag_len); 
                 j++ )
            {
                window_density += traces->traces[0][chr_index][j];
            }
        }
        
        new_prbs[i] = r->locations[i].seq_error*window_density;
        density_sum += new_prbs[i];
    }
    
    /* renormalize the read probabilities */
    /* 
     * if we skip the density at the mapped bp, it's possible for the sum to 
     * be zero. If this is the case, ignore the read ( and report at the end )
     * 
     *   BUG!!!! FIXME The 'report at the end' ( from above ) is not happening
     */
    if( density_sum > 0 )
    {
        for( i = 0; i < r->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* absolute value error */
            double normalized_density = new_prbs[i]/density_sum; 
            
            abs_error += MAX( normalized_density - r->locations[i].cond_prob,
                              -normalized_density  + r->locations[i].cond_prob);
            
            r->locations[i].cond_prob = normalized_density;
        }
    }      
    
    /* Free the tmp array */
    free( new_prbs );

    return abs_error;
}

/*****************************************************************************
 * 
 * CAGE specific functions 
 *
 *****************************************************************************/

void update_CAGE_trace_expectation_from_location(
    const traces_t* const traces, 
    const mapped_read_location* const loc )
{
    int chr_index = loc->chr;
    unsigned char flag = loc->flag;
    unsigned int start = loc->start_pos;
    unsigned int stop = loc->stop_pos;
    ML_PRB_TYPE cond_prob = loc->cond_prob;
    
    assert( cond_prob >= -0.0001 );
    assert( stop >= start );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_traces );
    assert( traces->trace_lengths[chr_index] >= stop );

    assert( cond_prob >= -0.0001 );
    
    /* Make sure the reference genome is correct */            
    assert( chr_index < traces->num_traces );
    
    /* update the trace */
    /* If the reads are paired */
    if( (flag&IS_PAIRED) != 0 )
    {
        fprintf( stderr, "FATAL: paired cage reads are not supported" );
        exit( -1 );                
    } 
    /* If the read is *not* paired */
    else {
        /* store the trace that we care about */
        TRACE_TYPE** trace;
        
        /* If this is in the fwd strnanded transcriptome */
        if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
        {
            trace = traces->traces[0];
        } 
        /* We are in the 3' ( negative ) transcriptome */
        else {
            trace = traces->traces[1];
        }
        
        trace[ chr_index ][ start ] += cond_prob; 
    }
    
    return;
}

double update_CAGE_mapped_read_prbs( 
    const traces_t* const traces, 
    const mapped_read* const r  )
{
    double abs_error = 0;
    
    /* allocate space to store the temporary values */
    ML_PRB_TYPE* new_prbs = malloc( sizeof(double)*(r->num_mappings) );
    
    /* Update the reads from the trace */
    double density_sum = 0;
    unsigned int i;
    for( i = 0; i < r->num_mappings; i++ )
    {
        /* calculate the mean density */
        /* We set this to 2*DBL_EPSILON to prevent the division by 0 */
        double window_density = 2*DBL_EPSILON;
        
        int chr_index = r->locations[i].chr;
        unsigned char flag = r->locations[i].flag;
        unsigned int start = r->locations[i].start_pos;
        
        /* If the reads are paired */
        unsigned int j = 0;
        if( (flag&IS_PAIRED) > 0 )
        {
            fprintf( stderr, "FATAL: paired cage reads are not supported" );
            exit( -1 );
        } 
        /* If the read is *not* paired */
        else {
            /* store the trace that we care about */
            TRACE_TYPE** trace;
            
            /* If this is in the fwd strnanded transcriptome */
            if( flag&FIRST_READ_WAS_REV_COMPLEMENTED )
            {
                trace = traces->traces[0];
            } 
            /* We are in the 3' ( negative ) transcriptome */
            else {
                trace = traces->traces[1];
            }
            
            for( j = start;
                 j < MIN( traces->trace_lengths[chr_index], 
                          start + WINDOW_SIZE ); 
                 j++ )
            {
                window_density += trace[chr_index][j];
            }
        }
        
        new_prbs[i] = r->locations[i].seq_error*window_density;
        density_sum += new_prbs[i];
    }
    
    /* renormalize the read probabilities */
    if( density_sum > 0 )
    {
        for( i = 0; i < r->num_mappings; i++ )
        {
            /** Calculate the error, to check for convergence */
            /* quadratic error */
            // abs_error += pow(new_prbs[i]/density_sum 
            //                 - r.locations[i].cond_prob, 2 ) ;
            /* absolute value error */
            double normalized_density = new_prbs[i]/density_sum; 
            
            abs_error += MAX( normalized_density - r->locations[i].cond_prob,
                              -normalized_density  + r->locations[i].cond_prob);
            
            r->locations[i].cond_prob = normalized_density;
        }
    }      
    
    /* Free the tmp array */
    free( new_prbs );   

    return abs_error;
}


/*
 * END Cage Mapping Code
 *
 ******************************************************************************/

/******************************************************************************
 *
 * Entry point(s) into the iterative mapping code
 *
 */



int
generic_update_mapping( struct mapped_reads_db* rdb, 
                        struct genome_data* genome,
                        enum assay_type_t assay_type      )
{
    const int max_num_iterations = 500;
    const int num_samples = 100;

    int error = 0;
    
    /* 
     * dispatch the correct update function. If no 
     * update code exists, then print a warning
     * and return without doing anything. 
     */

    void (*update_expectation)(
        const traces_t* const traces, 
        const mapped_read_location* const loc) 
        = NULL;
    
    double (*update_reads)( const traces_t* const traces, 
                            const mapped_read* const r  )
        = NULL;

    
    int trace_size = -1;

    switch( assay_type )
    {
    case CAGE:
        update_expectation = update_CAGE_trace_expectation_from_location;
        update_reads = update_CAGE_mapped_read_prbs;
        trace_size = 2;
        break;
    
    case CHIP_SEQ:
        update_expectation = update_chipseq_trace_expectation_from_location;
        update_reads = update_chipseq_mapped_read_prbs;
        trace_size = 1;
        break;
    
    default:
        fprintf( stderr, "WARNING     :  Can not iteratively map for assay type '%u'.\n", assay_type);
        fprintf( stderr, "WARNING     :      Returning marginal mappings.\n" );
        return 1;
        break;
    }

    /* initialize the trace that we will store the expectation in */
    traces_t* uniform_trace;
    init_traces( genome, &uniform_trace, trace_size );

    error = update_mapping (
        rdb, 
        uniform_trace,
        max_num_iterations,
        update_expectation,
        update_reads
    );
    
    int j;
    for( j = 0; j < trace_size; j++ )
    {
        char buffer[100];
        sprintf( buffer, "expectation_track_%i", j+1 );
        
        write_wiggle_from_traces( 
            uniform_trace, j, genome->chr_names, 
            "relaxed_mapping.wig", buffer, 1e-2 );
    }
    
    close_traces( uniform_trace );
    
    error = sample_random_traces(
        rdb, genome, trace_size, num_samples, max_num_iterations, 
        update_expectation, update_reads
    );

    return 0;
    
}


/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/
