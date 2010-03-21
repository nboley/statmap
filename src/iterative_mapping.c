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

#include "statmap.h"
#include "iterative_mapping.h"
#include "mapped_location.h"
#include "snp.h"

/******************************************************************************
 * Chr Trace Code
 *
 * Chr traces are just arrays that store expectations of cnts at each basepair. 
 * Below is code for manipulating them.
 *
 *****************************************************************************/

/* build an mmapped array to store the density */
TRACE_TYPE*
init_trace( size_t size  )
{
    void* chr_trace
        = mmap( NULL, sizeof(TRACE_TYPE)*size, 
                PROT_READ|PROT_WRITE, MAP_ANON|MAP_POPULATE|MAP_PRIVATE, 0, 0 );
    
    if( chr_trace == (void*) -1 )
    {
        fprintf(stderr, "Can not anonymous mmap of size '%u'\n", size );
        exit( -1 );
    }
    
    return (TRACE_TYPE*) chr_trace;
}

/* Build mmapped arrays for all of the chrs ***/
void
init_stranded_traces( const genome_data* const genome,
                      stranded_traces_t** traces )
{
    /* Allocate space for the struct */
    *traces = malloc( sizeof( stranded_traces_t ) );
    
    /* set the number of traces */
    (*traces)->num_traces = genome->num_chrs;

    /* Allocate space for the pointers to the chr's individual traces */
    (*traces)->fwd_traces = malloc( (*traces)->num_traces*sizeof(TRACE_TYPE*) );
    (*traces)->bkwd_traces = malloc( (*traces)->num_traces*sizeof(TRACE_TYPE*) );
    
    int i;
    for( i = 0; i < (*traces)->num_traces; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->trace_lengths[i] = genome->chr_lens[i];

        /* initialize the forward trace */
        (*traces)->fwd_traces[i] = init_trace( (*traces)->trace_lengths[i] );
        memset( (*traces)->fwd_traces + i, 0, 
                sizeof(TRACE_TYPE)*(*traces)->trace_lengths[i] );

        /* initialize the backwards trace */
        (*traces)->bkwd_traces[i] = init_trace( (*traces)->trace_lengths[i] );
        memset( (*traces)->bkwd_traces + i, 0, 
                sizeof(TRACE_TYPE)*(*traces)->trace_lengths[i] );
    }
    
    return;
}

void
close_stranded_traces( stranded_traces_t* traces )
{
    int i, error; 
    for( i = 0; i < traces->num_traces; i++ )
    {
        error = munmap( traces->fwd_traces[i], 
                        sizeof(TRACE_TYPE)*traces->trace_lengths[i] );
        if( error == -1 )
        {
            perror( "Problem closing mmaped trace" );
            assert( false );
            exit( -1 );
        }

        error = munmap( traces->bkwd_traces[i], 
                        sizeof(TRACE_TYPE)*traces->trace_lengths[i] );
        if( error == -1 )
        {
            perror( "Problem closing mmaped trace" );
            assert( false );
            exit( -1 );
        }
    }

    free( traces->trace_lengths );
    free( traces );

    return;
}

void
init_traces( const genome_data* const genome,
             traces_t** traces )
{
    /* Allocate space for the struct */
    *traces = malloc( sizeof( stranded_traces_t ) );
    
    /* set the number of traces */
    (*traces)->num_traces = genome->num_chrs;

    /* Allocate space for the pointers to the chr's individual traces */
    (*traces)->trace_lengths = malloc( 
        (*traces)->num_traces*sizeof(unsigned int) );
    (*traces)->traces = malloc( (*traces)->num_traces*sizeof(TRACE_TYPE*) );
    
    int i;
    for( i = 0; i < genome->num_chrs; i++ )
    {
        /* initialize the trace length, in bps */
        (*traces)->trace_lengths[i] = genome->chr_lens[i];

        /* initialize the forward trace */
        (*traces)->traces[i] = init_trace( (*traces)->trace_lengths[i] );
        memset( (*traces)->traces[i], 0, 
                sizeof(TRACE_TYPE)*(*traces)->trace_lengths[i] );
    }
    
    return;
}

void
close_traces( traces_t* traces )
{
    int i, error; 
    for( i = 0; i < traces->num_traces; i++ )
    {
        error = munmap( traces->traces[i], 
                        sizeof(TRACE_TYPE)*traces->trace_lengths[i] );
        if( error == -1 )
        {
            perror( "Problem closing mmaped trace" );
            assert( false );
            exit( -1 );
        }
    }

    free( traces->trace_lengths );
    free( traces->traces );
    free( traces );

    return;
}

void
zero_traces( traces_t* traces )
{
    /* Zero out the trace for the update */
    /* BUG FIXME TODO use memset! WTF? */
    long int i;
    #pragma omp parallel for num_threads( num_threads )
    for( i = 0; i < traces->num_traces; i++ )
    {
        unsigned int j;
        for( j = 0; j < traces->trace_lengths[i]; j++ )
        {
            traces->traces[i][j] = 0;
        }
    }
}


double
sum_traces( traces_t* traces )
{
    double sum = 0;

    int i;
    unsigned int j;
    for( i = 0; i < traces->num_traces; i++ )
        for( j = 0; j < traces->trace_lengths[i]; j++ )
            sum += traces->traces[i][j];
    
    return sum;
}


 /*
void
renormalize_traces( TRACE_TYPE** chr_traces, 
                   int num_chrs, 
                   unsigned int* chr_lens,
                   unsigned long num_reads )
{
    long int i;
    #pragma omp parallel for num_threads( num_threads )
    for( i = 0; i < num_chrs; i++ )
    {
        unsigned int j;
        for( j = 0; j < chr_lens[i]; j++ )
        {
            chr_traces[i][j] /= num_reads;
        }
    }
    
}
 */

void
write_wiggle_from_traces( traces_t* traces,
                          char** trace_names,
                          
                          const char* output_fname, 
                          const char* track_name,
                          
                          const double filter_threshold )
{    
    FILE* wfp = fopen( output_fname, "w" );
    /* Print out the header */
    fprintf( wfp, "track type=wiggle_0 name=%s\n", track_name );

    int i;
    unsigned int j;
    for( i = 0; i < traces->num_traces; i++ )
    {
        /* Print out the new chr start line */
        if( trace_names == NULL ) {
            fprintf( wfp, "variableStep chrom=%i\n", i );
        } else {
            fprintf( wfp, "variableStep chrom=%s\n", trace_names[i] );
        }
        
        
        for( j = 0; j < traces->trace_lengths[i]; j++ )
        {
            if( traces->traces[i][j] > filter_threshold )
                fprintf( wfp, "%i\t%e\n", j+1, traces->traces[i][j] );
        }
    }
    
    fclose( wfp );
}


/*
 * END Chr Trace Code
 *
 *****************************************************************************/

/*****************************************************************************
 * 
 * Mapped Reads Fns
 *
 * Methods for dealing with mapped short reads in the context of iterative 
 * updates.
 *
 *****************************************************************************/

/* use this for wiggles */
void
update_traces_from_read_densities( 
    mapped_reads_db* reads_db,
    traces_t* traces
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
        double cond_prob_sum = 0;
        for( j = 0; j < r.num_mappings; j++ )
        {
            int chr_index = r.locations[j].chr;
            unsigned int start = r.locations[j].start_pos;
            unsigned int stop = r.locations[j].stop_pos;
            ML_PRB_TYPE cond_prob = r.locations[j].cond_prob;
            cond_prob_sum += cond_prob;
            
            assert( cond_prob >= -0.0001 );
            assert( stop >= start );            
            assert( chr_index < traces->num_traces );
            assert( traces->trace_lengths[chr_index] >= stop );

            unsigned int k = 0;
            for( k = start; k < stop; k++ )
            {
                /* update the trace */
                #pragma omp atomic
                traces->traces[chr_index][k] 
                    += (1.0/(stop-start))*cond_prob;
            }
        }

        /* Make sure that the conditional probabilities sum to 1 */
        // fprintf( stderr, "%e\n", cond_prob_sum );
        // assert( cond_prob_sum > 0.95 && cond_prob_sum < 1.05 );
    }
    
    return;
}


void
update_traces_from_mapped_chipseq_reads( 
    mapped_reads_db* reads_db,
    traces_t* traces
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
        double cond_prob_sum = 0;
        for( j = 0; j < r.num_mappings; j++ )
        {
            int chr_index = r.locations[j].chr;
            unsigned char flag = r.locations[j].flag;
            unsigned int start = r.locations[j].start_pos;
            unsigned int stop = r.locations[j].stop_pos;
            ML_PRB_TYPE cond_prob = r.locations[j].cond_prob;
            cond_prob_sum += cond_prob;
            
            assert( cond_prob >= -0.0001 );
            assert( stop >= start );

            /* Make sure the reference genome is correct */
            //  FIXME - debugging code 
            /*
            printf( "%i:%u/%i - Chr Index: %i Loc: %u\n", 
                    (int)i, j, (int)reads_db->num_mmapped_reads, chr_index, stop  );
            */
            
            assert( chr_index < traces->num_traces );
            assert( traces->trace_lengths[chr_index] >= stop );

            /* update the trace */
            /* If the reads are paired */
            unsigned int k = 0;
            if( (flag&IS_PAIRED) != 0 )
            {
                for( k = start; k < stop; k++ )
                {
                    assert( chr_index < traces->num_traces );
                    assert( k < traces->trace_lengths[chr_index] );
                    
                    #pragma omp atomic
                    traces->traces[chr_index][k] 
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
                        traces->traces[chr_index][start + k] 
                            += k*cond_prob/( frag_len*frag_len  );
                    }
                    
                    if( stop >= k )
                    {
                        #pragma omp atomic
                        traces->traces[chr_index][stop - k] 
                            += k*cond_prob/( frag_len*frag_len  );
                    }
                }
            }
        }

        /* Make sure that the conditional probabilities sum to 1 */
        // fprintf( stderr, "%e\n", cond_prob_sum );
        // assert( cond_prob_sum > 0.95 && cond_prob_sum < 1.05 );
    }
    
    return;
}

double
update_mapped_chipseq_reads_from_trace( 
    mapped_reads_db* reads_db,
    traces_t* traces
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
        
        /* allocate space to store the temporary values */
        ML_PRB_TYPE* new_prbs = malloc( sizeof(double)*(r.num_mappings) );

        /* Update the reads from the trace */
        double density_sum = 0;
        unsigned int i;
        for( i = 0; i < r.num_mappings; i++ )
        {
            /* calculate the mean density */
            double window_density = 0;

            int chr_index = r.locations[i].chr;
            unsigned char flag = r.locations[i].flag;
            unsigned int start = r.locations[i].start_pos;
            unsigned int stop = r.locations[i].stop_pos;
            
            /* If the reads are paired */
            unsigned int j = 0;
            if( (flag&IS_PAIRED) > 1 )
            {
                for( j = start; j <= stop; j++ )
                {
                    assert( j > 0 && j < traces->trace_lengths[chr_index] );

                    #pragma omp atomic
                    window_density += traces->traces[chr_index][j];
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
                    window_density += traces->traces[chr_index][j];
                }
            }
            
            new_prbs[i] = r.locations[i].seq_error*window_density;
            density_sum += new_prbs[i];
        }
        
        /* renormalize the read probabilities */
        /* 
         * if we skip the density at the mapped bp, it's possible for the sum to be zero.
         * If this is the case, ignore the read ( and report at the end )
         * 
         *   BUG!!!! FIXME The 'report at the end' ( from above ) is not happening
         */
        if( density_sum > 0 )
        {
            for( i = 0; i < r.num_mappings; i++ )
            {
                /** Calculate the error, to check for convergence */
                /* quadratic error */
                // abs_error += pow(new_prbs[i]/density_sum 
                //                 - r.locations[i].cond_prob, 2 ) ;
                /* absolute value error */
                double normalized_density = new_prbs[i]/density_sum; 

                abs_error += MAX( normalized_density - r.locations[i].cond_prob,
                                 -normalized_density  + r.locations[i].cond_prob);
                
                r.locations[i].cond_prob = normalized_density;
            }
        }      

        /* Free the tmp array */
        free( new_prbs );
    }
    

    return abs_error;
}


int
update_chipseq_mapping( mapped_reads_db* rdb, 
                        genome_data* genome,
                        int max_num_iterations )
{
    clock_t start, stop;

    /*** Update the trace from the reads ***/
    
    // build the chr traces
    traces_t* traces;
    init_traces( genome,
                 &traces );
    
    double abs_error = 0;
    int num_iterations = 0;
    for( num_iterations = 0; 
         num_iterations < max_num_iterations; 
         num_iterations++ )
    {
        start = clock();
        
        /* Update the trace from the read probabilities  */
        update_traces_from_mapped_chipseq_reads( rdb, traces );
        
        /* Update the read probabilities from the trace */
        abs_error = update_mapped_chipseq_reads_from_trace( rdb, traces );
        
        stop = clock( );
        fprintf( stderr, "Iter %i: Error: %e \tUpdated trace in %.2f sec\tTrace Sum: %e\n", 
                 num_iterations, abs_error, 
                 ((double)stop-(double)start)/CLOCKS_PER_SEC,
                 sum_traces( traces )/rdb->num_mmapped_reads
            );
        
        if( abs_error < 1e-4 )
            break;
    }

    close_traces( traces );
    
    return 0;
}

int
update_mapping( mapped_reads_db* rdb, 
                genome_data* genome,
                int max_num_iterations,
                enum assay_type_t assay_type      )
{

    int error = 0;

    /* 
     * dispatch the correct update function. If no 
     * update code exists, then print a warning
     * and return without doing anything. 
     */

    switch( assay_type )
    {
    case CAGE:
        fprintf( stderr, "Cant currently iteratively update CAGE experiments." );
        break;
    case CHIP_SEQ:
        error = update_chipseq_mapping (
            rdb, genome, max_num_iterations 
        );
        break;
    default:
        fprintf( stderr, "WARNING     :  Can not iteratively map for assay type '%u'.\n", assay_type);
        fprintf( stderr, "WARNING     :  Returning marginal mappings.\n" );
        break;
    }

    return 0;
}

void
write_mapped_reads_to_wiggle( mapped_reads_db* rdb, 
                              genome_data* genome,
                              FILE* wfp )
{
    const double filter_threshold = 1e-6;

    /* Print out the header */
    fprintf( wfp, "track type=wiggle_0 name=%s\n", "cage_wig_track" );
    
    // build and update the chr traces
    traces_t* traces;
    init_traces( genome, &traces );


    update_traces_from_read_densities( rdb, traces );
    
    int i;
    unsigned int j;
    for( i = 0; i < traces->num_traces; i++ )
    {
        fprintf( wfp, "variableStep chrom=%s\n", genome->chr_names[i] );
        
        for( j = 0; j < traces->trace_lengths[i]; j++ )
        {
            if( traces->traces[i][j] >= filter_threshold )
                fprintf( wfp, "%i\t%e\n", j+1, traces->traces[i][j] );
        }
    }

    close_traces( traces );
}



/*
 * 
 * END Mapped Short Reads Fns
 *
 *****************************************************************************/
