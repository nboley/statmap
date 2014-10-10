#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
/* pow */
#include <math.h>

/* reading directory entries */
#include <dirent.h>

#include "wiggle.h"
#include "trace.h"
#include "genome.h"
#include "util.h"

#include "log.h"

static inline TRACE_TYPE 
sum( const TRACE_TYPE a, const TRACE_TYPE b )
{    return a + b; }

void
parse_meta_info( FILE* mi_fp, double** lhds )
{
    int rv;
    
    /* initialize the lhds array */
    *lhds = NULL;
    int num_lhds = 0;
    
    /* skip the header */
    rv = fscanf( mi_fp, "%*s" );
    assert( rv == 0 );
    
    /* store the maximum lhd for normalization */
    double max_log_lhd = -1e99;

    int i;
    for( i = 0; !feof( mi_fp ); i++ ) 
    {
        *lhds = realloc( *lhds, (i+1)*sizeof(double) );
        rv = fscanf( mi_fp, "%*i,%lf\n", (*lhds) + i );
        assert( rv == 1 );
        max_log_lhd = MAX( max_log_lhd, (*lhds)[i] );
    }
    num_lhds = i;

    /* normalize the lhds */
    double lhd_sum = 0;
    for( i = 0; i < num_lhds; i++ )
    {
        (*lhds)[i] = pow( 10, ((*lhds)[i] - max_log_lhd) );
        lhd_sum += (*lhds)[i];
    }

    /* normalize the sums to add to 1 */
    for( i = 0; i < num_lhds; i++ )
    {
        (*lhds)[i] /= lhd_sum;
    }
    
    return;
}

/* 
   Calculate the probability of observing >= cnt counts given that
   they are drawn from a 0.5 binomial 
*/
double calc_pvalue( double cnt )
{
    assert( 25 == NUM_BOOTSTRAP_SAMPLES );

    switch ( (int) cnt  )
    {
        /* generated from ../utilities/generate_binom_p_values.py */
    case 25:
        return 2.980232e-08;
    case 24:
        return 7.748604e-07;
    case 23:
        return 9.715557e-06;
    case 22:
        return 7.826090e-05;
    case 21:
        return 4.552603e-04;
    case 20:
        return 2.038658e-03;
    case 19:
        return 7.316649e-03;
    case 18:
        return 2.164263e-02;
    case 17:
        return 5.387607e-02;
    case 16:
        return 1.147615e-01;
    case 15:
        return 2.121781e-01;
    case 14:
        return 3.450190e-01;
    case 13:
        return 5.000000e-01;
    case 12:
        return 6.549810e-01;
    case 11:
        return 7.878219e-01;
    case 10:
        return 8.852385e-01;
    case 9:
        return 9.461239e-01;
    case 8:
        return 9.783574e-01;
    case 7:
        return 9.926834e-01;
    case 6:
        return 9.979613e-01;
    case 5:
        return 9.995447e-01;
    case 4:
        return 9.999217e-01;
    case 3:
        return 9.999903e-01;
    case 2:
        return 9.999992e-01;
    case 1:
        return 1.000000e+00;
    case 0:
        return 1.000000e+00;
    default:
        statmap_log( LOG_ERROR, "Illegal value '%i' encountered", (int) cnt );
        assert( false );
    }
    // silence the compiler
    return(-1);
}

double subtract_from_1( double cnt )
{
    return 1 - cnt;
}

void
update_ip_is_gt_nc_counts_from_traces( 
    struct trace_t* ip_trace, 
    struct trace_t* nc_trace, 
    struct trace_t* peaks_trace )
{
    int track, chr, segment, bp;
    for( track = 0; track < peaks_trace->num_tracks; track++ )
    {
        for( chr = 0; chr < peaks_trace->num_chrs; chr++ )
        {
            /* TODO can we assume the different traces will have the same layout
            of segments? For now, we do. */

            struct trace_segments_t* ip_segs = &(ip_trace->segments[track][chr]);
            struct trace_segments_t* nc_segs = &(nc_trace->segments[track][chr]);
            struct trace_segments_t* peak_segs
                = &(peaks_trace->segments[track][chr]);

            /* make sure all traces have the same set of segments */
            assert( ip_segs->num_segments == nc_segs->num_segments &&
                    nc_segs->num_segments == peak_segs->num_segments );

            for( segment = 0; segment < peak_segs->num_segments; segment++ )
            {
                struct trace_segment_t* ip_seg = ip_segs->segments + segment;
                struct trace_segment_t* nc_seg = nc_segs->segments + segment;
                struct trace_segment_t* peak_seg = peak_segs->segments + segment;

                /* make sure the segments match - simplifies things for now */
                assert( ip_seg->length == nc_seg->length &&
                        nc_seg->length == peak_seg->length );

                for( bp = 0; bp < peak_seg->length; bp++ )
                {
                    if( ip_seg->data[bp] > nc_seg->data[bp] ) {
                        peak_seg->data[bp] += 1;
                    } else if( ip_seg->data[bp] == nc_seg->data[bp] ) {
                        peak_seg->data[bp] += 0.5;
                    }
                }
            }
        }
    }
    
    return;
}

void
call_peaks_at_local_maxima( struct genome_data* genome, char* samples_dname ) 
{
    /* buffer for output/input filenames */
    char fname[500];
    
    /** find the likelihood **/
    /* This involves parsing hte directory name */
    /* we just move the pointer past 'sample', and then use atoi */
    int sample_id = atoi( samples_dname + 6 );
    if( 0 == sample_id ) {
        statmap_log( LOG_FATAL, "Could not convert directory name id" );
    }   

    char* track_names[2] = {"fwd_strand_peaks", "bkwd_strand_peaks"};     
    
    /* init the trace */
    struct trace_t* peaks_trace;
    /* two tracks correspond to pos and neg strands */
    init_full_trace( genome, &peaks_trace, 2, track_names, 1 );
    zero_traces( peaks_trace );
    
    int bi; /* bootstrap index */
    for( bi = 0; bi < NUM_BOOTSTRAP_SAMPLES; bi++ )
    {
        /* open the IP trace file */
        sprintf( fname, "%s/%s/bssample%i.ip.bin.trace", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, samples_dname, bi+1 );
        struct trace_t* ip_trace;
        load_trace_from_file( &ip_trace, fname );
        
        /* open the NC trace file */
        sprintf( fname, "%s/%s/bssample%i.nc.bin.trace", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, samples_dname, bi+1 );
        struct trace_t* nc_trace;
        load_trace_from_file( &nc_trace, fname );
        
        update_ip_is_gt_nc_counts_from_traces( ip_trace, nc_trace, peaks_trace );        
        
        close_traces( ip_trace );
        close_traces( nc_trace );
    }
    
    apply_to_trace( peaks_trace, calc_pvalue );
    apply_to_trace( peaks_trace, subtract_from_1 );
    
    sprintf( fname, "%s/sample%i.bin.trace", CALLED_PEAKS_OUTPUT_DIRECTORY, sample_id );
    write_trace_to_file( peaks_trace, fname );
    
    close_traces( peaks_trace );
    
    return;
}

void
call_peaks( struct genome_data* genome )
{
    /* open and parse the meta info ( the likelihoods ) */
    FILE* meta_info_fp = open_check_error( RELAXED_SAMPLES_META_INFO_FNAME, "r" );
    double* lhd_weights = NULL;
    parse_meta_info( meta_info_fp, &lhd_weights );
    fclose( meta_info_fp );

    DIR *dp;
    struct dirent *ep;

    /** Call peaks over bootstrap samples **/
    /* Each of these is *conditional* on the marginal distribution - next
       we will go back and combine them for a global peak list */
    dp = opendir ( BOOTSTRAP_SAMPLES_ALL_PATH );
    if( NULL == dp )
    {
        statmap_log( LOG_FATAL, "Could not open bootstrap samples directory." );
        exit( 1 );
    }

    int num_samples = 0;
    while( 0 != ( ep = readdir (dp) ) )
    {
        /* skip the references */
        if( ep->d_name[0] == '.' )
            continue;     
        
        num_samples += 1;
        
        call_peaks_at_local_maxima( genome, ep->d_name );
    }
    
    (void) closedir (dp);

    /** Combine all of the peaks */
    /* open a trace to aggregate over */
    char* track_names[2] = {"fwd_strand_peaks", "bkwd_strand_peaks"};
    struct trace_t* peaks_trace;
    init_full_trace( genome, &peaks_trace, 2, track_names, 1 );
    zero_traces( peaks_trace );    
    
    float lhds_sum = 0;
    
    int i;
    /* aggregate them into a 'final' wiggle */
    for( i = 0; i < num_samples; i++ )
    {
        char fname[500];
        sprintf( fname, "%s/sample%i.bin.trace", CALLED_PEAKS_OUTPUT_DIRECTORY, i+1 );

        struct trace_t* local_peaks;
        load_trace_from_file( &local_peaks, fname );
        multiply_trace_by_scalar( local_peaks, lhd_weights[i] );
        aggregate_over_trace_pairs( peaks_trace, local_peaks, sum );
        lhds_sum += lhd_weights[i];
        close_traces( local_peaks );
    }
    /* renormalize by the lhds */
    divide_trace_by_sum( peaks_trace, lhds_sum );
    
    
}
