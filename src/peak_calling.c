#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
/* pow */
#include <math.h>

/* reading directory entries */
#include <dirent.h>

#include "../src/wiggle.h"
#include "../src/trace.h"
#include "../src/genome.h"

/* 
   A bit of a hack - store the lhds in a gloally accessible variable so that 
   I can use the wiggle aggregation framework 
*/
double* global_lhd_weights;

void
parse_meta_info( FILE* mi_fp, double** lhds )
{
    /* initialize the lhds array */
    *lhds = NULL;
    int num_lhds = 0;
    
    /* skip the header */
    fscanf( mi_fp, "%*s" );
    
    /* store the maximum lhd for normalization */
    double max_log_lhd = -1e99;

    int i;
    for( i = 0; !feof( mi_fp ); i++ ) 
    {
        *lhds = realloc( *lhds, (i+1)*sizeof(double) );
        fscanf( mi_fp, "%*i,%lf\n", (*lhds) + i );
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

static FILE* 
open_check_error( char* fname, char* file_mode )
{
    FILE* tmp;
    tmp = fopen( fname, file_mode );
    if( tmp == NULL )
    {
        fprintf( stderr, "Error opening '%s\n'", fname );
        assert( false );
        exit( -1 );
    }
    return tmp;
}

float 
wig_lines_weighted_sum( 
    const struct wig_line_info* lines, const int ub, const int num_wigs  )
{
    assert( num_wigs > 0 );
    
    float sum = 0;
    int i;
    for( i = 0; i <= ub; i++ )
        sum += ( ( global_lhd_weights[ lines[i].file_index ] )*lines[i].value );
    
    return sum;
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
        fprintf( stderr, "Illegal value '%i' encountered.", (int) cnt );
        assert( false );
    }
}

double subtract_from_1( double cnt )
{
    return 1 - cnt;
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
    if( 0 == sample_id )
    {
        perror( "Could not convert directory name id" );
        exit( 1 );
    }   
     
    /* init the trace */
    struct trace_t* trace;
    /* two tracks correspond to pos and neg strands */
    init_traces( genome, &trace, 2 );
    zero_traces( trace );
    
    int bi; /* bootstrap index */
    for( bi = 0; bi < NUM_BOOTSTRAP_SAMPLES; bi++ )
    {
        /* open the IP wiggle file */
        sprintf( fname, "%s/%s/bssample%i.ip.wig", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, samples_dname, bi+1 );
        FILE* ip_wig = open_check_error( fname, "r" ); 
        
        /* open the NC wiggle file */
        sprintf( fname, "%s/%s/bssample%i.nc.wig", 
                 BOOTSTRAP_SAMPLES_ALL_PATH, samples_dname, bi+1 );
        FILE* nc_wig = open_check_error( fname, "r" ); 
        
        call_peaks_from_wiggles( ip_wig, nc_wig, genome, trace );
        
        fclose( ip_wig );
        fclose( nc_wig );
    }
    
    apply_to_trace( trace, calc_pvalue );
    apply_to_trace( trace, subtract_from_1 );
    
    char* track_names[2] = {"fwd_strand_peaks", "bkwd_strand_peaks"};
    
    sprintf( fname, "%s/sample%i.wig", CALLED_PEAKS_OUTPUT_DIRECTORY, sample_id );
    FILE* ofp = open_check_error( fname, "w" );
    write_wiggle_from_trace_to_stream( 
        trace, genome->chr_names, track_names, ofp, 0.0 );
    fclose( ofp );
    
    close_traces( trace );
    
    return;
}

void
call_peaks( struct genome_data* genome )
{
    /* open and parse the meta info ( the likelihoods ) */
    FILE* meta_info_fp = open_check_error( RELAXED_SAMPLES_META_INFO_FNAME, "r" );
    global_lhd_weights = NULL;
    parse_meta_info( meta_info_fp, &global_lhd_weights );
    fclose( meta_info_fp );

    DIR *dp;
    struct dirent *ep;

    /** Call peaks over bootstrap samples **/
    /* Each of these is *conditional* on the marginal distribution - next
       we will go back and combine them for a global peak list */
    dp = opendir ( BOOTSTRAP_SAMPLES_ALL_PATH );
    if( NULL == dp )
    {
        perror( "FATAL       : Could not open bootstrap samples directory." );
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
    /* open all of the wiggles */
    FILE** fps = malloc( sizeof(FILE*)*num_samples );
    int i;
    for( i = 0; i < num_samples; i++ )
    {
        char fname[500];
        sprintf( fname, "%s/sample%i.wig", CALLED_PEAKS_OUTPUT_DIRECTORY, i+1 );
        fps[i] = open_check_error( fname, "r" );
    }
    
    /* aggregate them into a 'final' wiggle */
    FILE* peaks_fp = open_check_error( JOINED_CALLED_PEAKS_FNAME, "w" );
    aggregate_over_wiggles( fps, num_samples, peaks_fp, 0, wig_lines_weighted_sum );
    fclose( peaks_fp );
    
    /* close the sample wig fp's */
    for( i = 0; i < num_samples; i++ )
        fclose( fps[i] );
    free( fps );
}
