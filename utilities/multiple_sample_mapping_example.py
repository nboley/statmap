import psycopg2
import gzip
import subprocess
import os
import sys
from time import sleep

# assumption: this script's location is /path/to/statmap_dir/utilities/
STATMAP_PATH = os.path.abspath( os.path.join(os.path.dirname(__file__), "..") )

# Build paths to tools
STATMAP_BIN = os.path.join( STATMAP_PATH, "bin/statmap" )
BUILD_MARGINAL_BIN = os.path.join( STATMAP_PATH, "utilities/build_marginal_mappings_trace_from_mpd_read_db" )
BUILD_MIN_TRACE_BIN = os.path.join( STATMAP_PATH, "bin/build_min_trace" )
AGGREGATE_OVER_TRACES_BIN = os.path.join( STATMAP_PATH, "bin/aggregate_over_traces" )
CONVERT_TRACE_INTO_WIGGLE_CMD = os.path.join( STATMAP_PATH, "bin/convert_trace_into_wiggle" )

# Path to genome file
GENOME_PATH = "/media/scratch/genomes/drosophila/drosophila_all.genome"
#GENOME_PATH = "/media/scratch/genomes/D.pseudoobscura/genome.20.D.pseudoobscura"

# Performance configuration
NUM_THREAD_PER_PROC = 4
MAX_NUM_PROC = 7

# value of Statmap -n parameter to use for iterative mapping
NUM_STARTING_SAMPLES = 2

def read_samples_file( fname ):
    '''
    Return a list of tuples of (sn, lab_id, fns) from tab-delimited
    samples file
    '''
    samples = []
    with open( fname ) as fp:
        for line in fp:
            # split sample line into fields
            sn, lab_id, fns = line.strip().split()
            sn = sn.replace( ";", "_" )
            lab_id = lab_id.replace( ";", "_" )
            
            # list of comma-sep filenames
            fns = fns.split(',')
            samples.append( (sn, lab_id, fns) )

    return samples

def strip_and_cat_files( fnames, ofp, chars_to_strip ):
    for fname in fnames:
        fp = gzip.open( fname )
        for line_num, line in enumerate(fp):
            if line_num%2 == 1: 
                ofp.write( line[chars_to_strip:] )
            else: 
                ofp.write( line )
            
            # if line_num > 1002: break
        
        fp.close()
    
    return

def fastq_fname_from_sample_name( sample_name, bs_id ):
    return sample_name + "_" + bs_id + ".trimmed.fastq"

def output_dir_name_from_sample_name( sample_name, bs_id ):
    return sample_name + "_" + bs_id + "_mapped"

def log_fname_from_sample_name( sample_name, bs_id ):
    return sample_name + "_" + bs_id + ".run.log"

def run_statmap( sample_name, bs_id ):

    statmap_cmd = " ".join(
        (
            STATMAP_BIN,
            "-g %s" % GENOME_PATH,
            "-r %s" % fastq_fname_from_sample_name( sample_name, bs_id ),
            "-t %i" % NUM_THREAD_PER_PROC,
            "-o %s" % output_dir_name_from_sample_name( sample_name, bs_id ),
            "-n %i" % NUM_STARTING_SAMPLES,
            "-a a"
        )
    )
    
    # build wiggle from marginals
    marginal_wig_fname = os.path.join(
        output_dir_name_from_sample_name( sample_name, bs_id ),
        "marginal.wig" )

    build_marginal_cmd = "%s %s 1 > %s" % ( 
        BUILD_MARGINAL_BIN,
        output_dir_name_from_sample_name( sample_name, bs_id ),
        marginal_wig_fname
    )

    """
    # build the minimum trace over each sample using build_min_trace
    bootstrap_samples_cmd = ";\n".join([
        "%s min %s %s %i" % (
            BUILD_MIN_TRACE_BIN,
            "min_trace%i.bin.trace" % sample_number, # NOTE: build_min_trace cd's into smo_dir
            output_dir_name_from_sample_name( sample_name, bs_id ),
            sample_number
        )
        for sample_number in range(1, NUM_STARTING_SAMPLES+1)
    ])

    # aggregate min bootstrap traces from each sample into single min trace
    aggregated_min_trace_fname = os.path.join(
        output_dir_name_from_sample_name( sample_name, bs_id ),
        "aggregated_min.trace" )
    aggregate_bootstraps_cmd = "%s min %s %s" % (
        AGGREGATE_OVER_TRACES_BIN,
        aggregated_min_trace_fname,
        " ".join([
            os.path.join(
                output_dir_name_from_sample_name( sample_name, bs_id ),
                "min_trace%i.bin.trace" % sample_number
            )
            for sample_number in range(1, NUM_STARTING_SAMPLES+1)
        ])
    )
    """
    
    # covert single aggregated min trace to wiggle
    aggregated_min_wig_fname = os.path.join(
        output_dir_name_from_sample_name( sample_name, bs_id ),
        "aggregated_min.trace.wig" )
    convert_aggregated_min_trace_to_wig_cmd = "%s %s > %s" % (
        CONVERT_TRACE_INTO_WIGGLE_CMD,
        aggregated_min_trace_fname,
        aggregated_min_wig_fname
    )
    
    log_fp = open( log_fname_from_sample_name( sample_name, bs_id ), "a" )

    big_cmd = "\n" + ";\n".join(
        (
            statmap_cmd,
            build_marginal_cmd, # marginal wiggle

            # aggregated mins wiggle
            bootstrap_samples_cmd,
            aggregate_bootstraps_cmd,
            convert_aggregated_min_trace_to_wig_cmd,
        )
    ) + "\n"

    proc = subprocess.Popen(
        big_cmd, shell=True,
        stderr=subprocess.STDOUT, stdout=log_fp )
    
    return proc, log_fp


def main():
    '''
    Starts multiple sample mapping runs to process all samples given in a
    tab-delimited samples file
    (such as one produced by generate_multiple_samples_file.py)
    '''
    # parse arguments
    if len(sys.argv) != 2:
        print "Usage: ./multiple_sample_mapping.py samples.txt"; sys.exit(1)

    sample_filename = sys.argv[1]

    # add samples from file to process queue
    proc_queue = []
    for sample_name, bs_id, fnames in read_samples_file( sample_filename ):
        fastq_fname = fastq_fname_from_sample_name( sample_name, bs_id )

        # write out the trimmed, merged fastq's
        ofp = open( fastq_fname, "w" )
        strip_and_cat_files( fnames, ofp, 9 )
        ofp.close()

        proc_queue.append( (sample_name, bs_id) )

    # map all samples, MAX_NUM_PROC at a time
    running_procs = []
    while len(proc_queue) > 0 or len( running_procs ) > 0:
        for proc, log_fp in running_procs:
            proc.poll()
            if proc.returncode > 0:
                print proc.returncode, log_fp.name
        
        # empty out the finished processes
        procs_to_remove = [ i for i, (proc, log_fp) in enumerate(running_procs)
                            if proc.returncode != None  ]
        for proc_i in reversed(sorted( procs_to_remove )):
            running_procs[ proc_i ][1].close()
            del running_procs[ proc_i ]
        
        # if there is still space, add new processes
        for loop in xrange( MAX_NUM_PROC - len( running_procs ) ):
            # if there's no processes ready, break
            if len( proc_queue ) == 0: break
            running_procs.append( run_statmap( *(proc_queue.pop())) )
        
        # sleep a bit, so we're not wasting resources in a tight loop
        sleep( 1 )

if __name__ == "__main__": main()
