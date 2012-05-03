import psycopg2
import gzip
import subprocess
import os
from time import sleep

STATMAP_BIN = "~/statmap/bin/statmap"
BUILD_MARGINAL_BIN = "~/statmap/utilities/build_marginal_mappings_trace_from_mpd_read_db"
BUILD_AGGREGATES_CMD = "python ~/statmap/utilities/build_aggregates.py"
CONVERT_TRACE_INTO_WIGGLE_CMD = "~/statmap/bin/convert_trace_into_wiggle"
PARSE_INTO_PROMOTERS_CMD = "python /users/jbbrown/parse_cage_into_promoters.py"

#GENOME_LOC = "/media/scratch/genomes/drosophila/Manuel/genome.20.drosophila"
GENOME_LOC = "/media/scratch/genomes/drosophila/Manuel_latest/genome.20.drosophila"
NUM_THREAD_PER_PROC = 6
MAX_NUM_PROC = 4

def get_filenames_by_sample():
    conn = psycopg2.connect("host=eel dbname=labtrack user=nboley")
    cur = conn.cursor()
    query = """
    SELECT sample_name, array_agg(file) as files 
    FROM ( 
             SELECT sample_name, unnest(file) as file 
               FROM solexa_lane_view 
              WHERE sample_name ~ 'CAGE' 
    ) as foo 
    WHERE file IS NOT NULL 
    GROUP BY sample_name;
    """
    
    cur.execute( query )
    # remove spaces from the sample names, and return  the results
    return [ ( sn.replace(" ", "_"), fns ) for sn, fns in cur.fetchall() ]


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

def fastq_fname_from_sample_name( sample_name ):
    return sample_name + ".trimmed.fastq"

def output_dir_name_from_sample_name( sample_name ):
    return sample_name + "_mapped"

def log_fname_from_sample_name( sample_name ):
    return sample_name + ".run.log"

def run_statmap( sample_name ):
    statmap_cmd = " ".join([STATMAP_BIN,
                "-g %s" % GENOME_LOC,
                "-r %s" % fastq_fname_from_sample_name( sample_name ),
                "-t %i" % NUM_THREAD_PER_PROC,
                "-o %s" %  output_dir_name_from_sample_name( sample_name ),
                "-n 5", "-a a"
                ])
    
    # build marginal wiggle
    marginal_wig_fname = os.path.join( output_dir_name_from_sample_name( sample_name ), "marginal.wig" )

    build_marginal_cmd = "%s %s 1 > %s" % ( 
        BUILD_MARGINAL_BIN,
        output_dir_name_from_sample_name( sample_name ),
        marginal_wig_fname
    )

    # build minimum aggregate wiggle 
    aggregate_traces_cmd = "%s %s min" % (
        BUILD_AGGREGATES_CMD,
        output_dir_name_from_sample_name( sample_name ),
    )

    # write min trace to wiggle
    agg_min_trace_fname = "min_trace.bin.trace"
    agg_min_wig_fname = os.path.join( output_dir_name_from_sample_name( sample_name ), "agg_min.wig" )
    min_trace_to_wig_cmd = "%s %s > %s" % (
        CONVERT_TRACE_INTO_WIGGLE_CMD,
        agg_min_trace_fname,
        agg_min_wig_fname
    )
    
    # parse into promoters
    parse_into_proms_cmd = "%s %s 250 9 25 + > %s" % (
        PARSE_INTO_PROMOTERS_CMD, marginal_wig_fname, \
        os.path.join( output_dir_name_from_sample_name( sample_name ), "marginal.gff" )
    )
    
    log_fp = open( log_fname_from_sample_name( sample_name ), "a" )

    big_cmd = "\n" + ";\n".join( 
        (
            statmap_cmd,
            build_marginal_cmd,
            aggregate_traces_cmd,
            min_trace_to_wig_cmd,
            parse_into_proms_cmd
        )
    ) + "\n"

    proc = subprocess.Popen( big_cmd, shell=True, stderr=subprocess.STDOUT, stdout=log_fp )
    
    return proc, log_fp


proc_queue = []
for sample_name, fnames in get_filenames_by_sample():
    fastq_fname = fastq_fname_from_sample_name( sample_name )

    # write out the trimmed, merged fastq's
    ofp = open( fastq_fname, "w" )
    strip_and_cat_files( fnames, ofp, 9 )
    ofp.close()
    
    proc_queue.append( sample_name )


running_procs = []
while len(proc_queue) > 0 or len( running_procs ) > 0:
    for proc, log_fp in running_procs:
        proc.poll()
        if proc.returncode > 0:
            print proc.returncode, log_fp.name
    
    # empty out the finished processes
    procs_to_remove = [ i for i, (proc, log_fp) in enumerate( running_procs ) if proc.returncode != None  ]
    for proc_i in reversed(sorted( procs_to_remove )):
        running_procs[ proc_i ][1].close()
        del running_procs[ proc_i ]
    
    # if there is still space, add new processes
    for loop in xrange( MAX_NUM_PROC - len( running_procs ) ):
        # if there's no processes ready, break
        if len( proc_queue ) == 0: break
        running_procs.append( run_statmap( proc_queue.pop() ) )
    
    # sleep a bit, so we're not wasting resources in a tight loop
    sleep( 1 )

