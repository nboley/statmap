import subprocess
from sys import stdout, stderr
from itertools import izip
import re
from rpy import r
import numpy
from math import log

VERBOSE = True

STATMAP_PATH = "/users/nboley/statmap/trunk/bin/statmap"
MAPPED_READS_INTO_SAM_PATH = "/users/nboley/statmap/trunk/bin/mapped_reads_into_sam"
OUTPUT_DIR = "mapped_output_cage"
GENOME_FNAME = "/media/scratch/genomes/drosophila/BDGP_5/dros_bdgp_5.fa"
BINARY_GENOME_PATH = "/media/scratch/genomes/drosophila/BDGP_5/dros_bdgp_5.20.genome"
QUAL_SHIFT = 50

rc = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'  }

def est_prb_error_from_code( qual_code  ):
    return pow( 10, (float(ord(qual_code)) - 50)/-10 )

def rc_string( seq ):
    return "".join( rc[bp] for bp in reversed(seq) )

def load_genome( ):
    fp = open( GENOME_FNAME )
    curr_chr = None
    genome = {}
    for line in fp:
        if line.startswith( ">" ): 
            curr_chr = line.strip()[1:]
            genome[curr_chr] = []
            continue
        assert curr_chr != None
        genome[curr_chr].append( line.strip().upper() )
    fp.close()
    for chr_name in genome.keys():
        genome[chr_name] = "".join( genome[chr_name] )
    return genome


def map_data( fname ):
    """Map the data in the fastq file, fname.

    """
    
    
    call = "%s -g %s -s -r %s -s -o %s" \
         % ( STATMAP_PATH, BINARY_GENOME_PATH, fname, OUTPUT_DIR )
    if VERBOSE:
        print >> stdout, re.sub( "\s+", " ", call)    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        raise ValueError, "Mapping failed."

    pass

def build_sam(  ):
    fp = open( "output.sam", "w+" )
    call = "%s %s" % ( MAPPED_READS_INTO_SAM_PATH, OUTPUT_DIR )
    if VERBOSE:
        print >> stdout, re.sub( "\s+", " ", call)    
    ret_code = subprocess.call( call, shell=True, stdout=fp, stderr=stderr )
    if ret_code != 0:
        raise ValueError, "Building sam failed."
    fp.seek( 0 )
    return fp

def parse_out_unique_basepairs( fp, genome ):
    """Iterate through training basepairs and their design characteristics.
    
    Load the sam file in output_reads_fname, filter by the reads that are
    unique, and then yield basepairs along with error probabilities and counts.
    """
    fp.seek( 0 )
    
    for line_num, line in enumerate(fp):
        data = line.strip().split("\t")
        cond_prob = float(data[-1].split(":")[-1])
        if cond_prob < 1.0:
            continue
        genome_str = genome[ data[2] ][ int(data[3]):(int(data[3])+len(data[9]))  ]
        if data[1] == '16':
            genome_str = rc_string( genome_str )
        for index, (ref_bp, seq_bp, qual) in \
                enumerate( izip( genome_str, data[9], data[10] ) ):
            yield index+1, ref_bp, seq_bp, ord(qual)
    
    return

def analyze_errors( sam_fp, genome  ):
    # first, build the qual correction model
    qual_agg_count = dict()
    qual_agg_mismatch = dict()
    max_pos = -1
    for pos, ref, seq, qual in parse_out_unique_basepairs( sam_fp, genome ):
        # we also keep track of the maximum sequence position, to 
        # simplify the dictionary creation
        max_pos = max( max_pos, pos )
        
        # make sure that the qual key exists and, if not, add it
        if not qual_agg_count.has_key( qual ):
            qual_agg_count[qual] = 0
            qual_agg_mismatch[qual] = 0
    
        # increment the qual cnt
        qual_agg_count[ qual ] += 1
        
        # if necessary, increment the mismatch counts
        if ref != seq:
            qual_agg_mismatch[qual] += 1
    
    # plot the error estimates, just to have some idea of what we've actually done
    qual_keys = qual_agg_count.keys()
    qual_keys.sort()
    qual_est = {}
    for key in qual_keys:
        qual_est[key] = float(qual_agg_mismatch[key])/qual_agg_count[key]
        print key, float(qual_agg_mismatch[key]), qual_agg_count[key]
    
    data1, data2 = zip(*[ \
                          ( log(est_prb_error_from_code(chr(key))), log(float(val)) ) \
                          for key, val in  qual_est.iteritems() \
                   ] )
    r.png( "qual_plot.png"  )
    r.plot( data1, data2, \
            main='Pred vs Observed Error Rates', \
            xlab='Log Pred', ylab='Log Obs' )
    r.dev_off()

    fp = open("quals.csv", "w")
    for x, y in zip(data1, data2):
        print >> fp, "%e,%e" % ( x,y )
    fp.close()
    
    # now, fit the position model conditional on the error code estiamtes
    loc_agg_count = [0.0]*max_pos
    loc_agg_mismatch = [0.0]*max_pos
    loc_exp_mismatch = [0.0]*max_pos
    for pos, ref, seq, qual in parse_out_unique_basepairs( sam_fp, genome ):
        # increment the qual cnt
        loc_agg_count[ pos-1 ] += 1
        loc_exp_mismatch[ pos-1 ] += \
            float(qual_agg_mismatch[qual])/qual_agg_count[qual]
        if ref != seq:
            loc_agg_mismatch[ pos-1 ] += 1

    fp = open("pos.csv", "w")
    for loop in xrange(len(loc_agg_count)):
        print >> fp, ",".join( map( str, (loop+1, loc_agg_count[loop], loc_agg_mismatch[loop], loc_exp_mismatch[loop]) ) )
    fp.close()
    
    return 


#map_data( "./error_rate_tests/cage.dros.lane1.trimmed.first1M.fastq" )
genome = load_genome( )
sam_fp = build_sam()

analyze_errors( sam_fp, genome )

"""
fp = open("mismatches.pos.table", "w")
print >> fp, "Pos\tn_mis\ttotal"
for i in range(1,27):
    print >> fp, "%i\t%.0f\t%.0f" % (i, pos_agg_mismatch[ i ], pos_agg_count[ i ])
fp.close()

fp = open("qual.pos.table", "w")
print >> fp, "Qual\tn_mis\ttotal"
qual_keys = qual_agg_count.keys()
qual_keys.sort()
for qual in qual_keys:
    print >> fp, "%e\t%.0f\t%.0f" \
        % (qual, qual_agg_mismatch[ qual ], qual_agg_count[ qual ])
fp.close()
"""

# r analysis code
"""
data = read.table("mismatches.table", header=TRUE)
plot( data$Pos, data$n_mis/data$total )



"""
