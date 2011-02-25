from math import exp
import random
import numpy
from scipy import signal
#import rpy
import sys

BLOCK_SAMPLE_SIZE = 10000
NUM_SAMPLES = 100000
SEQ_LEN = 35
MIX = 0.3


# Command line parameters
input_dnase = sys.argv[1]
input_chrm = sys.argv[2]
chrm_length = int(sys.argv[3])
chrm_name = sys.argv[4]
frag_length = open(sys.argv[5])
output_prefix = sys.argv[6]


# save the simulated reads
OUTPUT_IP = open('./data/' + output_prefix + '_IP.fastq','w')
OUTPUT_CONTROL = open('./data/' + output_prefix + '_CONTROL.fastq','w')


# save the true wiggle so that we can get statistics on how accurately we are detecting binding sites of different strengths
TRUTH_WIG = open('./data/' + output_prefix + '_IP_truth.wig','w')
CONTROL_WIG = open('./data/' + output_prefix + '_CONTROL_truth.wig','w')

# save the exact location of each read so we can get statistics on how accurately each read is placed
TRUTH_TXT = open('./data/' + output_prefix + '_IP_truth_fastq.txt','w')
CONTROL_TXT = open('./data/' + output_prefix + '_CONTROL_truth_fastq.txt','w')


bps = ['A', 'C', 'G', 'T', 'N' ]
comp = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }
bp_index = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }

bcd_motif = [
    [ 35.000000,   35-11.090456,  35-9.002181,   35-5.832550, -1000],
    [-11.508077,  -26.296853,  -13.735887,   0.000000, -1000],
    [-11.533022,  -9.851050,   -14.876660,   0.000000, -1000],
    [ 0.000000,   -6.228376,   -9.521876,   -9.307275, -1000],
    [-6.842543,   -7.893736,    0.000000,   -7.614652, -1000],
    [-5.025504,   -4.851005,    0.000000,   -5.449993, -1000],
    [-3.162324,    0.000000,   -0.148652,   -1.038180, -1000],
    [-2.535381,   -0.893484,    0.000000,   -1.141034, -1000]
]

tf_conc = 1e-14


def parse_frag_length( fid ):
    data = []
    mx = 0
    for line in fid:
        d = line.strip().split('\t')
        data.append(d)
        d[0] = int(d[0])
        if d[0] > mx:
            mx = d[0]
    frag_dist = numpy.zeros( mx+1 )
    for d in data:
        frag_dist[d[0]] = float(d[1])
    return frag_dist

def score_binding_sites( region, motif ):
    scores = [0]*800
    for i in xrange( 800, len(region) - 800 - len(motif) ):
        score1 = sum( bcd_motif[j][bp_index[bp]] \
                      for j, bp in enumerate(region[i:(i+len(motif))]) )
        score2 = sum( bcd_motif[j][bp_index[comp[bp]]] \
                      for j, bp in enumerate(region[i:(i+len(motif))][::-1]) )
        score = max( score1, score2 )
        scores.append( tf_conc*exp( score  )/( 1 + tf_conc*exp( score  ) )  )
    scores.extend([0]*(800+len(motif)))
    return numpy.array( scores )

def build_chipseq_region( fname  ):
    fp = open( fname )
    genome = []
    for line in fp:
        if line.startswith( ">" ): continue
        genome.append( line.strip().upper() )
    return ''.join(genome)

def parse_wig( fname, chr_len ):
    fp = open( fname )
    density_fwd = numpy.zeros( chr_len )
    density_bkwd = numpy.zeros( chr_len )
    density = density_fwd
    for lin_num, line in enumerate(fp):
        if lin_num%1000000 == 0:
            print >> sys.stderr, "Parsed Through Line: ", lin_num
        if line.startswith('track'):
            density = density_bkwd
        elif line.startswith('variable'):
            continue
        else:
            data = line.strip().split("\t")
            loc, value = int(data[0]), float(data[1])
            density_fwd[ loc ] += value
    fp.close()
    return density_fwd, density_bkwd

def parse_bed( fname, chr_len ):
    # assumes data is unstranded.  Fine while doing DNase, should be generalized later.
    # although probably not necessary for this paper
    fp = open( fname )
    density = numpy.zeros( chr_len )
    for line in fp:
        data = line.strip().split('\t')
        chrm, start, end, name, signal = data[0:5]  
        if chrm == chrm_name:
            density[int(start):int(end)+1] += float(signal)
    fp.close()
    return density
    
def make_cum_dist( array ):
    sum = 0
    for index, entry in enumerate( array ):
        sum += entry
        array[index] = sum
    array = array/sum
    return array

def block_sample_density( density ):
    """Block sample frm the trace.

    """
    length = len( density )
    new_density = numpy.zeros( length  )
    for loop in xrange( 1, (length / BLOCK_SAMPLE_SIZE) ):
        start_loc = random.randint( 0, length - BLOCK_SAMPLE_SIZE  )
        new_density[((loop-1)*BLOCK_SAMPLE_SIZE):(loop*BLOCK_SAMPLE_SIZE)] \
           = density[start_loc:(start_loc+BLOCK_SAMPLE_SIZE)]
    return new_density/new_density.sum()

def block_sample_wiggle( wiggle ):
    """Block sample from the trace.
    coded a little diffeerently to make sure the length comes out right.
    """
    length = len( wiggle )
    new_density = []
    for ind in xrange(0, int(numpy.ceil((length + BLOCK_SAMPLE_SIZE) / BLOCK_SAMPLE_SIZE)) ):
        start_loc = random.randint( 0, length - BLOCK_SAMPLE_SIZE  )
        new_density.extend(wiggle[start_loc:(start_loc+BLOCK_SAMPLE_SIZE)])
    new_density = numpy.asarray(new_density[0:length])
    return new_density/new_density.sum()
    

def smooth_signal( input_sig, window):
    rv = signal.convolve( input_sig, window )[ 0:len(input_sig) ]
    return rv/rv.sum()
    
def sample_from_density_mixture( densities, mix_params ):
    assert( round( sum( mix_params ), 5 ) == 1.0 )
    # mixture params cumsum
    mpcs = numpy.array( mix_params ).cumsum()
    cumsums = [ density.cumsum() for density in densities ]
    
    for loop in xrange( NUM_SAMPLES ):
        index = mpcs.searchsorted( random.random() )
        bp_index = cumsums[index].searchsorted( random.random() )
        yield bp_index
    return
        

def sample_from_density( density ):
    cumsum = density.cumsum()
    bp_index = cumsum.searchsorted( random.random() )
    return bp_index



def print_wig( wig, outfile ):
    for bp,score in enumerate(wig):
        print >>outfile, chrm_name + '\t' + str(bp) + '\t' + str(bp) + '\t' + str(score)
    outfile.close()
    return


if __name__ == "__main__":
    #fwd_den, bkwd_den = parse_wig( input_bed, chrm_length  )
    fl_den = parse_frag_length( frag_length )
    fwd_den, rev_den = parse_wig( input_dnase, chrm_length  )
    region = build_chipseq_region( input_chrm )
    bind_site_scores = score_binding_sites( region, bcd_motif )
    #print sum( score > 0.2 for score in bind_site_scores )
    #sys.exit()
    #print fwd_den
    #sys.exit()
    bs_density = block_sample_wiggle( fwd_den )
    #sys.exit()
    #print len(bs_density), len(bind_site_scores)
    joint_dist = bs_density*bind_site_scores
    joint_dist = joint_dist/joint_dist.sum()
    
    print_wig(joint_dist, TRUTH_WIG)
    print_wig(bs_density, CONTROL_WIG)
    
    # NOTE FROM BEN:    The following smoothing makes it conceptually tricky to 
    #                   select the strandedness of reads.  Instead we'll do this
    #                   explicitely at read-generation time (see below).  This 
    #                   will also make it explicite to play with the uniformity
    #                   assumption regarding fragmentation.  
    # smooth all of the signals
    #joint_dist = smooth_signal( joint_dist, signal.gaussian( 200, 50 ) )
    #bs_density = smooth_signal( bs_density, signal.gaussian( 200, 50 ) )

    #bind_site_scores_cum = make_cum_dist( bind_site_scores )
    DEBUG_READS = False
    
    if DEBUG_READS: poss = []
    for loop, pos in enumerate( sample_from_density_mixture( [joint_dist, bs_density], [MIX, 1-MIX]  ) ):

        # select a fragment length from the input distribution:
        curr_frag_length = sample_from_density(fl_den)

        # select a position on the fragment for the center of the binding site under uniformity:
        # (note that the above assumption, uniformity, is something we can mess with to check stability and robustness)
        b_pos = int(numpy.random.randint(0,curr_frag_length))
        frag_start = pos - b_pos
        frag_end = pos + b_pos + curr_frag_length
        if frag_start < 0:
            continue
        if frag_end >= chrm_length:
            continue

        # select strand with 50/50 chance:
        # Note that python's [) intervals mean that an artificial binding-site width of 1 has been induced,
        # and the frag_dist input by the user is off by -1 (is 1 too short). 
        #frag = region[frag_start:frag_end]
        #frag_flip = ''.join( [ comp[bp] for bp in region[ frag_end:frag_start:-1 ] ] )  
        if numpy.random.rand() > 0.5:
            seq = region[ frag_start:(frag_start+SEQ_LEN) ]
            print >>TRUTH_TXT, '@' + str(loop) + '\t' + seq + '\t' + chrm_name + '\t' + str(frag_start) + '\t' + str(frag_start+SEQ_LEN) + '\t' + '+'
        else:
            seq = ''.join( [ comp[bp] for bp in region[ frag_end:(frag_end-SEQ_LEN):-1 ] ] )  
            print >>TRUTH_TXT, '@' + str(loop) + '\t' + seq + '\t' + chrm_name + '\t' + str(frag_end) + '\t' + str(frag_end-SEQ_LEN) + '\t' + '-'

        print >>OUTPUT_IP, "@%i" % loop
        print >>OUTPUT_IP, seq
        print >>OUTPUT_IP, "+%i" % loop  
        print >>OUTPUT_IP, "h"*len( seq )

        if DEBUG_READS: poss.append( pos )

    for loop, pos in enumerate( sample_from_density_mixture( [joint_dist, bs_density], [0.0, 1.0]  ) ):

        # select a fragment length from the input distribution:
        curr_frag_length = sample_from_density(fl_den)

        # select a position on the fragment for the center of the binding site under uniformity:
        # (note that the above assumption, uniformity, is something we can mess with to check stability and robustness)
        b_pos = int(numpy.random.randint(0,curr_frag_length))
        frag_start = pos - b_pos
        frag_end = pos + b_pos + curr_frag_length
        if frag_start < 0:
            continue
        if frag_end >= chrm_length:
            continue

        # select strand with 50/50 chance:
        # Note that python's [) intervals mean that an artificial binding-site width of 1 has been induced,
        # and the frag_dist input by the user is off by -1 (is 1 too short). 
        #frag = region[frag_start:frag_end]
        #frag_flip = ''.join( [ comp[bp] for bp in region[ frag_end:frag_start:-1 ] ] )  
        if numpy.random.rand() > 0.5:
            seq = region[ frag_start:(frag_start+SEQ_LEN) ]
            print >>CONTROL_TXT, '@' + str(loop) + '.control' + '\t' + seq + '\t' + chrm_name + '\t' + str(frag_start) + '\t' + str(frag_start+SEQ_LEN) + '\t' + '+'
        else:
            seq = ''.join( [ comp[bp] for bp in region[ frag_end:(frag_end-SEQ_LEN):-1 ] ] )  
            print >>CONTROL_TXT, '@' + str(loop) + '.control' + '\t' + seq + '\t' + chrm_name + '\t' + str(frag_end) + '\t' + str(frag_end-SEQ_LEN) + '\t' + '-'

        print >>OUTPUT_CONTROL, "@%i" % loop
        print >>OUTPUT_CONTROL, seq
        print >>OUTPUT_CONTROL, "+%i" % loop  
        print >>OUTPUT_CONTROL, "h"*len( seq )
        

    
    if DEBUG_READS:
        sampled_density = numpy.zeros( len( fwd_den ) )    
        for pos in poss:
            sampled_density[ pos:(pos+SEQ_LEN)  ] += 1
        sampled_density = sampled_density/sampled_density.sum()
    
        rpy.r.png( "test.png" )
        rpy.r.plot( bs_density[1:30000], type='l', main='', xlab='', ylab='', ylim=(0,max(bs_density.max(), joint_dist.max())) )
        rpy.r.lines( joint_dist[1:30000], type='l', col='blue' )
        rpy.r.lines( sampled_density[1:30000], type='l', col='red' )
        rpy.r.dev_off()
    
    
    # make_cum_dist( bkwd_den )
