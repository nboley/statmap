from math import exp
import random
import numpy
from scipy import signal
import rpy
import sys

BLOCK_SAMPLE_SIZE = 100000/100
NUM_SAMPLES = 100000
SEQ_LEN = 35

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

tf_conc = 1e-12

class frag_len_generator_t():
    def __init__( self ):
        self.sizes = range( 400, 801 )
        
        # sapply( 400:800, function(x)pnorm( x, 600, 100 ) )
        self.weights = [
            0.00000000, 0.02329547, 0.02385176, 0.02441919, 0.02499790, 0.02558806,
            0.02618984, 0.02680342, 0.02742895, 0.02806661, 0.02871656, 0.02937898,
            0.03005404, 0.03074191, 0.03144276, 0.03215677, 0.03288412, 0.03362497,
            0.03437950, 0.03514789, 0.03593032, 0.03672696, 0.03753798, 0.03836357,
            0.03920390, 0.04005916, 0.04092951, 0.04181514, 0.04271622, 0.04363294,
            0.04456546, 0.04551398, 0.04647866, 0.04745968, 0.04845723, 0.04947147,
            0.05050258, 0.05155075, 0.05261614, 0.05369893, 0.05479929, 0.05591740,
            0.05705343, 0.05820756, 0.05937994, 0.06057076, 0.06178018, 0.06300836,
            0.06425549, 0.06552171, 0.06680720, 0.06811212, 0.06943662, 0.07078088,
            0.07214504, 0.07352926, 0.07493370, 0.07635851, 0.07780384, 0.07926984,
            0.08075666, 0.08226444, 0.08379332, 0.08534345, 0.08691496, 0.08850799,
            0.09012267, 0.09175914, 0.09341751, 0.09509792, 0.09680048, 0.09852533,
            0.10027257, 0.10204232, 0.10383468, 0.10564977, 0.10748770, 0.10934855,
            0.11123244, 0.11313945, 0.11506967, 0.11702320, 0.11900011, 0.12100048,
            0.12302440, 0.12507194, 0.12714315, 0.12923811, 0.13135688, 0.13349951,
            0.13566606, 0.13785657, 0.14007109, 0.14230965, 0.14457230, 0.14685906,
            0.14916995, 0.15150500, 0.15386423, 0.15624765, 0.15865525, 0.16108706,
            0.16354306, 0.16602325, 0.16852761, 0.17105613, 0.17360878, 0.17618554,
            0.17878638, 0.18141125, 0.18406013, 0.18673294, 0.18942965, 0.19215020,
            0.19489452, 0.19766254, 0.20045419, 0.20326939, 0.20610805, 0.20897009,
            0.21185540, 0.21476388, 0.21769544, 0.22064995, 0.22362729, 0.22662735,
            0.22965000, 0.23269509, 0.23576250, 0.23885207, 0.24196365, 0.24509709,
            0.24825223, 0.25142890, 0.25462691, 0.25784611, 0.26108630, 0.26434729,
            0.26762889, 0.27093090, 0.27425312, 0.27759532, 0.28095731, 0.28433885,
            0.28773972, 0.29115969, 0.29459852, 0.29805597, 0.30153179, 0.30502573,
            0.30853754, 0.31206695, 0.31561370, 0.31917751, 0.32275811, 0.32635522,
            0.32996855, 0.33359782, 0.33724273, 0.34090297, 0.34457826, 0.34826827,
            0.35197271, 0.35569125, 0.35942357, 0.36316935, 0.36692826, 0.37069998,
            0.37448417, 0.37828048, 0.38208858, 0.38590812, 0.38973875, 0.39358013,
            0.39743189, 0.40129367, 0.40516513, 0.40904588, 0.41293558, 0.41683384,
            0.42074029, 0.42465457, 0.42857628, 0.43250507, 0.43644054, 0.44038231,
            0.44433000, 0.44828321, 0.45224157, 0.45620469, 0.46017216, 0.46414361,
            0.46811863, 0.47209683, 0.47607782, 0.48006119, 0.48404656, 0.48803353,
            0.49202169, 0.49601064, 0.50000000, 0.50398936, 0.50797831, 0.51196647,
            0.51595344, 0.51993881, 0.52392218, 0.52790317, 0.53188137, 0.53585639,
            0.53982784, 0.54379531, 0.54775843, 0.55171679, 0.55567000, 0.55961769,
            0.56355946, 0.56749493, 0.57142372, 0.57534543, 0.57925971, 0.58316616,
            0.58706442, 0.59095412, 0.59483487, 0.59870633, 0.60256811, 0.60641987,
            0.61026125, 0.61409188, 0.61791142, 0.62171952, 0.62551583, 0.62930002,
            0.63307174, 0.63683065, 0.64057643, 0.64430875, 0.64802729, 0.65173173,
            0.65542174, 0.65909703, 0.66275727, 0.66640218, 0.67003145, 0.67364478,
            0.67724189, 0.68082249, 0.68438630, 0.68793305, 0.69146246, 0.69497427,
            0.69846821, 0.70194403, 0.70540148, 0.70884031, 0.71226028, 0.71566115,
            0.71904269, 0.72240468, 0.72574688, 0.72906910, 0.73237111, 0.73565271,
            0.73891370, 0.74215389, 0.74537309, 0.74857110, 0.75174777, 0.75490291,
            0.75803635, 0.76114793, 0.76423750, 0.76730491, 0.77035000, 0.77337265,
            0.77637271, 0.77935005, 0.78230456, 0.78523612, 0.78814460, 0.79102991,
            0.79389195, 0.79673061, 0.79954581, 0.80233746, 0.80510548, 0.80784980,
            0.81057035, 0.81326706, 0.81593987, 0.81858875, 0.82121362, 0.82381446,
            0.82639122, 0.82894387, 0.83147239, 0.83397675, 0.83645694, 0.83891294,
            0.84134475, 0.84375235, 0.84613577, 0.84849500, 0.85083005, 0.85314094,
            0.85542770, 0.85769035, 0.85992891, 0.86214343, 0.86433394, 0.86650049,
            0.86864312, 0.87076189, 0.87285685, 0.87492806, 0.87697560, 0.87899952,
            0.88099989, 0.88297680, 0.88493033, 0.88686055, 0.88876756, 0.89065145,
            0.89251230, 0.89435023, 0.89616532, 0.89795768, 0.89972743, 0.90147467,
            0.90319952, 1.00000000
        ]
    
    def sample_frag_len( self ):
        return self.sizes[ bisect.bisect_left( self.weights, random.random() ) ]

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
        
    # pick an index at random

if __name__ == "__main__":
    fwd_den, bkwd_den = parse_wig( "./data/random_dnase_sample_chr2L_short.wig", 49950  )
    
    region = build_chipseq_region( "./data/chr2L_short.fa" )
    bind_site_scores = score_binding_sites( region, bcd_motif )
    bs_density = block_sample_density( fwd_den )

    joint_dist = bs_density*bind_site_scores
    joint_dist = joint_dist/joint_dist.sum()
    
    # smooth all of the signals
    joint_dist = smooth_signal( joint_dist, signal.gaussian( 200, 50 ) )
    bs_density = smooth_signal( bs_density, signal.gaussian( 200, 50 ) )

    #bind_site_scores_cum = make_cum_dist( bind_site_scores )
    DEBUG_READS = True
    
    if DEBUG_READS: poss = []
    for loop, pos in enumerate( \
        sample_from_density_mixture( [joint_dist, bs_density], [0.1, 0.9]  ) ):
        seq = region[ pos:(pos+SEQ_LEN)  ]
        print "@%i" % loop
        print seq
        print "+"
        print "h"*len( seq )

        if DEBUG_READS: poss.append( pos )
        

    
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
