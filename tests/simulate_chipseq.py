import random
import array
from math import exp
import bisect
import rpy
import numpy
import sys
import os
import fnmatch
import re
import subprocess

import tests as sc # for simulation code

NUM_READS = 1000
NUM_SAMPLES = 100

bps = ['A', 'C', 'G', 'T' ]
comp = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }
bp_index = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }


bcd_motif = [
    [  0.0,  0.0, 0.0, 12.0 ],
    [ 12.0,  0.0, 0.0,  0.0 ],
    [ 12.0,  0.0, 0.0,  0.0 ],
    [  0.0,  0.0, 0.0, 12.0 ],
    [  0.0,  5.0, 0.0,  0.0 ],
    [  0.0,  5.0, 0.0,  0.0 ]
]

bcd_motif = [
    [ 35.000000,   35-11.090456,  35-9.002181,   35-5.832550],
    [-11.508077,  -26.296853,  -13.735887,   0.000000],
    [-11.533022,  -9.851050,   -14.876660,   0.000000],
    [ 0.000000,   -6.228376,   -9.521876,   -9.307275],
    [-6.842543,   -7.893736,    0.000000,   -7.614652],
    [-5.025504,   -4.851005,    0.000000,   -5.449993],
    [-3.162324,    0.000000,   -0.148652,   -1.038180],
    [-2.535381,   -0.893484,    0.000000,   -1.141034]
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


# build a random chr from which to simulate chipseq
def build_chipseq_region(  ):
    fp = open( "./data/bcd_region.fasta" )
    genome = []
    for line in fp:
        if line.startswith( ">" ): continue
        genome.append( line.strip().upper() )
    return ''.join(genome)

def score_binding_sites( region, motif ):
    scores = [0]*800
    for i in xrange( 800, len(region) - 800 - len(motif) ):
        score1 = sum( bcd_motif[j][bp_index[bp]] for j, bp in enumerate(region[i:(i+len(motif))]) )
        score2 = sum( bcd_motif[j][bp_index[comp[bp]]] for j, bp in enumerate(region[i:(i+len(motif))][::-1]) )
        scores.append( max( score1, score2 ) )
    scores.extend([0]*800)
    return array.array( 'f', scores )

def assign_bind_prbs( scores ):
    return array.array( 'f', [ tf_conc*exp( score  )/( 1 + tf_conc*exp( score  ) ) for score in scores] )

def build_cum_array( items ):
    new_list = []
    curr_sum = 0
    for item in items:
        curr_sum += item
        new_list.append( curr_sum )
    for i in xrange( len( new_list ) ):
        new_list[i] /= curr_sum
    new_list.append( 1.0 )
    return array.array( 'f', new_list )

def build_fragments( region, bind_prbs, num ):
    motif = bcd_motif        
    cum_array = build_cum_array( bind_prbs )
    fl_gen = frag_len_generator_t()
    
    fragments = []

    for i in xrange( num ):
        # choose a random binding site
        bs_index = bisect.bisect_left( cum_array, random.random() )
        # determine the fragment length for this read
        fl = fl_gen.sample_frag_len()
        # determine where, on the fragment, the bs is
        # ( we say that it is uniform, except for the edges )
        bind_loc = random.randint( 0, fl - len(motif)  )
        fl_start = bs_index - bind_loc
        fragments.append( ('chr2L', fl_start, fl_start+fl) )

    """
    """
    
    return fragments

def parse_wig( fname, genome ):
    fp = open( fname )
    density = numpy.zeros(len(genome.values()[0]))
    for line in fp:
        if line.startswith('track') or line.startswith('variable'):
            continue
        data = line.strip().split("\t")
        loc, value = int(data[0]), float(data[1])
        density[loc] += value
    fp.close()
    return density

def test_chipseq_region( num_mutations, wig_fname = 'tmp.wig', iterative=True ):
    region = build_chipseq_region( )
    bind_site_scores = score_binding_sites( region, bcd_motif )
    bind_prbs = assign_bind_prbs( bind_site_scores )
    fragments = build_fragments( region, bind_prbs, NUM_READS )
    
    genome = { 'chr2L': region }
    
    reads = sc.build_reads_from_fragments( genome, fragments, paired_end=True )
    
    genome_of = open("tmp.genome", "w")
    # write a second genome
    chr2L = genome['chr2L']
    mutated_chr2L = array.array( 'c', chr2L )

    # mutated num_mutations random indexes
    rand_indexes = random.sample( xrange(len(mutated_chr2L )), num_mutations )
    for rand_index in rand_indexes:
        curr_bp = mutated_chr2L[ rand_index ]
        valid_bps = [ bp for bp in bps if bp != curr_bp ]
        mutated_chr2L[rand_index] = random.choice( valid_bps )
    
    genome['chr2L'] = chr2L + mutated_chr2L.tostring()
    sc.write_genome_to_fasta( genome, genome_of, 1 )
    genome_of.close()
    
    reads_of_1 = open('tmp.1.fastq', "w");
    reads_of_2 = open('tmp.2.fastq', "w");
    fastq = sc.build_paired_end_fastq_from_seqs( reads, reads_of_1, reads_of_2 )
    reads_of_1.close()
    reads_of_2.close()

    call = "%s -g tmp.genome -1 tmp.1.fastq -2 tmp.2.fastq \
                             -p %.2f -m %.2f -n %i \
                             " % ( sc.STATMAP_PATH, -10, 2, NUM_SAMPLES )
    if iterative:
        call += " -a i"
        
    print re.sub( "\s+", " ", call)

    ret_code = subprocess.call( call, shell=True )
    # ret_code = ( os.system( call ) >> 8 )
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
        
    ## Plot everything

    return

    true_density = numpy.zeros( len( region ) )
    for chr, start, stop in fragments:
        true_density[start:stop] += 1
    true_density = true_density/true_density.sum()
    
    density = parse_wig( wig_fname, genome )
    density = density/density.sum()
    
    rpy.r.plot( density, col='green', ylim=(0,true_density.max()), \
                type='l', main='', ylab='', xlab='' )
    rpy.r.points( true_density, type='l' )
    rpy.r.points( density.max()*numpy.array(bind_prbs), type='l', col='red' )
    raw_input()

def build_all_wiggles( mutations, do_marginal=True ):
    all_mutations = [ ]
    all_mutations.extend( mutations )
    for mr in all_mutations:
        if do_marginal:
            test_chipseq_region( mr, "relaxed_%i_marginal.wig" % mr, iterative=False  )
        test_chipseq_region( mr, "relaxed_%i_relaxed.wig" % mr, iterative=True  )

def plot_all_wiggles( mutations ):
    # build and plot the true density
    region = build_chipseq_region( )
    genome = { 'chr2L': region + region }
    bind_site_scores = score_binding_sites( region, bcd_motif )
    bind_prbs = assign_bind_prbs( bind_site_scores )
    fragments = build_fragments( region, bind_prbs, NUM_READS )

    true_density = numpy.zeros( 2*len( region ) )
    for chr, start, stop in fragments:
        true_density[start:stop] += 1
    true_density = true_density/true_density.sum()
    
    rpy.r.png( 'chipseq_sim_output.png', width=1900, height=750 )

    rpy.r.par(mfrow=(2,1))

    # plot the real density, and an impossible one
    rpy.r.plot( true_density, col='black', \
                ylim=(0,1.2*true_density.max()), \
                type='l', main='', 
                ylab='Density', xlab='Bps from Region Start' )
    

    rpy.r.points( true_density, col='black', type='l' )
    
    wig_fname = "relaxed_%i_marginal.wig" % 0
    density = parse_wig( wig_fname, genome )
    density = density/density.sum()
    rpy.r.points( density, type='l', col='red', main='', xlab='', ylab='', lty=4 )
    
    wig_fname = "relaxed_%i_relaxed.wig" % 0
    density = parse_wig( wig_fname, genome )
    density = density/density.sum()
    rpy.r.points( density, type='l', col='red' )

    rpy.r.points( true_density.max()*numpy.array(bind_prbs), type='l', col='green' )

    if len(mutations) == 0:
        rpy.r.dev_off()
        return
    
    # plot the various mutation rates
    wig_fname = "relaxed_%i_marginal.wig" % mutations[0]
    density = parse_wig( wig_fname, genome )
    density = density/density.sum()
    rpy.r.plot( density, type='l', ylim=(0,1.2*true_density.max()), col='orange', ylab='Density', xlab='Bps from Region Start', lty=4 )
    
    wig_fname = "relaxed_%i_relaxed.wig" % mutations[0]
    density = parse_wig( wig_fname, genome )
    density = density/density.sum()
    rpy.r.points( density, type='l', col='orange' )
    
    for mr, color in zip( mutations[1:], \
                          ( 'blue', 'red', 'green', 'black') ):
        wig_fname = "relaxed_%i_marginal.wig" % mr
        print wig_fname
        density = parse_wig( wig_fname, genome )
        density = density/density.sum()
        rpy.r.points( density, type='l', col=color, main='', xlab='', ylab='', lty=4 )
        
        wig_fname = "relaxed_%i_relaxed.wig" % mr
        density = parse_wig( wig_fname, genome )
        density = density/density.sum()
        rpy.r.points( density, type='l', col=color )
    
    rpy.r.dev_off()

def plot_wig_bounds( dir, png_fname):
    region = build_chipseq_region( )
    genome = { 'chr2L': region + region }

    curr_dir = os.getcwd()
    rpy.r.png( os.path.join(curr_dir, png_fname), width=1900, height=750 )
    
    density = parse_wig( "./statmap_output/max_trace.wig", genome )
    density_max = density.max()
    density = density/density_max
    rpy.r.plot( density, type='l', col='black', main='', xlab='', ylab='', lty=4, ylim=(0, 1.2) )

    fnames = []
    os.chdir(dir)
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, '*.wig'):
            fnames.append( file )

    for fname in fnames:
        density = parse_wig( fname, genome )/density_max
        rpy.r.points( density, type='l', col='black', main='', xlab='', ylab='', lty=4 )

    os.chdir(curr_dir)
    
    density = parse_wig( "./statmap_output/min_trace.wig", genome )/density_max
    rpy.r.points( density, type='l', col='red', main='', xlab='', ylab='', lty=4 )

    density = parse_wig( "./statmap_output/relaxed_mapping.wig", genome )/density_max
    rpy.r.points( density, type='l', col='green', main='', xlab='', ylab='', lty=4 )
        
    rpy.r.dev_off()



if __name__ == '__main__':
    mutations = [2,] # [ 1, 10, 25, 100, 1000  ]
    build_all_wiggles( mutations, False )
    plot_wig_bounds( "./statmap_output/samples/", "relaxed_samples.png")
    plot_wig_bounds( "./statmap_output/starting_samples/", "starting_samples.png")
    if sc.CLEANUP:
        subprocess.call( "rm tmp.*", shell=True )
        subprocess.call( "rm ./statmap_output/ -rf", shell=True )
    sys.exit( -1 )

    mutations = [] # [ 1, 10, 25, 100, 1000  ]
    build_all_wiggles( mutations, False )
    plot_all_wiggles( mutations )
    if sc.CLEANUP:
        subprocess.call( "rm tmp.*", shell=True )
        subprocess.call( "rm ./statmap_output/ -rf", shell=True )






