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
import gzip

sys.path.insert(0, "../python_lib/" )
import trace

import tests as sc # for simulation code

NUM_READS = 1000
NUM_SAMPLES = 5

bps = ['A', 'C', 'G', 'T', 'N' ]
comp = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }
bp_index = { 'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4 }

def build_cage_region( length = 5000 ):
    pass

def build_random_cage_reads( num_mutations, DIRTY=True, are_paired_end=True, polymorphic=True ):
    region = build_chipseq_region( )
    
    # simulate the ip reads
    bind_site_scores = score_binding_sites( region, bcd_motif )
    bind_prbs = assign_bind_prbs( bind_site_scores )
    ip_fragments = build_fragments( region, bind_prbs, NUM_READS )
    
    # simulate the negative control reads
    # zero the boundaries to prevent framgnets from overlapping
    zeroed_bndry_size = 650
    uniform = array.array( 'f', [1.0]*len(bind_prbs) )
    uniform[:zeroed_bndry_size] = array.array( 'f', [0.0]*zeroed_bndry_size )
    uniform[-zeroed_bndry_size:] = array.array( 'f', [0.0]*zeroed_bndry_size )
    nc_fragments = build_fragments( region, uniform, NUM_READS )
    
    # write the binding site probs to a wiggle
    # so that we can plot them
    fp = open( "binding_site_occupancies.wig", "w" )
    fp.write("track type=wiggle_0 name=binding_site_occupancy\n")
    fp.write("variableStep chrom=chr2L_paternal\n")
    for pos, value in enumerate(bind_prbs):
        if value > 0:
            fp.write( "%i\t%e\n" % (pos+1, value) )
    fp.close()
    
    density = numpy.zeros(len(region))
    for chr, start, stop in ip_fragments:
        density[start:stop] += 1.0/(stop-start)
    fp = open( "true_read_coverage.wig", "w" )
    fp.write("track type=wiggle_0 name=fwd_true_read_coverage\n")
    fp.write("variableStep chrom=chr2L_paternal\n")
    for pos, value in enumerate(density):
        if value > 0:
            fp.write( "%i\t%e\n" % (pos+1, value/2) )
    fp.write("track type=wiggle_0 name=bkwd_true_read_coverage\n")
    fp.write("variableStep chrom=chr2L_paternal\n")
    for pos, value in enumerate(density):
        if value > 0:
            fp.write( "%i\t%e\n" % (pos+1, value/2) )
    fp.close()
    
    genome = { 'chr2L_paternal': region }
    
    if DIRTY:
        sample_file = gzip.open( './data/cage_error_strs.fastq.gz' )
        error_strs = sc.get_error_strs_from_fastq( sample_file )
        sample_file.close()

        ip_reads = sc.build_reads_from_fragments( 
            genome, ip_fragments, rev_comp=True, paired_end=are_paired_end )
        ip_mutated_reads = sc.mutate_reads( ip_reads, error_strs )
        
        nc_reads = sc.build_reads_from_fragments( 
            genome, nc_fragments, rev_comp=True, paired_end=are_paired_end )
        nc_mutated_reads = sc.mutate_reads( nc_reads, error_strs )
    else:
        ip_reads = sc.build_reads_from_fragments( 
            genome, ip_fragments, are_paired_end=paired_end )
        nc_reads = sc.build_reads_from_fragments( 
            genome, nc_fragments, are_paired_end=paired_end )
    
    genome_of = open("tmp.genome", "w")
    # write a second genome
    chr2L = genome['chr2L_paternal']
    mutated_chr2L = array.array( 'c', chr2L )
    
    # mutate num_mutations random indexes
    rand_indexes = random.sample( xrange(len(mutated_chr2L )), num_mutations )
    for rand_index in rand_indexes:
        curr_bp = mutated_chr2L[ rand_index ]
        valid_bps = [ bp for bp in bps if bp != curr_bp ]
        mutated_chr2L[rand_index] = random.choice( valid_bps )

    if polymorphic:
        genome['chr2L_maternal'] = mutated_chr2L.tostring()
    
    sc.write_genome_to_fasta( genome, genome_of )
    genome_of.close()
    
    reads_of_1 = open('tmp.1.fastq', "w");
    reads_of_2 = open('tmp.2.fastq', "w");

    reads_nc_of_1 = open('tmp.nc.1.fastq', "w");
    reads_nc_of_2 = open('tmp.nc.2.fastq', "w");

    
    if not paired_end and not DIRTY:
        sc.build_single_end_fastq_from_seqs( ip_reads, reads_of_1 )
        sc.build_single_end_fastq_from_seqs( nc_reads, reads_nc_of_1 )
    elif paired_end and not DIRTY:
        sc.build_paired_end_fastq_from_seqs( ip_reads, reads_of_1, reads_of_2 )
        sc.build_paired_end_fastq_from_seqs( nc_reads, reads_nc_of_1, reads_nc_of_2 )
    elif not paired_end and DIRTY:
        sc.build_single_end_fastq_from_mutated_reads( ip_mutated_reads, reads_of_1 )
        sc.build_single_end_fastq_from_mutated_reads( nc_mutated_reads, reads_nc_of_1 )
    elif paired_end and DIRTY:        
        sc.build_paired_end_fastq_from_mutated_reads( ip_mutated_reads,reads_of_1,reads_of_2 )
        sc.build_paired_end_fastq_from_mutated_reads( nc_mutated_reads,reads_nc_of_1,reads_nc_of_2 )
    
    reads_of_1.close()
    reads_of_2.close()
    
    return rand_indexes

def map_with_statmap( iterative=True, paired_end=True ):
    # build the index
    call = "%s tmp.genome 20 tmp.genome.bin" % sc.BUILD_INDEX_PATH
    print re.sub( "\s+", " ", call)
    ret_code = subprocess.call( call, shell=True )    
    if ret_code != 0:
        print "TEST FAILED - build_index call returned error code ", ret_code
        sys.exit( -1 )    
    
    if paired_end:
        call = "%s -g tmp.genome.bin -1 tmp.1.fastq -2 tmp.2.fastq \
                                     -3 tmp.nc.1.fastq -4 tmp.nc.2.fastq \
                                     -o smo_chipseq_sim -q 0 \
                                     -n %i -f ./data/fl_dist.txt\
                                 " % ( sc.STATMAP_PATH, NUM_SAMPLES )
    else:
        call = "%s -g tmp.genome.bin -r tmp.1.fastq -c tmp.nc.1.fastq \
                                     -o smo_chipseq_sim -q 0 \
                                     -n %i -f ./data/fl_dist.txt\
                                 " % ( sc.STATMAP_PATH, NUM_SAMPLES )
    
    if iterative:
        call += " -a i"
        
    # run statmap
    print re.sub( "\s+", " ", call)    
    ret_code = subprocess.call( call, shell=True )    
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
    
    # run the aggregation generation code
    call = [ "../utilities/build_aggregates.py", "smo_chipseq_sim" ]
    print >> sc.stdout, " ".join( call )
    ret_code = subprocess.call( call  )
    if ret_code != 0:
        print "TEST FAILED - aggregation call returned error code ", ret_code
        sys.exit( -1 )
        
    return

def call_peaks( ):
    ret_code = subprocess.call( "%s ./smo_chipseq_sim/" % sc.CALL_PEAKS_PATH, shell=True )    
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        print "%s ./smo_chipseq_sim/" % sc.CALL_PEAKS_PATH
        sys.exit( -1 )
    

def map_with_bowtie( paired_end = True ):
    # build the index
    # pie out to null to ignore output
    cmd = "bowtie-build -f tmp.genome bowtie_index/tmp.ebwt > /dev/null"
    subprocess.call( cmd, shell=True )
    # map the reads with bowtie
    if paired_end:
        cmd = "bowtie -a --tryhard -X 2500 --fr bowtie_index/tmp.ebwt \
               -1 tmp.1.fastq -2 tmp.2.fastq mapped_reads.bwtout"
    else:
        cmd = "bowtie -a --tryhard -X 2500 --fr bowtie_index/tmp.ebwt \
               tmp.1.fastq mapped_reads.bwtout"        
    subprocess.call( cmd, shell=True )
    
def plot_bootstrap_bounds( png_fname, paired_end, mut_indexes=[], polymorphic=True ):
    region = build_chipseq_region( )
    if polymorphic:
        genome = { 'chr2L_paternal': region, 'chr2L_maternal': region }
    else: 
        genome = { 'chr2L_paternal': region }
            
    # set up the plotting environment
    curr_dir = os.getcwd()
    rpy.r.png( os.path.join(curr_dir, png_fname), width=7.0, height=3.5*len(genome), units='in', res=300 )
    rpy.r("par(cex=0.47, mai=c(0.5,0.4,0.4,0.2), lwd=0.5, mfrow=c(%i,1))" % len(genome) )
    
    density_max = 1.0 # max( density[0].max(), density[1].max() )
    if polymorphic:
        rpy.r.plot( (), main='Paternal', xlab='', ylab='', lty=1, xlim=(0,5000), ylim=(-0.65, 0.65) )
        rpy.r.plot( (), main='Maternal', xlab='', ylab='', lty=1, xlim=(0,5000), ylim=(-0.65, 0.65) )
    else:
        rpy.r.plot( (), main='', xlab='', ylab='', lty=1, xlim=(0,5000), ylim=(-0.65, 0.65) )

    def plot_wiggle( wiggle_density, color, lty=1, lwd=0.5, norm_factor = 1.0  ):
        for index, key in enumerate(genome.keys()):
            rpy.r("par(mfg=c(%i,1))" % (index+1))
            if len( wiggle_density ) > 0 and wiggle_density[0].has_key( key ):
                rpy.r.points( wiggle_density[0][key]*norm_factor, type='l', col=color, main='', xlab='', ylab='', lty=lty, lwd=lwd )
            if len( wiggle_density ) > 1 and wiggle_density[1].has_key( key ):
                rpy.r.points( -1*wiggle_density[1][key]*norm_factor, type='l', col=color, main='', xlab='', ylab='', lty=lty, lwd=lwd )
            if len( wiggle_density ) > 2 and wiggle_density[2].has_key( key ):
                rpy.r.points( wiggle_density[2][key]*norm_factor, type='l', col=color, main='', xlab='', ylab='', lty=3, lwd=lwd/2.0 )
            if len( wiggle_density ) > 3 and wiggle_density[3].has_key( key ):
                rpy.r.points( -1*wiggle_density[3][key]*norm_factor, type='l', col=color, main='', xlab='', ylab='', lty=3, lwd=lwd/2.0 )

    def plot_trace( density, color, lty=1, lwd=0.5, norm_factor = 1.0  ):
        # plot the IP traces
        for track_name in density.keys():
            for index, chr_name in enumerate(genome.keys()):
                rpy.r("par(mfg=c(%i,1))" % (index+1))
                # change the line type for NC
                curr_lwd = lwd
                curr_lty = lty
                if track_name.startswith("NC"):
                    curr_lwd /= 2
                    curr_lty = 3
                # change the multiplicative factor for neg vs pos stranded data
                mult_factor = -1 if track_name.find( "bkwd" ) != -1 else 1
                rpy.r.points( mult_factor*density[track_name][chr_name]*norm_factor, 
                              type='l', col=color, main='', xlab='', ylab='', 
                              lty=lty, lwd=lwd )
    
    def plot_traces( dir, color, lty=1, lwd=0.5, filter=""  ):
        fnames = []
        os.chdir(dir)
        for file in os.listdir("./"):
            if fnmatch.fnmatch(file, '*%s.bin.trace' % filter):
                fnames.append( file )
                
        for fname in fnames:
            density = parse_trace( fname )
            plot_trace( density, color )
        os.chdir(curr_dir)

    plot_traces( "./smo_chipseq_sim/samples/", 'gray', lty=3, lwd=1, filter=".nc" )
    plot_traces( "./smo_chipseq_sim/samples/", 'black', lty=1, lwd=0.5, filter=".ip" )

    plot_traces( "./smo_chipseq_sim/bootstrap_samples/min_traces/", 'green', lty=1, lwd=0.5, filter=".ip" )
    plot_traces( "./smo_chipseq_sim/bootstrap_samples/max_traces/", 'orange', lty=1, lwd=0.5, filter=".ip" )
    
    rpy.r("par(mfg=c(2,1))")
    true_density = parse_wig( "true_read_coverage.wig", genome )    
    plot_wiggle( true_density, 'black', lty=3, lwd=3 )
    
    for index in xrange(len(genome)):
        rpy.r("par(mfg=c(%i,1))" % (index+1))
        for bp_index in mut_indexes:
            rpy.r.abline( v=bp_index, col='red', lty=2  )

    # parse bowtie out
    # density = parse_bwtout( "mapped_reads.bwtout", genome, paired_end )
    # plot_wiggle( density, 'dark green', lty=1, lwd=2 )
    
    for index in xrange(len(genome)):
        rpy.r("par(mfg=c(%i,1))" % (index+1))
        true_density = parse_wig( "binding_site_occupancies.wig", genome )
        for index, value in enumerate(true_density[0].values()[0]):
            if value > 0.50:
                rpy.r.abline( v=index, col='black', lty=2, lwd=1  )
    
    
    # marginal mapping
    """
    density = parse_trace( "./smo_chipseq_sim/marginal_mappings_fwd.wig", genome )
    plot_wiggle( density, 'blue' )

    density = parse_trace( "./smo_chipseq_sim/marginal_mappings_bkwd.wig", genome )
    plot_wiggle( density, 'blue' )
    """
    
    density = parse_trace( "./smo_chipseq_sim/max_trace.ip.bin.trace" )
    plot_trace( density, 'blue' )
    
    density = parse_trace( "./smo_chipseq_sim/min_trace.ip.bin.trace" )
    plot_trace( density, 'red' )

    # plot the called peak p-values

    nf = 0.4
    #density = parse_wig( "./smo_chipseq_sim/called_peaks/peaks.wig", genome )
    #plot_wiggle( density, 'red', lwd=2, norm_factor=nf )
    density = parse_trace( "./smo_chipseq_sim/called_peaks/sample1.bin.trace" )
    plot_trace( density, 'red', lwd=2, norm_factor=nf )
    for index in xrange(len(genome)):
        rpy.r("par(mfg=c(%i,1))" % (index+1))
        rpy.r.abline( h=nf*0.95, col='red', lwd=1, lty=2  )
        rpy.r.abline( h=-0.95*nf, col='red', lwd=1, lty=2  )

    
    """
    density = parse_trace( "./smo_chipseq_sim/relaxed_mapping.wig", genome )
    if len( density ) > 0:
        rpy.r.points( density[0]/density_max, type='l', col='green', lwd=1.5, main='', xlab='', ylab='', lty=1 )
    if len( density ) > 1:
        rpy.r.points( -1*density[1]/density_max, type='l', col='green', lwd=1.5, main='', xlab='', ylab='', lty=1 )
    """
        
    for index in xrange(len(genome)):
        rpy.r("par(mfg=c(%i,1))" % (index+1))
        for bp_index in mut_indexes:
            rpy.r.abline( v=bp_index, col='red', lty=2  )

    """
    # x=7700, y=0.5,
    rpy.r(""legend( x=50, y=0.49,
             legend=c("Statmap Upper Bound", "Statmap Lower Bound", 
                      "Statmap Local Maxima", "Statmap Bootstrap Upper Bounds", 
                      "Statmap Bootstrap Lower Bounds",
                      "Bowtie", "True Read Coverage", "Maternal SNPs", "Bind Sites w/ > 0.25 Average Occupancy"),
             col=c("Blue", "Red", "Black", "Purple", "Orange", "Green", "Black", "Red", "Black"),
             lwd=c(0.5,0.5,0.5,0.5,0.5,2,2,0.5,0.5), lty=c(1,1,1,1,1,3,3,2,2) )"" )
    """

    rpy.r.dev_off()


if __name__ == '__main__':
    paired_end=False
    NUM_MUTS = 3
    
    if False:
        test_cage_region( NUM_MUTS, "relaxed_%i_marginal.wig" % NUM_MUTS, iterative=False )
        if sc.CLEANUP:
            #subprocess.call( "rm tmp.*", shell=True )
            #subprocess.call( "rm ./smo_chipseq_sim/ -rf", shell=True )
            pass

    if True:
        mut_indexes = build_random_chipseq_reads( NUM_MUTS, are_paired_end=paired_end, polymorphic=True )
        #print mut_indexes
        map_with_bowtie( paired_end )
        map_with_statmap( paired_end=paired_end )
        call_peaks()
        #mut_indexes = [ 1327, 3755, 261 ]
        plot_bootstrap_bounds( "bootstrap_bnds.png", paired_end, mut_indexes, polymorphic=True )
        # plot_wig_bounds( "./smo_chipseq_sim/samples/", "relaxed_samples.png")
        # plot_wig_bounds( "./smo_chipseq_sim/starting_samples/", "starting_samples.png")
        if False and sc.CLEANUP:
            subprocess.call( "rm tmp.*", shell=True )
            # subprocess.call( "rm ./smo_chipseq_sim/ -rf", shell=True )

    # BROKEN
    """
    mutations = [ 1, 3, 10, 25, 100  ]
    build_all_wiggles( mutations, False )
    plot_all_wiggles( mutations )
    if sc.CLEANUP:
        subprocess.call( "rm tmp.*", shell=True )
        subprocess.call( "rm ./smo_chipseq_sim/ -rf", shell=True )
    """





