import sys
import numpy
import string
import random
from itertools import izip
from operator import itemgetter
import bisect 
import gzip    

import re
import array

import os
import shutil
import subprocess
import StringIO
import tempfile

import numpy

STATMAP_PATH = '../src/statmap'

#################################################################################################
### verbosity level information 
#
# whether or not to print statmap output
P_STATMAP_INPUT = False
if not P_STATMAP_INPUT:
    stdout = tempfile.TemporaryFile()
    stderr = tempfile.TemporaryFile()
else:
    stdout = sys.stdout
    stderr = sys.stderr

CLEANUP = True
    
### END verbosity level information  ############################################################


bps = [ 'a', 'A', 'c', 'C', 'g', 'G', 't', 'T', 'n', 'N' ]
bps_set = set( [ 'a', 'A', 'c', 'C', 'g', 'G', 't', 'T'] )
bp_complement = {'a': 't', 'A': 'T', 'c': 'g', 'C': 'G', 'g': 'c', 'G': 'C', 't': 'a', 'T': 'A'}

# build the translation table, for building the random genome
from_str = ''.join( str(i) for i in range( 10 ))
to_str = ''.join(bps)
translation_table = string.maketrans( from_str, to_str  )

# build a translation table for the rev complement
from_str = ''.join(bps[:-2])
to_str = ''.join( bp_complement[bp] for bp in from_str )
rev_comp_table = string.maketrans( from_str, to_str  )


#### Data Types ########

"Store mutated reads "
try:
    import collections
    mut_read_t = collections.namedtuple( "mutated_read", "mut_seq err_str seq")
except AttributeError:
    class mut_read_t(tuple):
        'mut_read_t(mut_seq err_str, seq)'
        
        __slots__ = ()
        
        _fields = ('mut_seq', 'err_str', 'seq')
        
        def __new__(_cls, x, y, z):
            return tuple.__new__(_cls, (x, y, z))
        
        @classmethod
        def _make(cls, iterable, new=tuple.__new__, len=len):
            'Make a new mut_read_t object from a sequence or iterable'
            result = new(cls, iterable)
            if len(result) != 3:
                raise TypeError('Expected 3 arguments, got %d' % len(result))
            return result
        
        def __repr__(self):
            return 'mut_read_t(mut_seq=%r, err_str=%r, seq=%r)' % self
        
        def _asdict(t):
            'Return a new dict which maps field names to their values'
            return {'mut_seq': t[0], 'err_str': t[1], 'seq':t[2]}
        
        def __getnewargs__(self):
            return tuple(self)
        
        mut_seq = property(itemgetter(0))
        err_str = property(itemgetter(1))        
        seq = property(itemgetter(2))


### SAM File parsing code

def iter_sam_reads( f ):
    # get the next line in the sam file
    next_line = f.readline().strip().split("\t")
    readname = next_line[0].split("/")[0]
    # read until we reach the end of the file
    while len(next_line) > 1:
        data = [ next_line, ]
        readname = next_line[0].split("/")[0]
        next_line = f.readline().strip().split("\t")
        # while the readnames are the same as the past readname
        while next_line[0].split("/")[0] == readname:
            # add the next line
            data.append( next_line )
            next_line = f.readline().strip().split("\t")
            # yield the split read lines
        yield data
    return

### SOLEXA Specific mutation routines ###############

## READ MUTATION Code ##############################

def calc_mutation_prob_from_solexa_char( quality_char ):
    numerical_value = ord(quality_char)- 64
    if numerical_value <= -5:
       numerical_value = -5 
    exp_value =10**(numerical_value/-10.0)
    return exp_value/(1+exp_value)

def calc_mutation_probs_from_solexa_string( quality_string ):
    return [ calc_mutation_prob_from_solexa_char(char) for char in quality_string ]

def mutate_solexa_read( read, quality_string ):
    rv = []
    for bp, char in izip( read, quality_string ):
        if random.random() > calc_mutation_prob_from_solexa_char( char ):
            rv.append( bp )
        else:
            rv.append( random.choice( bps[:-2] ) )
    return ''.join(rv).upper()

def get_error_strs_from_fastq( fastq_iterator ):
    error_lines = []
    for line_num, line in enumerate(fastq_iterator):
        # only look at error lines
        if line_num % 4 != 3: continue
        error_lines.append( 'h' + line.strip()[1:] )
    return error_lines

## END READ MUTATION Code ##########################

def calc_penalty( ref_seq, seq, error_str ):
    penalty = 0
    for bp, bp_ref, error_char in izip(seq, ref_seq, error_str):
        mut_prob = calc_mutation_prob_from_solexa_char( error_char )
        if bp.upper() == bp_ref.upper():
            penalty += math.log10( 1 - mut_prob )
        else:
            # accounting for uniform bp mutation with /3.0
            penalty += math.log10( mut_prob/3.0 ) 
    return penalty

def generate_solexa_error_strs_from_sample( population, num ):
    for loop in xrange( num ):
        sample = random.sample( population, 1 )
        yield sample[0]

### END SOLEXA Specific mutation routines ############


def build_random_genome( chr_sizes, chr_names, with_ns = False ):
    """Simulate a random uniform genome.
    
    Paramaters:
    chr_size: a list of chr_lengths, in bps
    chr_names: a list of chr names, corresponding with the chr sizes
    num_repeats: the number of times to repeat each chromosome
    with_ns: whether or not to put n's in a chromosome

    Returns:
    A dict of chr names, sequence strings
    """
    assert len( chr_sizes ) == len( chr_names )
        
    # set the upper index in the random numbers 
    # to determine whether or not we simulate N's
    if with_ns: ul = 9
    else: ul = 7
    
    rv = {}

    # iterate through the chrs
    for chr_name, chr_size in zip(chr_names, chr_sizes):
        r_bps = numpy.random.random_integers(0, ul, chr_size)
        r_bps = r_bps.astype('c').tostring()
        r_bps = r_bps.translate( translation_table )
        rv[chr_name] = r_bps
    
    return rv 

def write_genome_to_fasta( genome, fasta_fp, num_repeats=1 ):
    # the maximum length, in bp's, of each fasta line
    FA_LL = 50
    
    for chr_name in genome.keys():
        fasta_fp.write( ">%s\n" % chr_name)
        
        seq = genome[chr_name]*num_repeats
        chr_size = len(seq)

        for start, stop in izip( \
            xrange(0,chr_size,FA_LL), \
            xrange(FA_LL,chr_size+FA_LL+1,FA_LL) \
            ):
            fasta_fp.write( seq[start:stop] )
            fasta_fp.write( "\n" )
    
    return

def sample_uniformily_from_genome( genome, nsamples=100, frag_len=200 ):
    # store a list of the true fragments
    truth = []
    # calculate the genome lengths sum, to make
    # sure that our samples are uniform in the 
    # chr length
    chr_names = genome.keys()
    chr_names.sort()
    chr_lens = [ len(genome[key]) for key in chr_names ]
    # subtract frag len  to account for the missing sample range at the end
    chr_lens_cdr = [ sum(chr_lens[:i])-frag_len for i in xrange(len(genome)+1) ]
    chr_lens_cdr[0] = 0
    chr_lens_cdr = numpy.array( chr_lens_cdr )
    for loop in xrange( nsamples ):
        rnd_num = random.randint( 0, chr_lens_cdr[-1] - 1  )
        chr_index = chr_lens_cdr.searchsorted( rnd_num  ) - 1
        rnd_bp = random.randint( 0, chr_lens[chr_index] - frag_len - 1 )
        truth.append(  ( chr_names[chr_index], rnd_bp, rnd_bp+frag_len ) )
    
    return truth

def build_snps_from_genome( genome, num ):
    snps = []
    # calculate the genome lengths sum, to make
    # sure that our samples are uniform in the 
    # chr length
    chr_names = genome.keys()
    chr_names.sort()
    chr_lens = [ len(genome[key]) for key in chr_names ]
    chr_lens_cdr = [ sum(chr_lens[:i]) for i in xrange(len(genome)+1) ]
    chr_lens_cdr[0] = 0
    chr_lens_cdr = numpy.array( chr_lens_cdr )
    # sample all of the chr's
    samples = set()
    for loop in xrange( num ):
        rnd_num = random.randint( 0, chr_lens_cdr[-1] - 1  )
        chr_index = chr_lens_cdr.searchsorted( rnd_num  ) - 1
        # make sure that we dont get duplicates
        rnd_bp = random.randint( 0, chr_lens[chr_index] - 1 )
        while ( chr_index, rnd_bp ) in samples:
            rnd_bp = random.randint( 0, chr_lens[chr_index] - 1 )
        samples.add( ( chr_index, rnd_bp )  )
        ref_bp = genome[chr_names[chr_index]][rnd_bp].upper()
        bps = ['A', 'C', 'G', 'T']
        bps.remove( ref_bp.upper() )
        alt_bp = random.choice( bps )
        snps.append( ( chr_names[chr_index], rnd_bp, ref_bp, alt_bp )  )
        assert alt_bp.upper() != ref_bp.upper()
        assert genome[chr_names[chr_index]].upper() != alt_bp.upper()
    return snps

def write_snps_to_snpcov_file( snps, genome, fname  ):
    chr_names = genome.keys()
    chr_names.sort()
    
    fp = open( fname, "w" )
    for snp in snps:
        fp.write("%s	%i	ID	0	0	0	0	0	%s	%s	%s	H	H	50	50	Test	RSIM	BICKEL\n" \
                 % ( snp[0], snp[1], snp[2], snp[2], snp[3]  ) )
    fp.close()
    return 

def build_reads_from_fragments( \
    genome, fragments, read_len=35, rev_comp=True, paired_end=False ):
    reads = []
    for chr, start, stop in fragments:
        if paired_end:
            read_1 = genome[chr][start:(start+read_len)]
            read_2 = genome[chr][(stop-read_len):stop]
            if random.random() > 0.5:
                read_1 = read_1.translate( rev_comp_table )[::-1]
            else:
                read_2 = read_2.translate( rev_comp_table )[::-1]
            reads.append( ( read_1, read_2 ) )
        else:
            if random.random() > 0.5 or not rev_comp:
                read = genome[chr][start:(start+read_len)]
            elif rev_comp:
                read = genome[chr][(stop-read_len):(stop)]
            
            if rev_comp and random.random() > 0.5:
                read = read.translate( rev_comp_table )[::-1]
            reads.append( read )
    return reads

def build_single_end_fastq_from_mutated_reads( samples_iter, of=sys.stdout ):
    for sample_num, (sample, error_str, true_seq) in enumerate( samples_iter ):
        of.write("@%s\n" % sample_num )
        of.write(sample + "\n")
        of.write("+%s\n" % sample_num )
        of.write(error_str + "\n")

def build_single_end_fastq_from_seqs( samples_iter, of=sys.stdout ):
    for sample_num, seq in enumerate( samples_iter ):
        error_str = 'h'*len(seq)
        of.write("@%s\n" % sample_num )
        of.write(seq + "\n")
        of.write("+%s\n" % sample_num )
        of.write(error_str + "\n")

def build_paired_end_fastq_from_mutated_reads( mut_reads_iter, of1, of2 ):
    for sample_num, (sample, error_str, true_seq) in enumerate( mut_reads_iter ):
        sample_1, sample_2 = sample
        error_str_1, error_str_2 = error_str
        # write the first pair of the read
        of1.write("@%s/1\n" % sample_num )
        of1.write(sample_1 + "\n")
        of1.write("+%s/1\n" % sample_num )
        of1.write(error_str_1 + "\n")
        
        # write the second pair of the read
        of2.write("@%s/2\n" % sample_num )
        of2.write(sample_2 + "\n")
        of2.write("+%s/2\n" % sample_num )
        of2.write(error_str_2 + "\n")
    return

def build_paired_end_fastq_from_seqs( sample_iter, of1, of2 ):
    for sample_num, (sample_1, sample_2) in enumerate( sample_iter ):
        # write the first pair of the read
        of1.write("@%s/1\n" % sample_num )
        of1.write(sample_1 + "\n")
        of1.write("+%s/1\n" % sample_num )
        of1.write('h'*len(sample_1) + "\n")
        
        # write the second pair of the read
        of2.write("@%s/2\n" % sample_num )
        of2.write(sample_2 + "\n")
        of2.write("+%s/2\n" % sample_num )
        of2.write('h'*len(sample_2) + "\n")
    return

def mutate_reads( reads, error_strs ):
    # set paired end vs single enmd read parameters
    paired_end = (  isinstance( reads[0], tuple ) )
    n_error_strs = 2 if paired_end else 1

    read_len = len( reads[0][0] ) if paired_end else len(reads[0])
    error_str_len = len( error_strs[0] ) if paired_end else len(error_strs)
    multiplier = (1 + read_len/error_str_len)
    
    mutated_reads = []
    for read in reads:
        read_error_str = random.sample( error_strs, n_error_strs )
        if paired_end:
            err_str1 = (read_error_str[0]*multiplier)[:read_len]
            p1 = mutate_solexa_read( read[0], err_str1 )
            err_str2 = (read_error_str[1]*multiplier)[:read_len]
            p2 = mutate_solexa_read( read[1], err_str2 )[:read_len]
            mutated_reads.append( mut_read_t((p1, p2), (err_str1, err_str2), read ) )
        else:
            read_error_str = (read_error_str[0]*multiplier)[:read_len]
            mr = mutate_solexa_read( read, read_error_str )
            mutated_reads.append( mut_read_t(mr, read_error_str, read) )
    return mutated_reads

def build_expected_map_locations_from_repeated_genome( \
    truth, chr_len, num_repeats, paired_reads ):
    if paired_reads == True:
        for start, stop in truth:            
            read_starts = [ start + i*chr_len for i in xrange(num_repeats) ]
            read_stops = [ start + i*chr_len for i in xrange(num_repeats) ]
            
            truth_locs = set()
            for s in read_starts:
                for e in read_ends:
                    truth_locs.add( (s,e) )
                    truth_locs.add( (e,s) )
        return truth_locs
    
    if paired_reads == False:
        for start, stop in truth:            
            read_starts = [ start + i*chr_len for i in xrange(num_repeats) ]
            read_stops = [ start + i*chr_len for i in xrange(num_repeats) ]
        return set( read_starts ), set( read_stops )
    
###
# Test to make sure that we are correctly finding reverse complemented subsequences. These
# should all be short reads that we can map uniquely. We will test this over a variety of
# sequence lengths. 
def test_sequence_finding( read_len, rev_comp = False ):
    rl = read_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    fragments = sample_uniformily_from_genome( r_genome, nsamples=100, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=rev_comp, paired_end=False )
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_of = open("tmp.genome", "w")
    write_genome_to_fasta( r_genome, genome_of, 1 )
    genome_of.close()
    
    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_seqs( reads, reads_of )
    reads_of.close()

    call = "%s -g tmp.genome -r tmp.fastq \
                             -p %.2f -m %.2f \
                             -t 1 " % ( STATMAP_PATH, -10, 2 )
    
    print >> stdout, re.sub( "\s+", " ", call)
    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    # ret_code = ( os.system( call ) >> 8 )
    if ret_code != 0:
        raise ValueError, "TEST FAILED: statmap call returned error code '%s'" % str( ret_code )
        
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./statmap_output/mapped_reads.sam" )
    total_num_reads = sum( 1 for line in sam_fp )
    sam_fp.seek(0)

    if len(fragments) != total_num_reads:
        raise ValueError, "Mapping returned too few reads."
    
    
    for reads_data, truth in izip( iter_sam_reads( sam_fp ), fragments ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if len(reads_data) != 1:
            raise ValueError, "Mapping returned too many results."
        
        loc = ( reads_data[0][2], int(reads_data[0][3]) )
        
        # make sure the chr and start locations are identical
        if loc[0] != truth[0] \
           or loc[1] != truth[1]:
            raise ValueError, \
                "Truth (%s, %i) and Mapped Location (%s, %i, %i) are not equivalent" \
                % ( loc[0], loc[1], truth[0], truth[1], truth[2]  )
    
    sam_fp.close()
    
    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome")
        os.remove("./tmp.fastq")
    shutil.rmtree("./statmap_output/")


def test_fivep_sequence_finding( ):
    rls = [ 15, 25, 50, 75  ]
    for rl in rls:
        test_sequence_finding( rl, False )
        print "PASS: Forward Mapping %i BP Test. ( Statmap appears to be mapping 5', perfect reads correctly )" % rl

def test_threep_sequence_finding( ):
    rls = [ 15, 75  ]
    for rl in rls:
        test_sequence_finding( rl, True ) 
        print "PASS: Reverse Mapping %i BP test. ( Statmap appears to be mapping 3', perfect reads correctly )" % rl

def test_short_sequences():
    failed_lengths = []
    rls = [ 4, 5, 7, 8, 9, 11, 12, 13, 14  ]
    for rl in rls:
        try:
            test_sequence_finding( rl, False )
        except Exception, inst:
            failed_lengths.append( rl )
    if len( failed_lengths ) == 0:
        print "PASS: Short Read Mapping Passed All Index Lengths ( %s )" % ' '.join( map( str, rls ) ) 
    else:
        print "FAIL: Short Read Mapping Failed for Index Lengths: %s ( of %s )" \
            % ( ' '.join( map( str, failed_lengths ) ), ' '.join( map( str, rls ) ) ) 
        sys.exit( -1 )
###
# Test to make sure that we are correctly finding paired end reads. 
def test_paired_end_reads( read_len ):
    rl = read_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    assert 2*rl < 200 # make sure the fragments are long enough
    fragments = sample_uniformily_from_genome( r_genome, nsamples=100, frag_len=200 )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=True )

    # rev complement the second read in the pair
    #for i, read in enumerate(reads):
    #    print read[0], read[1].translate( rev_comp_table )
    #    reads[i] = ( read[0], read[1].translate( rev_comp_table ) )
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_of = open("tmp.genome", "w")
    write_genome_to_fasta( r_genome, genome_of, 1 )
    genome_of.close()
    
    # build and write the reads
    reads_of_1 = open("tmp.1.fastq", "w")
    reads_of_2 = open("tmp.2.fastq", "w")
    build_paired_end_fastq_from_seqs( reads, reads_of_1, reads_of_2 )
    reads_of_1.close()
    reads_of_2.close()

    call = "%s -g tmp.genome -1 tmp.1.fastq -2 tmp.2.fastq  \
                               -p %.2f -m %.2f \
                               -t 1 " % ( STATMAP_PATH, -10, 2 )
        
    print >> stdout, re.sub( "\s+", " ", call)
    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./statmap_output/mapped_reads.sam" )
    # we divide by two to account for the double write
    total_num_reads = sum( 1 for line in sam_fp )/2
    sam_fp.seek(0)

    if len(fragments) != total_num_reads:
        raise ValueError, "Mapping returned too few reads (%i/%i)." % ( total_num_reads, len(fragments) )
        
    for reads_data, truth in izip( iter_sam_reads( sam_fp ), fragments ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if len(reads_data) != 2:
            print reads_data
            raise ValueError, "Mapping returned too many results."
        
        loc = ( reads_data[0][2], int(reads_data[0][3]) )
        
        # make sure the chr and start locations are identical
        if loc[0] != truth[0] \
           or loc[1] != truth[1]:
            raise ValueError, \
                "Truth (%s, %i) and Mapped Location (%s, %i, %i) are not equivalent" \
                % ( loc[0], loc[1], truth[0], truth[1], truth[2]  )
    
    sam_fp.close()

    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome")
        os.remove("./tmp.1.fastq")
        os.remove("./tmp.2.fastq")
    shutil.rmtree("./statmap_output/")

def test_paired_end_sequence_finding( ):
    rls = [ 25, 75  ]
    for rl in rls:
        test_paired_end_reads( rl ) 
        print "PASS: Paired End Mapping %i BP Test. ( Statmap appears to be mapping randomly oriented, paired end perfect reads correctly )" % rl

### Test to make sure that duplicated reads are dealt with correctly ###
def test_duplicated_reads( read_len, n_chrs, n_dups, gen_len, n_threads ):
    rl = read_len
    GENOME_LEN = gen_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [GENOME_LEN]*n_chrs, map( str, range( 1, n_chrs+1 ) ) )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    fragments = sample_uniformily_from_genome( r_genome, nsamples=100, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_of = open("tmp.genome", "w")
    write_genome_to_fasta( r_genome, genome_of, n_dups )
    genome_of.close()
    
    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_seqs( reads, reads_of )
    reads_of.close()

    call = "%s -g tmp.genome -r tmp.fastq \
                             -p %.2f -m %.2f \
                             -t %i " % ( STATMAP_PATH, -10, 2, n_threads )
        
    print >> stdout, re.sub( "\s+", " ", call)
    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./statmap_output/mapped_reads.sam" )
    total_num_reads = sum( 1 for line in sam_fp )
    sam_fp.seek(0)

    if len(fragments)*n_dups != total_num_reads:
        raise ValueError, "Mapping returned too few reads."
    
    
    for reads_data, truth in izip( iter_sam_reads( sam_fp ), fragments ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if len(reads_data) != n_dups:
            raise ValueError, "Mapping returned incorrect number of results."
        
        
        locs = [ (reads_data[i][2], int(reads_data[i][3])) for i in xrange(len(reads_data)) ]
        
        # make sure the chr and start locations are identical
        for i, loc in enumerate( locs ):
            if loc[0] != truth[0] \
               or loc[1]%GENOME_LEN != truth[1]:
                raise ValueError, \
                    "Mapped Location (%s, %i) and Truth (%s, %i, %i) are not equivalent" \
                    % ( loc[0], loc[1]%GENOME_LEN, truth[0], truth[1], truth[2]  )
    sam_fp.close()

    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome")
        os.remove("./tmp.fastq")
    shutil.rmtree("./statmap_output/")

def test_repeat_sequence_finding( ):
    rls = [ 50, 75  ]
    for rl in rls:
        test_duplicated_reads( rl, n_chrs=3, n_dups=5, gen_len=1000, n_threads=1 ) 
        print "PASS: Multi-Chr and Repeated Chr Mapping %i BP Test. ( Statmap appears to be mapping multiple genome and chr with heavy perfect repeats correctly )" % rl

def test_lots_of_repeat_sequence_finding( ):
    rls = [ 16, ]
    for rl in rls:
        # setting n_threads to -1 makes it deafult to the number of avialable cores
        test_duplicated_reads( rl, n_chrs=1, n_dups=4000, gen_len=100, n_threads=-1 ) 
        print "PASS: lots of repeat seqs ( %i ) %i BP Test. ( This tests a genome with lots and lots of repeats ( pbly mostly corner cases )  )" % ( 4000, rl )


### Test to make sure that duplicated reads are dealt with correctly ###
def test_dirty_reads( read_len, min_penalty=-30, n_threads=1 ):
    rl = read_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. 
    fragments = sample_uniformily_from_genome( r_genome, nsamples=100, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )

    # mutate the reads by their error strings
    sample_file = gzip.open( './data/dirty_error_strs.fastq.gz' )
    error_strs = get_error_strs_from_fastq( sample_file )
    sample_file.close()
    mutated_reads = mutate_reads( reads, error_strs )
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_of = open("tmp.genome", "w")
    write_genome_to_fasta( r_genome, genome_of, 1 )
    genome_of.close()
    
    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_mutated_reads( mutated_reads, reads_of )
    reads_of.close()

    call = "%s -g tmp.genome -r tmp.fastq \
                             -p %.2f -m 20 \
                             -t %i " % ( STATMAP_PATH, min_penalty, n_threads )
        
    print >> stdout, re.sub( "\s+", " ", call)
    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./statmap_output/mapped_reads.sam" )
    mapped_read_ids = set( ( line.split("\t")[0].strip() for line in sam_fp ) )
    total_num_reads = len( mapped_read_ids ) 
    sam_fp.seek(0)

    # find the unmappable reads
    unmappable_fp = open( "./statmap_output/reads.unpaired.unmappable" )
    num_unmappable_reads = sum( 1 for line in unmappable_fp )/4
    unmappable_fp.seek(0)
    unmappable_reads_set = set( line.strip()[1:] for i, line in enumerate(unmappable_fp) if i%4 == 0  )
    unmappable_fp.close()

    all_read_ids = set(map(str, range(100)))
    all_read_ids = all_read_ids - mapped_read_ids
    all_read_ids = all_read_ids - unmappable_reads_set
    all_read_ids = list(all_read_ids)
    
    all_read_ids.sort( )

    if len(fragments) != total_num_reads + num_unmappable_reads:
        raise ValueError, "Mapping returned too few reads %i/( %i + %i ). NOT { %s }" \
            % ( len(fragments), total_num_reads, num_unmappable_reads, ','.join( all_read_ids  ) )

    # build a dictionary of mapped reads
    mapped_reads_dict = dict( (data[0][0], data) for data in iter_sam_reads(sam_fp) )
    
    unmapped_reads = set( map( str, xrange( 100 ) ) ).difference( unmappable_reads_set.union( set(mapped_reads_dict.keys()) ) )
    if len( unmapped_reads ) != 0:
          raise ValueError, " We are missing '%s' from the set of mappable reads. " % str( unmapped_reads )

    for mapped_reads in iter_sam_reads(sam_fp):
        while mutated_reads[index].mut_seq in unmappable_reads_strs:
            print "Unmappable %i" % cnted_unmappable_reads
            cnted_unmappable_reads += 1
            index += 1
        
        mutated_read = mutated_reads[index]
        truth = fragments[index]

        index += 1
        
        locs = [ (mapped_reads[i][2], int(mapped_reads[i][3])) for i in xrange(len(mapped_reads)) ]
        
        # make sure the chr and start locations are identical
        for i, loc in enumerate( locs ):
            if loc[0] != truth[0] \
               or loc[1] != truth[1]:
                print "INDEX: ", index, index-num_unmappable_reads
                print mutated_read
                print "Truth:", truth
                print r_genome[truth[0]][truth[1]:(truth[1]+rl)].upper()
                print mapped_reads
                print r_genome[loc[0]][loc[1]:(loc[1]+rl)].upper()
                raise ValueError, \
                    "Mapped Location (%s, %i) and Truth (%s, %i, %i) are not equivalent" \
                    % ( loc[0], loc[1], truth[0], truth[1], truth[2]  )
    
    sam_fp.close()

    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome")
        os.remove("./tmp.fastq")
    shutil.rmtree("./statmap_output/")

def test_mutated_read_finding( ):
    rls = [ 50, 75  ]
    for rl in rls:
        test_dirty_reads( rl, n_threads=1, min_penalty=-30 ) 
        print "PASS: Dirty Read (-30 penalty) Mapping %i BP Test. ( Statmap appears to be mapping fwd strand single reads with heavy errors correctly )" % rl
        # FIXME - do the work to fix these tests
        #test_dirty_reads( rl, n_threads=1, min_penalty=-1 ) 
        #print "PASS: Dirty Read (-1 penalty) Mapping %i BP Test. ( Statmap appears to be mapping fwd strand single reads with heavy errors correctly )" % rl

def test_multithreaded_mapping( ):
    rls = [ 50, 75  ]
    for rl in rls:
        try: 
            test_dirty_reads( rl ) 
        except:
            print "FAIL: Multi-Threaded Read Mapping %i BP Test Failed." % rl
            raise
        else:
            print "PASS: Multi-Threaded Read Mapping %i BP Test. ( Statmap appears to be mapping correctly with multiple threads )" % rl
###
# Test to make sure that we are correctly indexing and finding snps
def test_snps( read_len, num_snps = 10 ):
    rl = read_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # build the snps
    snps = build_snps_from_genome( r_genome, num_snps )
    snps.sort( key = itemgetter(1) )
    try:
        for snp in snps:
            assert r_genome[ snp[0] ][snp[1]].upper() != snp[3].upper()
            assert r_genome[ snp[0] ][snp[1]].upper() == snp[2].upper()
    except:
        print r_genome[ snp[0] ][snp[1]], snp
        raise
    
    write_snps_to_snpcov_file( snps, r_genome, "tmp.snpcov"  )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    fragments = sample_uniformily_from_genome( r_genome, nsamples=500, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )
 
    ref_snps = [0]*num_snps
    alt_snps = [0]*num_snps
    # mutate the snps
    for index, (read, fragment) in enumerate(zip(reads, fragments)):
        for snp_index, snp in enumerate(snps):
            # if this read covers the snp
            if snp[0] == fragment[0] \
               and snp[1] >= fragment[1] \
               and snp[1] < fragment[2]:
                if random.random() > 0.5:
                    tmp = array.array('c', reads[index])
                    tmp[snp[1]-fragment[1]] = snp[3]
                    assert reads[index] != tmp.tostring()
                    reads[index] = tmp.tostring()
                    alt_snps[snp_index] += 1
                else:
                    ref_snps[snp_index] += 1
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_of = open("tmp.genome", "w")
    write_genome_to_fasta( r_genome, genome_of, 1 )
    genome_of.close()
    
    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_seqs( reads, reads_of )
    reads_of.close()

    call = "%s -g tmp.genome -r tmp.fastq -s tmp.snpcov \
                             -p %.2f -m %.2f \
                             -t 1 " % ( STATMAP_PATH, -10.0, 0.50 )
        
    print >> stdout, re.sub( "\s+", " ", call)
    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    # ret_code = ( os.system( call ) >> 8 )
    if ret_code != 0:
        print "TEST FAILED - statmap call returned error code ", ret_code
        sys.exit( -1 )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./statmap_output/mapped_reads.sam" )
    total_num_reads = sum( 1 for line in sam_fp )
    sam_fp.seek(0)
    
    if len(fragments) != total_num_reads:
        raise ValueError, "Mapping returned too few reads. ( %i / %i )" % ( total_num_reads, len(fragments) )
            
    for loop, (reads_data, truth) in enumerate( izip( iter_sam_reads( sam_fp ), fragments ) ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if len(reads_data) != 1:
            raise ValueError, "Mapping returned too many results."

        if int(reads_data[0][0]) != loop:
            raise ValueError, "Key %i (%i) is off ( a read was probably skipped )" % ( loop, int(reads_data[0][0]) )
        
        loc = ( reads_data[0][2], int(reads_data[0][3]) )
        
        # make sure the chr and start locations are identical
        if loc[0] != truth[0] \
           or loc[1] != truth[1]:
            print reads_data
            raise ValueError, \
                "Truth (%s, %i) and Mapped Location (%s, %i, %i) are not equivalent" \
                % ( loc[0], loc[1], truth[0], truth[1], truth[2]  )
    
    sam_fp.close()

    ###### Test the snp file to make sure that each of the snps appears #############
    
    snp_fp = open( "./statmap_output/updated_snp_cnts.snp" )
    lines = list( snp_fp )
    snp_fp.close()
    for line_num, (line, ref_cnt, alt_cnt) in enumerate(zip( lines[1:], ref_snps, alt_snps)):
        data = map( int, line.strip().split("\t")[-2:] )
        if data[0] != ref_cnt or data[1] != alt_cnt:
            print line.strip(), data, [ref_cnt, alt_cnt]
            raise ValueError, "The wrong number of snps was discovered ( %s vs %s ) " % ( data, [ref_cnt, alt_cnt] )

    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome")
        os.remove("./tmp.fastq")
        os.remove("./tmp.snpcov")
    shutil.rmtree("./statmap_output/")
        

def test_snp_finding( ):
    rls = [ 25, ]
    for rl in rls:
        test_snps( rl )
        print "PASS: SNP Mapping %i BP Test. ( Statmap appears to be mapping perfect SNPs correctly )" % rl


if False:
    num_repeats = 1
    num_chrs = 1
    frag_len = 200
    paired = False
    chr_sizes = [450]*num_chrs
    chr_names = [ "chr%i" % i for i in xrange( num_chrs )  ]
    
    assert all( frag_len < chr_len for chr_len in chr_sizes )
    
    r_genome = build_random_genome( chr_sizes, chr_names )
    truth = sample_uniformily_from_genome( r_genome, nsamples=1000, frag_len=frag_len )
    reads = build_reads_from_fragments( r_genome, truth, paired_end=paired )
    
    sample_file = gzip.open( 'clean_error_strs.fastq.gz' )
    error_strs = get_error_strs_from_fastq( sample_file )
    sample_file.close()
    mutated_reads = mutate_reads( reads, error_strs )
    
    test_dirty_reads( mutated_reads, r_genome, truth, \
                      num_chr_repeats=1, \
                      min_penalty=-10, max_penalty_spread=2 )

if __name__ == '__main__':
    RUN_SLOW_TESTS = True

    print "Starting test_fivep_sequence_finding()"
    test_fivep_sequence_finding()
    print "Starting test_threep_sequence_finding()"
    test_threep_sequence_finding()
    print "Starting test_paired_end_sequence_finding()"
    test_paired_end_sequence_finding( )
    print "Starting test_repeat_sequence_finding()"
    test_repeat_sequence_finding()
    print "Starting test_mutated_read_finding()"
    test_mutated_read_finding()
    print "Starting test_multithreaded_mapping()"
    test_multithreaded_mapping( )
    print "Starting test_snp_finding()"
    test_snp_finding()
    
    # We skip this test because statmap can't currently
    # index reads less than 12 basepairs ( and it shouldn't: 
    #     we should be building a hash table for such reads )
    # test_short_sequences()
    
    if RUN_SLOW_TESTS:
        print "[SLOW] Starting test_lots_of_repeat_sequence_finding()"
        test_lots_of_repeat_sequence_finding( )
