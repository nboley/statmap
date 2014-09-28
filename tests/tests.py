#!/usr/bin/python

import sys
import numpy
import string
import random
from itertools import izip
from operator import itemgetter
import bisect 
import gzip    

from math import ceil

import re
import array

import os
import shutil
import subprocess
import StringIO
import tempfile

import numpy

# add python_lib to the path relative to the location of tests.py
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from error_model import load_error_data

STATMAP_PATH = '../bin/statmap'
BUILD_SAM_PATH = '../utilities/write_sam_from_mapped_reads_db.py'
BUILD_INDEX_PATH = '../utilities/build_index.py'
CALL_PEAKS_PATH = '../utilities/call_peaks.py'

################################################################################
### verbosity level information 
#
# whether or not to print statmap output
P_STATMAP_OUTPUT = True
if not P_STATMAP_OUTPUT:
    stdout = tempfile.TemporaryFile()
    stderr = tempfile.TemporaryFile()
else:
    stdout = sys.stdout
    stderr = sys.stderr

CLEANUP = False
    
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

def reverse_complement(seq):
    return seq.translate( rev_comp_table )[::-1]

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

def map_with_statmap( genome_fnames, read_fnames, output_dir, indexed_seq_len,
        mapping_metaparameter=None, num_threads=1, search_type="m", assay=None,
        num_samples=0, softclip_len=0 ):

    if not P_STATMAP_OUTPUT:
        stderr.seek( 0 )
        stderr.truncate()
        stdout.seek( 0 )
        stdout.truncate()
    
    # build the input fnames str
    assert len( read_fnames ) in (1,2)
    read_fname_str = None
    if 1 == len( read_fnames ):
        read_fname_str = "-r " + read_fnames[0]
    else:
        read_fname_str = "-1 " + read_fnames[0] + " -2 " + read_fnames[1]
    
    # build_index
    call = "%s %i tmp.genome.fa.bin %s" % ( BUILD_INDEX_PATH, indexed_seq_len,
            ' '.join(genome_fnames) )

    print >> stderr, "========", re.sub( "\s+", " ", call)
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        raise ValueError, "TEST FAILED: build_index call returned error code '%s'" \
            % str( ret_code )

    # run statmap
    call = "%s -g tmp.genome.fa.bin %s -o %s -t %i" \
        % ( STATMAP_PATH, read_fname_str, output_dir, num_threads )

    if mapping_metaparameter != None:
        call += " -p %.2f " % mapping_metaparameter

    if assay != None:
        call += " -a " + assay

    # Setting num_samples only makes sense if we also specify the assay, since
    # the iterative mapping code is assay-specific
    if num_samples > 0 and assay != None:
        call += " -n " + str(num_samples)

    # add the search type
    if search_type != None:
        call += " -s " + search_type

    if softclip_len > 0:
        call += " -S " + str(softclip_len)

    print >> stderr, "========", re.sub( "\s+", " ", call)    
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        raise ValueError, "TEST FAILED: statmap call returned error code '%s'" % str( ret_code )
    
    # build the sam file
    call = "%s %s > %s" % ( BUILD_SAM_PATH, output_dir, os.path.join(output_dir, "mapped_reads.sam") )
    print >> stderr, "========", re.sub( "\s+", " ", call)
    ret_code = subprocess.call( call, shell=True, stdout=stdout, stderr=stderr )
    if ret_code != 0:
        raise ValueError, "TEST FAILED: build_sam_from_mapped_reads call returned error code '%s'" \
            % str( ret_code )


### SAM File parsing code

def skip_sam_header( fp ):
    """Skips header in fp of SAM file, returning fp pointing to first read"""
    # skip header (all lines start with '@')
    while True:
        fpos = fp.tell()
        next_line = fp.readline()
        if next_line.startswith('@'):
            continue
        else:
            # restore file pointer to start of first read
            fp.seek(fpos)
            break

def convert_from_one_indexed(reads):
    """
    Subtract 1 from the location in all reads from SAM
    SAM is 1-indexed, so we convert to 0-indexed to make it cleaner to compare
    to generated sequence
    """
    for read in reads:
        read[3] = str( int(read[3]) - 1 )
    return reads

def iter_sam_reads( f ):

    # skip header (all lines start with '@')
    skip_sam_header( f )

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
        yield convert_from_one_indexed(data)
    return

def count_lines_in_sam( sam_fp ):
    """
    Counts number of lines, excluding the header (== number of mapped reads)
    in a SAM file
    """
    skip_sam_header( sam_fp )

    num_lines = sum( 1 for line in sam_fp )

    # reset sam_fp to beginning for next function
    sam_fp.seek(0)

    return num_lines

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

def write_genome_to_multiple_fastas( genome, file_prefix, num_repeats=1 ):
    # the maximum length, in bp's, of each fasta line
    FA_LL = 50
    files_out = []
    
    for chr_name in genome.keys():
        if file_prefix:
            fasta_fname = file_prefix + "_" + chr_name + ".fa"
        else:
            fasta_fname = chr_name + ".fa"

        fasta_fp = open( fasta_fname, "w" )
        
        fasta_fp.write( ">%s\n" % chr_name)
        
        seq = genome[chr_name]*num_repeats
        chr_size = len(seq)

        for start, stop in izip( \
                xrange(0,chr_size,FA_LL), \
                xrange(FA_LL,chr_size+FA_LL+1,FA_LL) \
            ):
            fasta_fp.write( seq[start:stop] )
            fasta_fp.write( "\n" )
        
        files_out.append( fasta_fname )
        fasta_fp.close()
    
    return files_out

def sample_uniformly_from_genome( genome, nsamples=100, frag_len=200 ):
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

def reverse_complement( seq ):
    return seq.translate( rev_comp_table )[::-1]

def build_reads_from_fragments(
    genome, fragments, read_len=35, rev_comp=True, paired_end=False ):
    reads = []
    for chr, start, stop in fragments:
        if paired_end:
            if random.random() > 0.5 or not rev_comp:
                read_1 = genome[chr][start:(start+read_len)]
                read_2 = genome[chr][(stop-read_len):stop]
                read_2 = reverse_complement( read_2 )
            else:
                read_1 = genome[chr][(stop-read_len):stop]
                read_1 = reverse_complement( read_1 )
                read_2 = genome[chr][start:(start+read_len)]
            reads.append( ( read_1, read_2 ) )
        else:
            if random.random() > 0.5 or not rev_comp:
                read = genome[chr][start:(start+read_len)]
            elif rev_comp:
                read = genome[chr][(stop-read_len):(stop)]
                read = reverse_complement( read )
            
            reads.append( read )
    return reads

def build_single_end_fastq_from_mutated_reads( samples_iter, of=sys.stdout ):
    for sample_num, (sample, error_str, true_seq) in enumerate( samples_iter ):
        of.write("@%s\n" % sample_num )
        of.write(sample + "\n")
        of.write("+%s\n" % sample_num )
        of.write(error_str + "\n")

#def build_single_end_fastq_from_seqs( samples_iter, of=sys.stdout,
#        untemplated_gs_perc=0.0, random_prefix_len=0 ):
#    for sample_num, seq in enumerate( samples_iter ):
#        num_untemplated_gs = random.randint(1,3) \
#            if random.random() < untemplated_gs_perc else 0
#        error_str = 'h'*len(seq)
#        of.write("@%s\n" % sample_num )
#        of.write('g'*num_untemplated_gs + seq[:(len(seq)-num_untemplated_gs)] + "\n")
#        of.write("+%s\n" % sample_num )
#        of.write(error_str + "\n")

def build_single_end_fastq_from_seqs( samples_iter, of=sys.stdout,
        untemplated_gs_perc=0.0, random_prefix_len=0 ):

    # (for now) only do one type of prepending transformation
    assert( untemplated_gs_perc == 0 or random_prefix_len == 0 )

    for sample_num, seq in enumerate( samples_iter ):

        if random_prefix_len > 0:
            prefix = ''.join( [ random.choice(bps)
                    for n in range(random_prefix_len) ] )
        else:
            num_untemplated_gs = random.randint(1,3) \
                    if random.random() < untemplated_gs_perc else 0
            prefix = 'g'*num_untemplated_gs
            # trim off start of seq so all reads are the same length
            seq = seq[num_untemplated_gs:]

        seq = prefix + seq

        error_str = 'h'*len(seq)

        of.write("@%s\n" % sample_num )
        of.write(seq + "\n")
        of.write("+%s\n" % sample_num )
        of.write(error_str + "\n")

def build_paired_end_fastq_from_mutated_reads( 
        mut_reads_iter, of1, of2, num_untemplated_gs=0 ):
    for sample_num, (sample, error_str, true_seq) in enumerate( mut_reads_iter ):
        sample_1, sample_2 = sample
        error_str_1, error_str_2 = error_str
        # write the first pair of the read
        of1.write("@%s/1\n" % sample_num )
        of1.write('g'*num_untemplated_gs 
                  + sample_1[:(len(sample_1)-num_untemplated_gs)] + "\n")
        of1.write("+%s/1\n" % sample_num )
        of1.write(error_str_1 + "\n")
        
        # write the second pair of the read
        of2.write("@%s/2\n" % sample_num )
        of2.write('g'*num_untemplated_gs 
                  + sample_2[:(len(sample_2)-num_untemplated_gs)] + "\n")
        of2.write("+%s/2\n" % sample_num )
        of2.write(error_str_2 + "\n")
    return

def build_paired_end_fastq_from_seqs( sample_iter, of1, of2, num_untemplated_gs=0 ):
    for sample_num, (sample_1, sample_2) in enumerate( sample_iter ):
        # write the first pair of the read
        of1.write("@%s/1\n" % sample_num )
        of1.write('g'*num_untemplated_gs + sample_1[:(len(sample_1)-num_untemplated_gs)] + "\n")
        of1.write("+%s/1\n" % sample_num )
        of1.write('h'*len(sample_1) + "\n")
        
        # write the second pair of the read
        of2.write("@%s/2\n" % sample_num )
        of2.write('g'*num_untemplated_gs + sample_2[:(len(sample_2)-num_untemplated_gs)] + "\n")
        of2.write("+%s/2\n" % sample_num )
        of2.write('h'*len(sample_2) + "\n")
    return

def mutate_reads( reads, error_strs ):
    # set paired end vs single enmd read parameters
    paired_end = (  isinstance( reads[0], tuple ) )
    n_error_strs = 2 if paired_end else 1

    mutated_reads = []
    for read in reads:
        read_error_str = random.sample( error_strs, n_error_strs )

        read_len = len( read[0] ) if paired_end else len(read)
        error_str_len = len( read_error_str[0] ) if paired_end else len(read_error_str)
        multiplier = (1 + float(read_len)/error_str_len )
        if int(multiplier) < multiplier:
            multiplier = int(multiplier) + 1
        multiplier = int(multiplier)
        
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
    
        assert 1 == len( set( map( len, mutated_reads[-1]  ) ) ) 
    
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

def check_sequence_match( mapped_read, truth, genome ):
    # Get the start and length of the fragment in the genome that this read
    # maps to (including soft clipping, introns, etc.)
    cigar_re = r'([0-9]+)([MNS])'
    cigar_string = mapped_read[5]
    cigar_entries = re.findall( cigar_re, cigar_string )

    mapping_seq = mapped_read[2]
    mapping_start = int(mapped_read[1])

    genome_pos = mapping_start
    mapping_pos = 0

    for entry in cigar_entries:
        entry_len = int( entry[0] )
        entry_op = entry[1]

        if entry_op == 'M':
            genome_segment = genome[truth[0]][genome_pos:genome_pos+entry_len]
            mapping_segment = mapping_seq[mapping_pos:mapping_pos+entry_len]

            if genome_segment.upper() != mapping_segment.upper():
                print "Original:", genome_segment
                print "Mapped  :", mapping_segment
                raise ValueError, "Mapped sequence in read_id %i does not match the genome" \
                        % (mapped_read[4])

            genome_pos += entry_len
            mapping_pos += entry_len

        elif entry_op == 'N':
            genome_pos += entry_len # skip the intron (in the genome sequence)
        elif entry_op == 'S':
            mapping_pos += entry_len # skip the soft clipped bases (in the read)
        else:
            print "ERROR : Invalid CIGAR entry op '%c'" % entry_op
            sys.exit(1)

    return
    
###
# Test to make sure that we are correctly finding reverse complemented
# subsequences. These should all be short reads that we can map uniquely. We
# will test this over a variety of sequence lengths. 
def test_sequence_finding( read_len, rev_comp = False, indexed_seq_len=None,
        untemplated_gs_perc=0.0, search_type="m", mapping_metaparameter=None,
        num_samples=0, assay=None, random_prefix_len=0 ):

    output_directory = "smo_test_sequence_finding_%i_rev_comp_%s_%s_%s" % ( \
        read_len, str(rev_comp), indexed_seq_len, search_type )

    # If no index sequence length is provided, use the length of the read
    indexed_seq_len = indexed_seq_len or read_len
    
    # if we are testing untemplated g's, increase the number of samples so that
    # we get some beginning overlap
    nsamples = 10000 if untemplated_gs_perc > 0 else 100
    
    ###### Prepare the data for the test #######################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # sample uniformly from the genome. This gives us the sequences that we
    # need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    fragments = sample_uniformly_from_genome( r_genome, nsamples=nsamples,
            frag_len=read_len )
    reads = build_reads_from_fragments( r_genome, fragments, read_len=read_len,
            rev_comp=rev_comp, paired_end=False )

    ###### Write out the test files, and run statmap ###########################
    # write genome
    genome_fnames = ( "tmp.genome.fa", )
    with open( genome_fnames[0], "w" ) as genome_of:
        write_genome_to_fasta( r_genome, genome_of, 1 )
    
    # build and write the reads
    read_fnames = ( "tmp.fastq", )
    with open(read_fnames[0], "w") as reads_of:
        build_single_end_fastq_from_seqs( reads, reads_of, untemplated_gs_perc,
                random_prefix_len )

    # map the data
    if assay == None:
        assay = None if untemplated_gs_perc == 0.0 else 'a'

    softclip_len = random_prefix_len

    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len, assay=assay, search_type=search_type,
            num_samples=num_samples,
            mapping_metaparameter=mapping_metaparameter,
            softclip_len=softclip_len )
    
    ###### Test the sam file to make sure that each of the reads appears #######
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    total_num_reads = count_lines_in_sam( sam_fp )

    if len(fragments) > total_num_reads:
        missing_reads = set( xrange(nsamples) ).difference(
                set( [ int(read[0][0]) for read in iter_sam_reads(sam_fp) ] ))
        print "Mapping returned the wrong number of reads ( %i vs expected %i )." % ( total_num_reads, len(fragments) )
        print "Missing reads:", ','.join( [ str(read_id) for read_id in missing_reads ] )
        sys.exit(1)

    sam_fp.seek(0)
    
    for reads_data, truth in izip( iter_sam_reads( sam_fp ), fragments ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if untemplated_gs_perc == 0.0 and len(reads_data) != 1:
            raise ValueError, "Mapping for readid %s returned too many results (read has %i mappings, expected 1)." \
                    % ( reads_data[0][0], len(reads_data) )
        
        locs = [ ( read_data[2], int(read_data[3]), read_data[9],
                   int(read_data[1]), int(read_data[0]), read_data[5] )
                   for read_data in reads_data ]
        
        # loc := (chr, start, sequence, flag, read_id, cigar )
        # truth := (chr, start, stop)

        for i, loc in enumerate(locs):
            # make sure the chr and start locations match
            if loc[0] != truth[0]:
                raise ValueError, "Mapping for read_id %i mapped to the wrong chromosome - expected %s, got %s" \
                        % (loc[4], loc[0], truth[0])

            # TODO proper tests for untemplated G reads

            if loc[1] != truth[1]:
                raise ValueError, "Mapping for read_id %i mapped to the wrong location - expected %i, got %i" \
                        % (loc[4], truth[1], loc[1])

            check_sequence_match( loc, truth, r_genome )

    sam_fp.close()
    
    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.fastq")
        shutil.rmtree(output_directory)


def test_fivep_sequence_finding( ):
    rls = [ 25, 50, 75  ]
    for rl in rls:
        test_sequence_finding( rl, False )
        print "PASS: Forward Mapping %i BP Test. ( Statmap appears to be mapping 5', perfect reads correctly )" % rl

def test_untemplated_g_finding( ):
    #rls = [ 25, 50, 75  ]
    rls = [ 25, ]
    for rl in rls:
        test_sequence_finding( rl, False, rl-4, untemplated_gs_perc=0.25,
                search_type='e' )
        print "PASS: Untemplated Gs %i BP Test. ( Statmap appears to be mapping 5', perfect reads correctly )" % rl

def test_threep_sequence_finding( ):
    rls = [ 75  ]
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

def test_build_index():
    rls = [ 25, 50, 75  ]
    for rl in rls:
        test_sequence_finding( rl, False, rl )
        print "PASS: Build Index w/ Forward Mapping %i BP Test. ( Statmap appears to be externally building/loading index correctly )" % rl

def test_iterative_mapping():
    """Run Statmap on some basic data and do iterative mapping just to see if it works"""
    # TODO why is this so slow?
    rls = [ 25, ]
    for rl in rls:
        test_dirty_reads( rl, assay='a', num_samples=1 ) # CAGE
        print "PASS: Iterative mapping code (%i BP reads) is at least running..." \
                % rl

def test_short_index_probe():
    rls = [ 25, 50, 75  ]
    for rl in rls:
        test_sequence_finding( rl, False, 20 )
        print "PASS: Short Index Probe (20 bp) w/ Forward Mapping %i BP Test. ( Statmap appears to be correctly mapping reads with an index with a shorter index probe )" % rl
        test_sequence_finding( rl, True, 20 )
        print "PASS: Short Index Probe (20 bp) w/ Backward Mapping %i BP Test. ( Statmap appears to be correctly mapping reads with an index with a shorter index probe )" % rl


###
# Test to make sure that we are correctly finding paired end reads. 
def test_paired_end_reads( read_len ):
    output_directory = "smo_test_paired_end_reads_%i" % ( read_len )    

    rl = read_len
    nsamples=100

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [1000,], ["1",] )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    assert 2*rl < 200 # make sure the fragments are long enough
    fragments = sample_uniformly_from_genome( r_genome, nsamples=nsamples, frag_len=200 )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=True )
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_fnames = ("tmp.genome.fa",)
    with open("tmp.genome.fa", "w") as genome_of:
        write_genome_to_fasta( r_genome, genome_of, 1 )
    
    # build and write the reads
    read_fnames = ( "tmp.1.fastq", "tmp.2.fastq" )
    reads_of1 = open( read_fnames[0], "w" )
    reads_of2 = open( read_fnames[1], "w" )
    build_paired_end_fastq_from_seqs( reads, reads_of1, reads_of2 )
    reads_of1.close()
    reads_of2.close()

    # map the reads - indexed_seq_len defaults to read_len
    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len=read_len )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    # we divide by two to account for the double write
    total_num_reads = count_lines_in_sam(sam_fp) / 2

    if len(fragments) != total_num_reads:
        missing_reads = set( xrange(nsamples) ).difference(
                set( [ int(read[0][0]) for read in iter_sam_reads(sam_fp) ] ))
        print "Mapping returned the wrong number of reads ( %i vs expected %i )." % ( total_num_reads, len(fragments) )
        print "Missing reads:", ','.join( [ str(read_id) for read_id in missing_reads ] )
        sys.exit(1)
        
    for reads_data, truth in izip( iter_sam_reads( sam_fp ), fragments ):
        # FIXME BUG - make sure that there arent false errors ( possible, but unlikely )
        if len(reads_data) != 2:
            print >> stderr, reads_data
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
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.1.fastq")
        os.remove("./tmp.2.fastq")
        shutil.rmtree(output_directory)

def test_paired_end_sequence_finding( ):
    rls = [ 25, 75  ]
    for rl in rls:
        test_paired_end_reads( rl ) 
        print "PASS: Paired End Mapping %i BP Test. ( Statmap appears to be mapping randomly oriented, paired end perfect reads correctly )" % rl

### Test to make sure that duplicated reads are dealt with correctly ###
def test_duplicated_reads( read_len, n_chrs, n_dups, gen_len, n_threads,
        n_reads=100, indexed_seq_len=None ):
    output_directory = "smo_test_duplicated_reads_%i_%i_%i_%i_%i_%i" % ( read_len, n_chrs, n_dups, gen_len, n_threads, n_reads )
    
    # If no indexed_seq_len explicitly set, use read_len
    indexed_seq_len = indexed_seq_len or read_len

    rl = read_len
    GENOME_LEN = gen_len

    ###### Prepare the data for the test ############################################
    # build a random genome
    r_genome = build_random_genome( [GENOME_LEN]*n_chrs, map( str, range( 1, n_chrs+1 ) ) )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in the
    # correct direction. ( ie 5 prime )
    fragments = sample_uniformly_from_genome( r_genome, nsamples=n_reads, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )
    # note: if we do rev_comp, statmap still correctly maps the read, but our comparison back to
    # the original genome will fail (incorrectly). Since comparison back to the genome is an 
    # important test for diploid mapping (to make sure contigs actually make sense), we don't rev_comp.
    
    ###### Write out the test files, and run statmap ################################
    # write genome
    genome_fnames = ( "tmp.genome.fa", )
    with open( genome_fnames[0] , "w" ) as genome_of:
        write_genome_to_fasta( r_genome, genome_of, n_dups )
    
    # build and write the reads
    read_fnames = ( "tmp.fastq", )
    with open(read_fnames[0], "w") as reads_of:
        build_single_end_fastq_from_seqs( reads, reads_of )

    # map with mismatches, since we won't be
    # able to bootstrap error scores from a
    # perfectly repeated genome
    # (no unique mappers)
    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len, num_threads=n_threads, search_type='m' )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./%s/mapped_reads.sam"  % output_directory )
    total_num_reads = count_lines_in_sam(sam_fp)

    if len(fragments)*n_dups != total_num_reads:
        raise ValueError, "Mapping returned too few reads (%s/%s)" % \
                (total_num_reads, len(fragments)*n_dups)
    
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
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.fastq")
        shutil.rmtree(output_directory)

def test_repeat_sequence_finding( ):
    rls = [ 50, 75  ]
    for rl in rls:
        test_duplicated_reads( rl, n_chrs=3, n_dups=5, gen_len=1000, n_threads=1 ) 
        print "PASS: Multi-Chr and Repeated Chr Mapping %i BP Test. ( Statmap appears to be mapping multiple genome and chr with heavy perfect repeats correctly )" % rl

def test_lots_of_repeat_sequence_finding( ):
    rls = [ 25, ]
    for rl in rls:
        # setting n_threads to -1 makes it deafult to the number of avialable cores
        test_duplicated_reads( rl, n_chrs=1, n_dups=4000, gen_len=100, n_threads=-1, n_reads=100 ) 
        print "PASS: lots of repeat seqs ( %i ) %i BP Test. ( This tests a genome with lots and lots of repeats ( pbly mostly corner cases )  )" % ( 4000, rl )


### Test to make sure that duplicated reads are dealt with correctly ###
def test_dirty_reads( read_len, n_threads=1, nreads=100, separate_fastas=False,
        assay=None, num_samples=0 ):
    output_directory = "smo_test_dirty_reads_%i_%i_%i" \
            % ( read_len, n_threads, nreads)

    if separate_fastas:
        output_directory += "_separate_fastas"
    
    rl = read_len

    ###### Prepare the data for the test ######################################
    # build a random genome
    r_genome = build_random_genome( [1000, 1000, 10000], ["1","2","3"] )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. 
    fragments = sample_uniformly_from_genome( r_genome, nsamples=nreads, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )
    
    # mutate the reads by their error strings
    sample_file = gzip.open( './data/dirty_error_strs.fastq.gz' )
    error_strs = get_error_strs_from_fastq( sample_file )
    sample_file.close()
    mutated_reads = mutate_reads( reads, error_strs )
    
    ###### Write out the test files, and run statmap ##########################
    # write genome
    if separate_fastas:
        genome_fnames = write_genome_to_multiple_fastas(r_genome, "tmp.genome")
    else:
        genome_fnames = ( "tmp.genome.fa", )
        with open( genome_fnames[0], "w" ) as genome_of:
            write_genome_to_fasta( r_genome, genome_of )

    # build and write the reads
    read_fnames = ( "tmp.fastq", )
    with open(read_fnames[0], "w") as reads_of:
        build_single_end_fastq_from_mutated_reads( mutated_reads, reads_of )

    map_with_statmap( genome_fnames, read_fnames, output_directory,
                      indexed_seq_len = read_len - 2,
                      num_threads=n_threads,
                      assay=assay, num_samples=num_samples )
    
    ###### Test the sam file to make sure that each of the reads appears ############
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    mapped_read_ids = set( ( line.split("\t")[0].strip() 
                             for line in sam_fp if not line.startswith('@') ) )
    total_num_reads = len( mapped_read_ids ) 
    sam_fp.seek(0)

    # find the unmappable reads
    unmappable_fp = open( "./%s/reads.unpaired.unmappable" % output_directory )
    num_unmappable_reads = sum( 1 for line in unmappable_fp )/4
    unmappable_fp.seek(0)
    unmappable_reads_set = set( line.strip()[1:] for i, line in enumerate(unmappable_fp) if i%4 == 0 )
    unmappable_fp.close()

    # find the nonmapping reads
    nonmapping_fp = open( "./%s/reads.unpaired.nonmapping" % output_directory )
    num_nonmapping_reads = sum(1 for line in nonmapping_fp)/4
    nonmapping_fp.seek(0)
    nonmapping_reads_set = set( line.strip()[1:] for i, line in enumerate(nonmapping_fp) if i%4 == 0 )
    nonmapping_fp.close()

    all_read_ids = set(map(str, range(nreads)))
    all_read_ids = all_read_ids - mapped_read_ids
    all_read_ids = all_read_ids - unmappable_reads_set
    all_read_ids = list(all_read_ids)
    
    all_read_ids.sort( )

    if len(fragments) != total_num_reads + num_unmappable_reads + num_nonmapping_reads:
        print "Mapping returned too few reads %i/( %i + %i ). NOT { %s }" \
            % ( len(fragments), total_num_reads, num_unmappable_reads, ','.join( all_read_ids  ) )
        sys.exit(1)
    
    # build a dictionary of mapped reads
    mapped_reads_dict = dict( (data[0][0], data) for data in iter_sam_reads(sam_fp) )
    
    unmapped_reads = set( map( str, xrange( nreads ) ) ).difference(
            nonmapping_reads_set.union(
                unmappable_reads_set.union(
                    set(mapped_reads_dict.keys())
        ) ))

    if len( unmapped_reads ) != 0:
        for key, entry in enumerate(mutated_reads):
            print key, entry
        for key, read in [ ( key, mutated_reads[ int(key) ] ) for key in unmapped_reads ]:
            print key, read
            print key, fragments[int(key)]
            print key, r_genome[fragments[int(key)][0]][fragments[int(key)][1]:fragments[int(key)][2]]
            
            genome_seq = r_genome[fragments[int(key)][0]][fragments[int(key)][1]:fragments[int(key)][2]]
            mut_seq = mutated_reads[int(key)][0]
            for bp1, bp2 in zip( genome_seq.upper(), mut_seq  ):
                print bp1, bp2
            
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
        for genome_fname in genome_fnames:
            os.remove(genome_fname)
        for read_fname in read_fnames:
            os.remove(read_fnames)
        shutil.rmtree(output_directory)

def test_mutated_read_finding( ):
    rls = [ 50, 75  ]
    for rl in rls:
        test_dirty_reads( rl ) 
        print "PASS: Dirty Read (-30 penalty) Mapping %i BP Test. ( Statmap appears to be mapping fwd strand single reads with heavy errors correctly )" % rl

def test_multithreaded_mapping( ):
    rls = [ 50, 75  ]
    for rl in rls:
        try: 
            test_dirty_reads( rl, n_threads=2 ) 
        except:
            print "FAIL: Multi-Threaded Read Mapping %i BP Test Failed." % rl
            raise
        else:
            print "PASS: Multi-Threaded Read Mapping %i BP Test. ( Statmap appears to be mapping correctly with multiple threads )" % rl

def test_multi_fasta_mapping( ):
    rl = 50
    try: 
        test_dirty_reads( rl, n_threads=2, separate_fastas=True ) 
    except:
        print "FAIL: Multi-Fasta Read Mapping %i BP Test Failed." % rl
        raise
    else:
        print "PASS: Multi-Fasta Read Mapping %i BP Test. ( Statmap appears to be mapping correctly from a genome with multiple fasta files )" % rl

def build_diploid_genome( seq_len, chr_name="chr1", gen_len=1000, n_dups=1,
        n_mut=10 ):
    '''
    Generates paternal and maternal genomes with a .map file for testing diploid mapping
    '''
    paternal_chr_name = chr_name + "_paternal"
    maternal_chr_name = chr_name + "_maternal"

    # Initialize diploid map file
    mapf = []
    mapf.append("#REF\tPAT\tMAT") # MAP header
    mapf.append('{0}\t{1}\t{2}'.format(0, 1, 1)) # sequence start

    # build a random paternal genome
    initial_genome = build_random_genome( [gen_len,], [chr_name,], )

    # create array copies of sequence to mutate
    mutated_paternal_chr = array.array( 'c', initial_genome[chr_name] )
    mutated_maternal_chr = array.array( 'c', initial_genome[chr_name] )

    # Mutate sequence, building map file along the way

    # indices to mutate; sorted so we can add entries to the map file
    # xrange(1, n) so we don't mutate the first bp
    muts = sorted(
        random.sample( xrange(1, len(initial_genome[chr_name])), n_mut )
    )

    # mut is the bp originally chosen to be the site of a mutation
    # keep track of the offsets from mut caused by indels
    p_offset = 0
    m_offset = 0
    for mut in muts:
        mut_type = random.choice( ['snp', 'indel'] )

        if mut_type == 'snp':
            # mutate (maternal - doesn't matter) sequence
            curr_bp = mutated_maternal_chr[ mut+m_offset ]
            valid_bps = [ bp for bp in bps_set if bp != curr_bp ]
            mutated_maternal_chr[ mut+m_offset] = random.choice( valid_bps )

        elif mut_type == 'indel':
            # insert random num bp's on paternal or maternal
            insertion_site = random.choice( ['paternal', 'maternal'] )
            insertion_len = random.choice( xrange(1, 6) )

            # insert random bps, add entries to map file, and update offsets
            # +1 because sequence array is 0-indexed, but .map files are 1-indexed
            random_bps = [ random.choice( bps[:8] ) for x in xrange( insertion_len ) ]
            if insertion_site == 'paternal':
                for random_bp in random_bps:
                    mutated_paternal_chr.insert( mut+p_offset, random_bp )
                mapf.append('{0}\t{1}\t{2}'.format(0, mut+p_offset+1, 0))
                p_offset += insertion_len
                mapf.append('{0}\t{1}\t{2}'.format(0, mut+p_offset+1, mut+m_offset+1))
            elif insertion_site == 'maternal':
                for random_bp in random_bps:
                    mutated_maternal_chr.insert( mut+m_offset, random_bp )
                mapf.append('{0}\t{1}\t{2}'.format(0, 0, mut+m_offset+1))
                m_offset += insertion_len
                mapf.append('{0}\t{1}\t{2}'.format(0, mut+p_offset+1, mut+m_offset+1))

    # build diploid genome (dictionary of sequences keyed by chr names)
    genome = {
        paternal_chr_name : mutated_paternal_chr.tostring(),
        maternal_chr_name : mutated_maternal_chr.tostring()
    }

    # write genome to multiple fasta files, with filename prefixes matching chr_names
    # duplicate genome n_dups times to simulate heavily repeated sequence
    output_filenames = write_genome_to_multiple_fastas( genome, "", n_dups )
    # write .map file
    map_fname = chr_name + "_test.map" # filename has to be chrname_x to be parsed by build_index
    with open( map_fname, "w" ) as f:
        f.write( '\n'.join(mapf) )
    output_filenames.append(map_fname)

    return genome, output_filenames

def map_diploid_genome( genome, genome_fnames, read_len, nreads=1000 ):
    '''
    Given a diploid genome, randomly sample reads and map them with statmap
    '''
    # sample reads uniformly from both genomes to get reads
    fragments = sample_uniformly_from_genome( genome, nsamples=nreads,
            frag_len=read_len )
    # note: if we do rev_comp, statmap still correctly maps the read, but our
    # comparison back to the original genome will fail (incorrectly). Since
    # comparison back to the genome is an important test for diploid mapping
    # (to make sure contigs actually make sense), we don't rev_comp.
    reads = build_reads_from_fragments(
            genome, fragments, read_len=read_len,
            rev_comp=False, paired_end=False
        )

    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_seqs( reads, reads_of )
    reads_of.close()

    # map the data
    output_directory = "smo_test_diploid_mapping_%i_%i_%i" % (
        read_len, nreads, len(genome) )
    read_fnames = [ "tmp.fastq" ]
    map_with_statmap( genome_fnames, read_fnames,
            output_directory, indexed_seq_len=read_len, search_type='e' )

    # test the sam file to make sure that each of the reads appears
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )

    # make sure the chr and start locations are identical
    for mapped_reads, truth in izip( iter_sam_reads(sam_fp), fragments ):
        locs = [ (
                    mapped_reads[i][2],         # chr
                    int(mapped_reads[i][3]),    # start pos
                    mapped_reads[i][9],         # sequence
                 )
                for i in xrange(len(mapped_reads)) ]

        # mapping reads to one diploid chr should return a maximum of 2 results
        # for each read
        if len(mapped_reads) > 2:
            print "Truth: ", truth
            print \
"A read sampled from a single diploid chromosome returned more than 2 mappings."
            print \
"You are probably indexing some portion of the genome more than once."
            raise ValueError, \
"Mapped_reads for read_id %i contains %i mappings; should have a maximum of 2" \
                % ( int(mapped_reads[0][0]), len(mapped_reads) )

        # compare mapped reads with original genome to make sure they match
        found_read = False
        for loc in locs:
            # we may randomly choose a read from one chr that is in fact shared on both
            # this is not an error - but we do want to make sure we got the read we wanted
            # so, check that at least one of the reads matches Truth
            # if chr_name and start_bp match
            if truth[0] == loc[0] and truth[1] == loc[1]: 
                # check original sequence
                if genome[loc[0]][truth[1]:truth[2]] == loc[2]:
                    found_read = True

#            since we add SNPs, it is valid to map a read sampled from one
#            chr to the opposite chr and have some mismatches (caused by snps
#            or insertions on the opposite chr). Therefore this check is
#            invalid.
#
#            # compare mapped read back to original genome
#            if genome[loc[0]][loc[1]:loc[1]+read_len].upper() != loc[2].upper():
#                print "Truth  : ", truth
#                print "Loc    : ", loc
#                print "Genome : %s" % genome[loc[0]][loc[1]:loc[1]+read_len]
#                print "Read   : %s" % loc[2]
#                raise ValueError, \
#                    "Error: Readid %i sequence failed to match at Genome (%s, %i) and Read (%s, %i)" \
#                    % ( int(mapped_reads[0][0]), truth[0], truth[1], loc[0], loc[1] )

        if found_read == False:
            print "Truth: ", truth
            print mapped_reads
            print "Original read not found in mapped locations - you probably failed to map the prior read."
            raise ValueError, \
                "Mapped locations at read id %i and Truth (%s, %i, %i) are not equivalent" \
                % ( int(mapped_reads[0][0]), truth[0], truth[1], truth[2] )

    sam_fp.close()

    # Cleanup the created files
    if CLEANUP:
        for fn in (read_fnames + output_filenames):
            os.remove(fn)
        shutil.rmtree(output_directory)

def map_duplicated_diploid_genome( genome, genome_fnames, read_len, genome_len,
        n_dups, nreads=100, indexed_seq_len=None, num_threads=1 ):
    '''
    Given a diploid genome, randomly sample reads and map with statmap.
    Test output of SAM from basis with expectation of duplicates
    '''
    output_directory = "smo_test_diploid_mapping_%i" % (read_len)

    indexed_seq_len = indexed_seq_len or read_len

    # sample reads uniformly from both genomes to get reads
    fragments = sample_uniformly_from_genome( genome, nsamples=nreads,
            frag_len=read_len )
    reads = build_reads_from_fragments( genome, fragments, read_len=read_len,
            rev_comp=False, paired_end=False)

    # build and write the reads
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_seqs( reads, reads_of )
    reads_of.close()

    # map the data
    read_fnames = ( "tmp.fastq", )
    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len, num_threads=num_threads )

    # test the sam file to make sure that each of the reads appears
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    total_num_reads = count_lines_in_sam( sam_fp )

    # Since the paternal and maternal chrs of the reference genome are
    # identical, each read will be repeated n_dups times * 2, since there will
    # be a paternal and maternal copy for every sampled read
    if len(fragments)*n_dups*2 != total_num_reads:
        raise ValueError, "Mapping returned the wrong number of reads (expected %i, got %i)" \
                % ( len(fragments)*n_dups*2, total_num_reads )

    for mapped_reads, truth in izip( iter_sam_reads(sam_fp), fragments ):

        if len(mapped_reads) != n_dups*2:
            raise ValueError, "Mapping returned incorrect number of results (expected %i, got %i)" \
                    % ( n_dups*2, len(mapped_reads) )

        locs = [ (
                    mapped_reads[i][2],
                    int(mapped_reads[i][3]),
                    mapped_reads[i][9],
                 )
                    for i in xrange(len(mapped_reads)) ]

        found_read = False
        # compare each mapped location to the original genome
        for loc in locs:

            # check to make sure we got the original read from fragments
            # if chr_name and start_bp match
            if truth[0] == loc[0] and truth[1] == loc[1]:
                # check original sequence
                if genome[loc[0]][truth[1]:truth[2]] == loc[2]:
                    found_read = True

            # compare mapped read back to original genome
            # need to mod by genome_len, since we took the original genome and duplicated it many times
            if genome[loc[0]][loc[1]%genome_len:loc[1]%genome_len+read_len].upper() != loc[2].upper():
                print "Truth  : ", truth
                print "Loc    : ", loc
                print "Genome : %s" % genome[loc[0]][loc[1]%genome_len:loc[1]%genome_len+read_len]
                print "Read   : %s" % loc[2]
                raise ValueError, \
                    "Mapped perfect read does not match sequence in reference genome."

        if found_read == False:
            print "Truth: ", truth
            print mapped_reads
            print "Mapped locations not equivalent - you probably failed to map the prior read."
            raise ValueError, \
                "Mapped locations at read id %i and Truth (%s, %i, %i) are not equivalent" \
                % ( int(mapped_reads[0][0]), truth[0], truth[1], truth[2] )

    sam_fp.close()

    # Cleanup the created files
    if CLEANUP:
        for fn in (read_fnames + output_filenames):
            os.remove(fn)
        shutil.rmtree(output_directory)

def test_diploid_genome():
    rls = [ 20, 50, 75 ]
    for rl in rls:
        genome, genome_fnames = build_diploid_genome( rl )
        map_diploid_genome( genome, genome_fnames, rl )
        print "PASS: Diploid genome Mapping %i BP Test. ( Statmap appears to be mapping diploid genomes correctly )" % rl

def test_diploid_genome_with_multiple_chrs():
    '''
    Make sure the machinery for handling multiple diploid chrs is working by
    mapping two separately generated diploid chrs
    '''
    rls = [ 20, 50, 75 ]
    for rl in rls:
        g1, of1 = build_diploid_genome( rl, "chr1" )
        g2, of2 = build_diploid_genome( rl, "chr2" )
        genome = dict( g1.items() + g2.items() )
        output_filenames = of1 + of2
        map_diploid_genome( genome, output_filenames, rl )
        print ( "PASS: Diploid genome with multiple chrs %i BP Test. "
                "( Statmap appears to be mapping diploid genomes with multiple "
                "chromosomes correctly )" % rl )

def test_lots_of_diploid_repeat_sequence_finding():
    rls = [ 25, ]
    n_dups = 4000
    genome_len=100

    for rl in rls:
        genome, output_filenames = build_diploid_genome(
                rl, gen_len=genome_len, n_dups=n_dups, n_mut=0 )
        map_duplicated_diploid_genome( genome, output_filenames, rl,
                genome_len=genome_len, n_dups=n_dups, nreads=100 )
        print (
"PASS: Highly repetitive diploid genome (pslocs integration)  %i BP Test. "
"( Statmap appears to be mapping highly repetitive diploid genomes correctly )"
            % rl )

def test_paired_end_diploid_repeat_sequence_finding( rl=20, n_dups=50 ):
    '''
    Make sure we are integrating diploid mapping, pslocs, and paired end reads properly
    '''
    output_directory = "smo_test_paired_end_diploid_repeat_sequence_finding_%i" % rl

    # build diploid genome
    genome, genome_fnames = build_diploid_genome( rl, gen_len=100,
            n_dups=n_dups, n_mut=0 )

    # build the paired end reads we want to sample whole chunks of the genome,
    # then build_read_from_fragments will take care of simulating paired end
    # reads
    assert 2*rl < 50 # make sure the fragments are long enough
    fragments = sample_uniformly_from_genome( genome, nsamples=100, frag_len=50)
    reads = build_reads_from_fragments(
        genome, fragments, read_len=rl, rev_comp=False, paired_end=True )

    # write the reads
    read_fnames = ( "tmp.1.fastq", "tmp.2.fastq" )
    reads_of_1 = open(read_fnames[0], "w")
    reads_of_2 = open(read_fnames[1], "w")
    build_paired_end_fastq_from_seqs( reads, reads_of_1, reads_of_2 )
    reads_of_1.close()
    reads_of_2.close()

    # map with statmap
    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len=rl )

    # test the sam file
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    # we divide by two since there are two lines for each paired end read
    total_num_reads = count_lines_in_sam(sam_fp) / 2

    # Every fragment will map twice in a perfect diploid genome (since the
    # paternal and maternal chromosomes are identical). Since we are simulating
    # duplicated reads by repeatedly copying the genome, each fragment will
    # also map an exponential number of times across each chromosome -
    # each first pair will map to all of the following second pairs for that
    # fragment, and so on.
    if len(fragments)*(n_dups**2)*2 != total_num_reads:
        raise ValueError, "Mapping returned the wrong number of reads (expected %i, got %i)" \
                % ( len(fragments)*(n_dups**2)*2, total_num_reads )

    # Cleanup the created files
    if CLEANUP:
        for fn in (read_fnames + output_filenames):
            os.remove(fn)
        shutil.rmtree(output_directory)

def test_multiple_indexable_subtemplates():
    rls = [ 50, 75 ]
    for rl in rls:
        test_sequence_finding( read_len=rl, indexed_seq_len=rl/2 )
        print "PASS: Multiple indexable subtemplates %i BP test." % rl
        # - 5 so the find_optimal_subseq_offset code should choose different offsets
        test_sequence_finding( read_len=rl, indexed_seq_len=rl/2 - 5 )
        print "PASS: Multiple indexable subtemplates (multiple offsets) %i BP test." % rl

def test_multiple_indexable_subtemplate_for_threep():
    rls = [ 50, 75 ]
    for rl in rls:
        test_sequence_finding( read_len=rl, indexed_seq_len=rl/2, rev_comp=True )
        print "PASS: Multiple indexable subtemplates 3' %i BP test." % rl

# note - for repeat sequence tests, n_dups must be greater than
# PSEUDO_LOC_MIN_SIZE in config.h to trigger the use of pseudo locations to store
# the locations

def test_multiple_indexable_subtemplates_for_repeat_sequences():
    rls = [ 50, 75 ]
    for rl in rls:
        test_duplicated_reads( read_len=rl,
                n_chrs=1, n_dups=100, gen_len=100, n_threads=-1, n_reads=100,
                indexed_seq_len = rl / 2  ) 
        print "PASS: Multiple indexable subtemplates in a highly repeated genome %i BP test" % rl

def test_multiple_indexable_subtemplates_for_diploid_mapping():
    n_dups = 4000
    genome_len=100

    rls = [ 50 ]
    for rl in rls:
        genome, output_filenames = build_diploid_genome(
                rl, gen_len=genome_len, n_dups=n_dups, n_mut=0 )
        map_duplicated_diploid_genome( genome, output_filenames, rl,
                genome_len=genome_len, n_dups=n_dups, nreads=100,
                indexed_seq_len=rl/2, num_threads=-1 )
        print "PASS: Multiple indexable subtemplates in a highly repeated diploid genome %i BP test " % rl

def sam_lines_match( line1, line2 ):
    """
    Compare two lines from a Statmap SAM file, controlling for case and
    numeric representation
    """
    for i, (f1, f2) in enumerate(
            izip( re.split( '\s', line1.strip() ),
                  re.split( '\s', line2.strip() ) )):
        # compare specific fields
        if i == 9: # SEQ field
            # control for case
            if f1.upper() != f2.upper():
                print "SEQ entries do not match"
                return False
        elif i == 12 or i == 13: # XP or XQ fields # control for numeric representation
            if float(f1.split(':')[-1]) != float(f2.split(':')[-1]):
                print "XP or XQ fields do not match"
                return False
        else:
            if f1 != f2:
                return False

    return True

def test_sam_output():
    try:
        import pysam
    except ImportError:
        print \
"WARNING : pysam not found - needed for optional test to check Statmap's SAM output. Skipping test_sam_output()"
        return

    output_directory = "smo_test_sam_output"

    n_reads=1000
    rl=20

    genome = build_random_genome( [1000, 1000], ["1", "2"] )

    fragments = sample_uniformly_from_genome( genome, nsamples=n_reads, frag_len=rl )
    reads = build_reads_from_fragments( 
        genome, fragments, read_len=rl, rev_comp=False, paired_end=False )

    ###### Write out the test files, and run statmap ###########################

    # write genome
    genome_fnames = ( "tmp.genome.fa", )
    with open( genome_fnames[0], "w") as genome_of:
        write_genome_to_fasta( genome, genome_of, 1 )
    
    # build and write the reads
    read_fnames = ( "tmp.fastq", )
    with open( read_fnames[0], "w" ) as reads_of:
        build_single_end_fastq_from_seqs( reads, reads_of )

    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len=rl )

    # Load the SAM file with pysam
    sam_fname = "./%s/mapped_reads.sam" % output_directory
    sam_copy_fname = sam_fname + ".copy"
    try:
        samfile = pysam.Samfile( sam_fname, "r" )
        samcopy = pysam.Samfile( sam_copy_fname, "wh", template=samfile )

        # copy reads from samfile into samcopy
        for read in samfile.fetch():
            samcopy.write(read)

        samfile.close()
        samcopy.close()
    except Exception, e:
        print e
        print "FAIL: Error parsing Statmap output SAM file with pysam"
        sys.exit(1)

    # compare the sam files
    sam = open( sam_fname )
    sam_copy = open( sam_copy_fname )

    for i, (l1, l2) in enumerate( izip( sam, sam_copy ) ):
        if not sam_lines_match( l1, l2 ):
            print "FAIL: Statmap output SAM file and pysam parsed copy do not match on line %i" \
                    % i
            print "Original   :", l1
            print "pysam copy :", l2
            sys.exit(1)

    sam.close()
    sam_copy.close()

    ###### Cleanup the created files ###############################################
    if CLEANUP:
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.fastq")
        shutil.rmtree(output_directory)

def test_softclipped_read_finding():
    rl = 20
    test_sequence_finding( rl, rev_comp = True, random_prefix_len=3 )
    print "PASS: Finding soft clipped reads"

def test_more_than_two_indexable_subtemplates():
    # no gaps
    test_sequence_finding( read_len=60, indexed_seq_len=20 )
    # gaps
    test_sequence_finding( read_len=80, indexed_seq_len=20 )
    print "PASS: Mapping perfect reads with more than two indexable subtemplates"

def test_overlapping_index_probes():
    pass

def main( RUN_SLOW_TESTS ):
    #print "Starting test_untemplated_g_finding()"
    #test_untemplated_g_finding()
    #sys.exit(1)
    print "Starting test_fivep_sequence_finding()"
    test_fivep_sequence_finding()
    print "Starting test_threep_sequence_finding()"
    test_threep_sequence_finding()
    #print "Start test_softclipped_read_finding()"
    #test_softclipped_read_finding()
    print "Starting test_paired_end_sequence_finding()"
    test_paired_end_sequence_finding( )
    print "Starting test_repeat_sequence_finding()"
    test_repeat_sequence_finding()
    print "Starting test_sam_output()"
    test_sam_output()
    print "Starting test_mutated_read_finding()"
    test_mutated_read_finding()
    print "Starting test_multithreaded_mapping()"
    test_multithreaded_mapping( )
    print "Starting test_multi_fasta_mapping()"
    test_multi_fasta_mapping( )
    print "Starting test_build_index()"
    test_build_index( )
    print "Starting test_iterative_mapping()"
    test_iterative_mapping()
    #print "Starting test_untemplated_g_finding()"
    #test_untemplated_g_finding()
    print "Starting test_diploid_genome()"
    test_diploid_genome()
    print "Starting test_diploid_genome_with_multiple_chrs()"
    test_diploid_genome_with_multiple_chrs()

    print "Start test_multiple_indexable_subtemplates()"
    test_multiple_indexable_subtemplates()
    print "Start test_more_than_two_indexable_subtemplates()"
    test_more_than_two_indexable_subtemplates()
    print "Start test_multiple_indexable_subtemplate_for_threep()"
    test_multiple_indexable_subtemplate_for_threep()
    print "Start test_multiple_indexable_subtemplates_for_repeat_sequences()"
    test_multiple_indexable_subtemplates_for_repeat_sequences()

    #print "Start test_multiple_indexable_subtemplates_for_diploid_mapping()"
    #test_multiple_indexable_subtemplates_for_diploid_mapping()
    return 

    if RUN_SLOW_TESTS:
        print "[SLOW] Starting test_lots_of_repeat_sequence_finding()"
        test_lots_of_repeat_sequence_finding( )
        print "[SLOW] Start test_lots_of_diploid_repeat_sequence_finding()"
        test_lots_of_diploid_repeat_sequence_finding()
        print "[SLOW] Start test_paired_end_diploid_repeat_sequence_finding()"
        test_paired_end_diploid_repeat_sequence_finding()

    #print "Starting test_index_probe()"
    #test_short_index_probe()

    # We skip this test because statmap can't currently
    # index reads less than 12 basepairs ( and it shouldn't: 
    #     we should be building a hash table for such reads )
    # test_short_sequences()

if __name__ == '__main__':
    try:
        main( RUN_SLOW_TESTS=True )
    except Exception, inst:
        if not P_STATMAP_OUTPUT:
            stderr.seek(0)
            print stderr.read()
        
        raise
