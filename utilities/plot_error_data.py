#!/usr/bin/python

'''
Plot error data from FASTQ reads
'''

import sys
import os
from subprocess import call

from itertools import izip

import matplotlib as mpl
#mpl.use('Agg') # use if X is not installed/for headless plot generation
import matplotlib.pyplot as plt

# Information about quality scores (Phred scores) for reads from different
# sequencers

FORMAT_TYPES = [
    { 'QUAL_SHIFT': 33, 'ARE_LOG_ODDS': False },
    { 'QUAL_SHIFT': 64, 'ARE_LOG_ODDS': False },
    { 'QUAL_SHIFT': 59, 'ARE_LOG_ODDS': True },
    { 'QUAL_SHIFT': 50, 'ARE_LOG_ODDS': False },
]

FILE_TYPES = {
    'SANGER_FQ':            FORMAT_TYPES[0],
    'ILLUMINA_v13_FQ':      FORMAT_TYPES[1],
    'SOLEXA_v14_FQ':        FORMAT_TYPES[1],
    'ILLUMINA_v15_FQ':      FORMAT_TYPES[1],
    'SOLEXA_LOG_ODDS_FQ':   FORMAT_TYPES[2],
    'TEST_SUITE_FORMAT':    FORMAT_TYPES[1],
    'MARKS_SOLEXA':         FORMAT_TYPES[3],
}

def usage():
    print "USAGE: ./plot_error_data.py statmap_output_directory"; sys.exit(1)

def determine_read_file_type_from_error_scores(error_scores):
    '''
    error_scores - a list of read quality score strings

    Return the filetype
    '''
    min_qual, max_qual = 255, 0

    # compute max and min quality scores over all reads
    for es in error_scores:
        for c in es:
            max_qual = max( max_qual, ord(c) )
            min_qual = min( min_qual, ord(c) )

    print "Calculated max and min quality scores as {} and {}".format(
            max_qual, min_qual )

    # Set filetype - this logic is from statmap.c:164
    filetype = None
    if max_qual == 73:
        filetype = "SANGER_FQ"
    elif max_qual > 73 and max_qual <= 90:
        filetype = "MARKS_SOLEXA"
    else:
        if max_qual > 73:
            if max_qual < 104:
                print "WARNING : maximum input score not achieved ({})".format(
                        max_qual )
            if min_qual < 64:
                filetype = "SOLEXA_LOG_ODDS_FQ"
            elif min_qual < 66:
                filetype = "ILLUMINA_v13_FQ"
            elif min_qual < 70:
                filetype = "ILLUMINA_v15_FQ"
            elif min_qual < 90:
                filetype = "MARKS_SOLEXA"
            elif min_qual == 104:
                filetype = "TEST_SUITE_FORMAT"
            else:
                print "Could not automatically determine filetype. Exiting."
                sys.exit(1)

    print "Automatically determined filetype: {}".format( filetype )
    return filetype

def build_genome_from_fastas( fastas ):
    '''
    Given a list of fasta files, build a dictionary of chr: sequence
    '''
    # genome is a dictionary of chr_name : sequence
    genome = {}
    chrid = None

    for fasta in fastas:
        fp = open( fasta )
        for line in fp:
            if line.startswith(">"): # chromosome identifier
                # start new chr
                chrid = line.split(">")[1].strip()
                genome[chrid] = []
            else:
                # add line of fasta to list of lines
                genome[chrid].append( line.strip() )

    # join lists of strings into strings
    # (this is faster than repeatedly appending to strings)
    for k, v in genome.iteritems():
        genome[k] = ''.join(v)

    return genome

class ErrorDataStruct():
    '''
    Stores error information from a single error_data struct
    '''

    def __init__(self, fp=None):

        # initialize instance vars
        self.num_unique_reads           = 0
        self.max_read_length            = 0
        self.position_mismatch_cnts     = []
        self.qual_score_cnts            = []
        self.qual_score_mismatch_cnts   = []

        if fp:
            self._load_from_fp(fp)

    def _load_from_fp(self, fp):
        """
        Load error data struct given a fp that is on the first line of an
        error data struct description
        """
        # make sure we are at the start of a new struct
        line = fp.readline().strip() # remove newline
        if line.startswith("num_unique_reads"):
            self.num_unique_reads = int(line.split()[1])
        else:
            raise RuntimeError("Error Log fp misaligned")

        # read the next line for the max_read_length
        line = fp.readline().strip()
        self.max_read_length = int(line.split()[1])

        # read the position_mismatch_cnts
        for i in range(self.max_read_length):
            line = fp.readline().strip()
            self.position_mismatch_cnts.append(
                    float(line.split()[1])
                )

        # read the quality score counts
        max_num_qual_scores = 256
        for i in range(max_num_qual_scores):
            line = fp.readline().strip()
            self.qual_score_cnts.append(
                    float(line.split()[2])
                )
            self.qual_score_mismatch_cnts.append(
                    float(line.split()[1])
                )

class ErrorStatsLog():
    '''
    A list of ErrorDataStructs parsed from error_stats.log
    '''
    
    def __init__(self, fp):
        # parse a list of ErrorDataStructs from error_stats.log
        self.data = []
        self._load_from_file(fp)

    def _load_from_file(self, fp):
        '''
        Read error_stats.log, store each set of error data as an ErrorDataStruct
        '''
        # loop over file, loading structs in one at a time
        # after every struct, read the next line in to test for EOF
        while True:
            pos = fp.tell()
            line = fp.readline() # returns "" on EOF
            if line:
                # return to original position
                fp.seek(pos)
                self.data.append( ErrorDataStruct(fp) )
            else:
                break

class MappedReads:
    '''
    Store relevant information about mapped reads so we can compare read
    quality scores with our error model
    '''

    def __init__(self, fp):
        self.reads = []
        # Load mapped reads from SAM
        for line in fp:
            # store chr, start_bp, read, quality score
            fields = line.split()
            self.reads.append(
                (
                    fields[2], fields[3],   # chr, start_bp
                    fields[9], fields[10]   # read, quality score
                )
            )

def plot_error_stats( esl, mpdrds ):
    for ed in esl.data:
        # Compute loc_error_rates
        loc_error_rates = [
                (pmc / ed.num_unique_reads)
                for pmc in ed.position_mismatch_cnts
            ]
        # Plot loc_error_rates
        plt.subplot(211)
        plt.plot( loc_error_rates )
        plt.xlabel("Loc in read"); plt.ylabel("P(mismatch)")

        # Compute qual_score error rates
        qual_score_error_rates = [
                (qsm / qsc ) if qsc != 0 else 0
                for qsm, qsc in
                izip( ed.qual_score_mismatch_cnts,
                      ed.qual_score_cnts )
            ]
        # Plot qual_score_error_rates
        plt.subplot(212)
        plt.plot( qual_score_error_rates )
        plt.xlabel("Quality Score"); plt.ylabel("P(mismatch)")

def main():
    if len(sys.argv) != 2: usage()

    # change into statmap output directory
    try:
        os.chdir(os.path.abspath(sys.argv[1]))
    except:
        print "Could not chdir to Statmap output directory %s" % \
                os.path.abspath(sys.argv[1])
        usage()

    # load error stats and mapped reads
    with open("error_stats.log") as error_stats_fp:
        esl = ErrorStatsLog( error_stats_fp )

    # check to see if mapped_reads.sam exists; if not, create it
    if "mapped_reads.sam" not in os.listdir("."):
        print "mapped_reads.sam not found, building it..."
        buildsam_cmd = "%s %s > %s" % (
                "~/statmap/bin/mapped_reads_into_sam",
                ".",
                "mapped_reads.sam"
            )
        rv = call( buildsam_cmd, shell=True )
        if rv != 0:
            raise RuntimeError("Could not build SAM file from mapped reads")

    with open("mapped_reads.sam") as mapped_reads_fp:
        mpdrds = MappedReads( mapped_reads_fp )

    # plot data interactively
    plot_error_stats( esl, mpdrds )
    plt.show()

if __name__ == "__main__": main()
