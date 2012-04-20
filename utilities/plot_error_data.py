#!/usr/bin/python

'''
Plot error data from FASTQ reads
'''

import sys

import matplotlib as mpl
#mpl.use('Agg') # use if X is not installed/for headless plot generation
import matplotlib.pyplot as plt

# Information about quality scores (Phred scores) for reads from different
# sequencers

#QUALITY_SCORE_RANGE = range()

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

def load_error_scores_from_file( read_fname ) :
    '''
    Given a fastq filename, return a list of the quality score strings for
    each read
    '''
    error_scores = []
    with open(read_fname) as fp:
        for i, line in enumerate(fp):
            if i%4 == 3:
                # be sure to remove newlines
                error_scores.append(line.strip())

    return error_scores

def process_error_scores( error_scores ):
    '''
    Loops through error_scores, computing file type, location error probs, and qual score error probs
    '''
    min_qual, max_qual = 255, 0

    # iterate over quality scores
    for es in error_scores:
        for c in es:
            max_qual = max( max_qual, ord(c) )
            min_qual = min( min_qual, ord(c) )

def categorize_files(args):
    files = { "sam": [], "fasta": [] }
    # Associate by filename extension (not perfect)
    for arg in args:
        if arg.endswith( ".sam" ):
            files["sam"].append(arg)
        elif arg.endswith( (".fa", ".fasta") ):
            files["fasta"].append(arg)
        else:
            print "Unrecognized file type: %s" % arg

    return files

def build_genome_from_fastas( fastas ):
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

def build_stats_from_sam_and_genome( sam, genome ):
    # set up data structures
    qual_score_cnts = [0]*255
    qual_score_mismatch_cnts = [0]*255
    # update length dynamically - should match the length of the longest read
    position_mismatch_cnts = []
    num_reads = 0

    sam_fp = open( sam )
    for line in sam_fp:
        fields = line.split()
        chr_name, start_bp, read, quality = (
            fields[2], int(fields[3]), fields[9], fields[10] )

        # if necessary, resize position_mismatch_cnts
        if len(read) > len(position_mismatch_cnts):
            for index in range( len(read) - len(position_mismatch_cnts) ):
                position_mismatch_cnts.append(0)

        for i, bp in enumerate(read):

            print "READ:   %s" % read
            print "GENOME: %s" % genome[chr_name][start_bp:start_bp+len(read)]
            # if bp in the read does not match sequence in the genome
            if( bp.upper() != genome[chr_name][i].upper() ):
                print "MISMATCH: %s, %s" % ( bp.upper(), genome[chr_name][i].upper() )
                qual_score_mismatch_cnts[ ord(bp) ] += 1
                position_mismatch_cnts[ i ] += 1

            qual_score_cnts[ ord(bp) ] += 1

        num_reads += 1

    print num_reads # may not be accurate for diploid - should do what iter_sam does in tests.py
    print qual_score_cnts
    print qual_score_mismatch_cnts
    print position_mismatch_cnts


class ErrorDataStruct():

    def __init__(self):
        self.num_unique_reads       = 0
        self.max_read_length        = 0
        self.loc_error_rates        = []
        self.qual_score_error_rates = []

    def __str__(self):
        rep = []
        rep.append("Num Unique Reads:\t%i" % self.num_unique_reads)
        rep.append("Max Read Length:\t%i" % self.max_read_length)
        rep.append("Loc Error Rates:")
        for i, ler in enumerate(self.loc_error_rates):
            rep.append("%i\t%i" % (i, ler))
        rep.append("Qual Score Error Rates:")
        for i, qser in enumerate(self.qual_score_error_rates):
            rep.append("%i\t%i" % (i, qser))

        return '\n'.join(rep)
        

class ErrorStatsLog():
    
    def __init__(self, filename):
        # init list of ErrorDataStruct
        self.data = []
        # open filename and load a list of ErrorDataStruct
        with open(filename) as fp:
            self._load_from_file( fp )

    def _load_from_file( self, fp ):
        '''
        Read error_stats.log, store each set of error data as an ErrorDataStruct
        '''
        index = -1
        curr_err_type = None
        for line in fp:
            if line.startswith("Num Unique Reads:"):
                # start of a new struct
                self.data.append( ErrorDataStruct() )
                index += 1
                # set num_unique_reads
                self.data[index].num_unique_reads = int(line.split()[-1].strip())
            elif line.startswith("Max Read Length:"):
                self.data[index].max_read_length = int(line.split()[-1].strip())
            elif line.startswith("Loc Error Rates:"):
                curr_err_type = "loc"
            elif line.startswith("Qual Score Error Rates:"):
                curr_err_type = "qual"
            else:
                if curr_err_type == "loc":
                    self.data[index].loc_error_rates.append(
                        float(line.split()[-1].strip()) )
                elif curr_err_type == "qual":
                    self.data[index].qual_score_error_rates.append(
                        float(line.split()[-1].strip()) )

    def __str__(self):
        return '\n'.join( map(str, self.data) ) 

    def plot(self):
        '''
        Plot two subplots for each ErrorDataStruct: one for
        loc_error_rates and one for qual_score_error_rates
        '''
        for i, eds in enumerate(self.data):
            # subplot indices are 1-based

            # loc_error_rates
            plt.subplot(len(self.data), 2, (i+1)*2-1)
            # plot style options: http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
            plt.plot( eds.loc_error_rates )
            plt.grid(True)
            plt.xlabel("Loc in read")
            plt.title("# Reads: %i" % eds.num_unique_reads)

            # qual_score_error_rates
            plt.subplot(len(self.data), 2, (i+1)*2)
            plt.plot( eds.qual_score_error_rates )
            plt.grid(True)
            plt.xlabel("Quality score value")

        plt.show()


def main():
    # load error_stats.log
    assert len(sys.argv) == 2
    esl = ErrorStatsLog( sys.argv[1] )
    esl.plot()

if __name__ == "__main__": main()
