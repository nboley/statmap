#!/usr/bin/python

'''
Plot error data from FASTQ reads
'''

import sys
import os

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
    '''
    A list of ErrorDataStructs parsed from error_stats.log
    '''
    
    def __init__(self, fp):
        self.data = []
        self.error_interval = 0

        # parse a list of ErrorDataStruct from error_stats.log
        self._load_from_file( fp )
        # set error_interval
        self.error_interval = self.data[0].num_unique_reads
        # remove unused error_stats

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

class ErrorPlot:

    def __init__(self, esl, mpdrds):
        self.current = 0
        self.esl = esl
        self.mpdrds = mpdrds

        # set up GUI callbacks
        fig = plt.figure()
        self.key_cid = fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.scroll_cid = fig.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.redraw()

    def redraw(self):
        # NOTE: matplotlib subplot indices are 1-based
        # NOTE: plot style options: http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
        edstruct = self.esl.data[self.current]

        plt.clf()
        #plt.title("# Reads: %i" % edstruct.num_unique_reads)
        print "# Reads: %i" % edstruct.num_unique_reads

        # Plot loc_error_rates
        plt.subplot(211)
        plt.plot( edstruct.loc_error_rates ) # hold=False
        plt.xlabel("Loc in read")
        plt.ylabel("P(mismatch)")
        #plt.ylim([0,1.0])

        # Plot qual_score_error_rates
        plt.subplot(212)
        plt.plot( edstruct.qual_score_error_rates )
        plt.xlabel("Quality Score")
        plt.ylabel("P(mismatch)")
        #plt.ylim([0,1.0])

    def on_key(self, event):
        #print 'you pressed', event.key
        if event.key == 'right' or event.key == 'up':
            self.current = min( self.current + 1, len(self.esl.data) - 1 )
        elif event.key == 'left' or event.key == 'down':
            self.current = max( self.current - 1, 0 )
        elif event.key == 'c':
            plt.clf()
        self.redraw()
    
    def on_scroll(self, event):
        #print 'you scrolled', event.button, event.step
        if event.button == 'left':
            self.current = min( self.current + 1, len(self.esl.data) - 1 )
        elif event.button == 'down':
            self.current = max( self.current - 1, 0 )
        self.redraw()


def usage():
    print "USAGE: ./plot_error_data.py statmap_output_directory/"
    sys.exit(1)

def main():

    if len(sys.argv) != 2: usage()

    # change into statmap output directory
    try:
        os.chdir(os.path.abspath(sys.argv[1]))
    except:
        print "Could not chdir to statmap_output_directory"
        usage()

    # load error stats and mapped reads
    with open("error_stats.log") as error_stats_fp:
        esl = ErrorStatsLog( error_stats_fp )
    with open("mapped_reads.sam") as mapped_reads_fp:
        mpdrds = MappedReads( mapped_reads_fp )

    # plot data interactively
    ep = ErrorPlot( esl, mpdrds )
    plt.show()

if __name__ == "__main__": main()
