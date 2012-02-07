#!/usr/bin/python

import sys

# Start by generating a list of possible permutations of known primers

# 5'Linker (#4)
KNOWN_PRIMERS = (
    'ACACAGCAG',  # 5'L-N6-Phos_#4
    'ACACAGCAGG', # 5'L-GN5-Phos_#4
    'CTGCTGAGT',  # 5'L-low_#4
    )

BASES = ('A', 'T', 'C', 'G')

# Consider line to match even if 1 character is off
# Any mismatch - configurable to just purine>purine, pyrimidine>pyrimidine
# Permute KNOWN_PRIMERS to compare against line
POSSIBLES = []
for primer in KNOWN_PRIMERS:
    for i in xrange(len(primer)):
        POSSIBLES.extend([primer[:i] + B + primer[i+1:] for B in BASES])

def trim_length(line):
    '''
    Return number of characters to trim from the start of a line
    Matches start of line against known linker primers, with small margin of error (1 chr)
    Returns 0 if no match
    '''
    # Match beginning of line against possibles
    # Lazy evaluation - takes the first match it finds
    for possible in POSSIBLES:
        if possible == line[:len(possible)]:
            return len(possible)
    
    # no match - trim 0 characters
    return 0

def usage():
    print "Usage: ./trim_linker_primers_from_reads.py reads.fastq"
    sys.exit(1)

if __name__ == '__main__':
    if len( sys.argv ) not in (2,):
        usage()

    try:
        fp = open( sys.argv[1] )
    except Exception:
        usage()

    trim_l = None
    for index, line in enumerate(fp):
        if index%4 == 1:
            trim_l = trim_length(line.upper())
            #print "TRIM LENGTH: %s" % trim_l
            sys.stdout.write( line[trim_l:] )
        elif index%4 == 3:
            sys.stdout.write( line[trim_l:] )
        else:
            sys.stdout.write( line )
