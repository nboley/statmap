#!/usr/bin/python

import sys

def find_longest_matching_primer(read, primers, mismatches):
    '''
    Given a list of primers, find the longest one that matches the beginning of
    read.
    '''
    primers.sort(key=lambda s: len(s), reverse=True)
    for primer in primers:
        if match(primer, read, mismatches):
            return primer

    return ""   # empty string for no matching primer

def match(primer, read, mismatches):
    '''
    Returns True is primer matches the beginning of the read, allowing mismatches
    '''
    mm_sofar = 0
    for i, c in enumerate(primer):
        if primer[i] != read[i]:
            mm_sofar += 1
        if mm_sofar > mismatches:
            return False

    return True

def print_with_strict(s, trim, strict):
    '''
    Print s if conditions based on trim and strict pass.

    trim is the length of characters to trim; if it is 0, there was no match
    with a primer.

    If strict is set, we don't print reads that don't begin with a primer.
    '''
    if strict and trim == 0:
        pass
    else:
        sys.stdout.write(s)

def load_linkers_from_file(f):
    linkers = []
    with open( f) as fp:
        for line in fp:
            if line[0] == ">" or line[0] == ";": # FASTA
                pass # Comment
            else:
                linkers.append(line.strip())

    return linkers

def usage():
    print '''
    Usage: ./trim_linker_primers_from_reads.py linkers.fasta reads.fastq [ -m NUM -s ]

        linkers.fasta - list of linker primers to trim
        reads.fastq   - reads to trim
        -m NUM        - number of mismatches to allow; defaults to 1
        -s            - strict mode: throws out reads that DO NOT start with a
                        linker sequence. Not used by default.
    '''
    sys.exit(1)

if __name__ == '__main__':

    # Parse command line arguments
    if len( sys.argv ) not in (3,4,5,6):
        usage()

    linkers_fn = sys.argv[1]
    reads_fn = sys.argv[2]
    try:
        mismatches = int(sys.argv[sys.argv.index("-m")+1])
    except ValueError:
        mismatches = 1 # defaults to 1

    if "-s" in sys.argv:
        strict = True
    else:
        strict = False

    # Load linkers
    linkers = load_linkers_from_file(linkers_fn)

    # Process reads
    try:
        fp = open( reads_fn )
    except Exception:
        usage()

    trim = None
    for index, line in enumerate(fp):
        if index%4 == 1: # Sequence
            trim = len(find_longest_matching_primer(line, linkers, mismatches))
            print_with_strict( line[trim:], trim, strict )
        elif index%4 == 3: # Quality score
            print_with_strict( line[trim:], trim, strict )
        else:
            print_with_strict( line, trim, strict)
