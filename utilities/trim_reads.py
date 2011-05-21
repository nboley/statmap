#!/usr/bin/python

import sys

# global trim paramaters
s_tl = 0
e_tl = 0

def usage():
    print "Usage: ./trim_reads.py reads.fastq beginning_trim_length [ end_trim_length ]"
    sys.exit(1)

if __name__ == '__main__':
    if len( sys.argv ) not in (3,4):
        usage()
    
    try:
        fp = open( sys.argv[1] )
        s_tl = int( sys.argv[2] )
        if len( sys.argv ) == 4:
            e_tl = int( sys.argv[3] )
    except Exception:
        usage()

    # account for the new line
    e_tl += 1
    
    for index, line in enumerate(fp):
        if index%2 == 1:
            sys.stdout.write( line[s_tl:-e_tl] + "\n" )
        else:
            sys.stdout.write( line )
