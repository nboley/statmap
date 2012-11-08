#!/usr/bin/env python
"""
Merge multiple wiggles, summing values from the same index and thresholding by
1.

For now, this can only handle wiggles that use variableStep formatting

Wiggle spec: https://cgwb.nci.nih.gov/goldenPath/help/wiggle.html
"""

import sys
import os
import re
import operator

from collections import defaultdict

# \S matches any non-whitespace character
contig_re = r'variableStep chrom=(?P<chrom>\S+)' 

FILTER_THRESHOLD = 1

def merge_wiggle( wig_fname, merged_dict ):
    wig_fp = open( wig_fname )

    current_contig = None
    for line in wig_fp:
        if line.startswith("track"):
            pass
        elif line.startswith("variableStep"):
            current_contig = re.match( contig_re, line ).group('chrom')
        elif line.startswith("fixedStep"):
            raise NotImplementedError, \
"This script doesn't support merging wiggles with fixedStep formatting (yet)."
        else:
            # parse the entry into pos and dataValue
            pos = int(line.split()[0])
            dataValue = float(line.split()[1])
            merged_dict[current_contig][pos] += dataValue

    wig_fp.close()

def output_merged_wiggle( out_fname, merged_dict ):
    out_fp = open( out_fname, "w" )

    # track line is just the output filename (we don't try anything clever)
    track_line = "track type=wiggle_0 name=%s" \
            % os.path.basename(out_fname)
    out_fp.write( track_line + "\n" )

    # Don't bother keeping chrs in their input order - doesn't matter (may be
    # different across input wigs, don't assume)
    for contig, entries in merged_dict.iteritems():
        # write declaration line for each contig
        out_fp.write("variableStep chrom=%s\n" % contig)
        # write data lines
        sorted_entries = sorted(entries.iteritems(), key=operator.itemgetter(0))
        for pos, value in sorted_entries:
            # only print entry in value above the threshold
            if value >= FILTER_THRESHOLD:
                out_fp.write( "%i\t%e\n" % (pos, value) )

    out_fp.close()

def main():
    if len(sys.argv) < 4:
        print "Usage: ./merge_wiggles.py merged.wig 1.wig 2.wig [3.wig ...]"
        sys.exit(1)

    merged_wig_fname = sys.argv[1]

    # For each chromsome, store the merged values for each chromosome position
    # merged_dict = {
    #     'chr1': { 1: 0.1, 2: 0.1, ... },
    #     ...
    # }
    merged_dict = defaultdict( lambda: defaultdict(int) )

    for wig_fname in sys.argv[2:]:
        merge_wiggle( wig_fname, merged_dict )

    output_merged_wiggle( merged_wig_fname, merged_dict )

if __name__ == "__main__":
    main()
