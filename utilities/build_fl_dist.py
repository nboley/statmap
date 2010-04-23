#!/usr/bin/python

import sys

import numpy
from scipy.stats import distributions

def usage():
    print """./build_fl_dist.py mean std

Builds a fragment length distribution file with a normal distribution, truncated at +- 2 std.
"""

def generate_norm_density( mean, std ):
    density = numpy.array( [distributions.norm( mean, std ).pdf(loop) for loop in xrange( max(0,mean-2*std), mean+2*std+1 ) ] )
    density = density/density.sum()
    for loop, den in zip(xrange( max(0,mean-2*std), mean+2*std+1 ), density ):
        print "%i\t%f" % ( loop, den )

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
    else:
        mean = int(sys.argv[1])
        std = int(sys.argv[2])
        generate_norm_density( mean, std )
