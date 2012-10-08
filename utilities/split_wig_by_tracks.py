#!/usr/bin/env python

import sys

track_name = sys.argv[1]

fp = None
fps = [ open( track_name + ".plus.wig", "w" ), open( track_name + ".minus.wig", "w" ) ]
for line in sys.stdin:
    if line.startswith( 'track' ):
        if fp != None:
            fp.close()
        fp = fps.pop(0)
        fp.write( "track type=wiggle_0 name=%s_plus\n" % track_name )
        continue
    fp.write( line )

fp.close()
