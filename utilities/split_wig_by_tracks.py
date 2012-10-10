#!/usr/bin/env python

import sys

track_name = sys.argv[1]

fp = None
fps = [ (open( track_name + ".plus.wig", "w" ), "plus"),
        (open( track_name + ".minus.wig", "w" ), "minus") ]
for line in sys.stdin:
    if line.startswith( 'track' ):
        if fp != None:
            fp.close()
        fp, strand = fps.pop(0)
        fp.write( "track type=wiggle_0 name=%s_%s\n" % (track_name, strand) )
        continue
    fp.write( line )

fp.close()
