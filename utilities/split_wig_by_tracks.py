import sys

track_name = sys.argv[1]

fp = None
fps = [ open( track_name + ".plus.wig", "w" ), open( track_name + ".minus.wig", "w" ) ]
strands = [ "plus", "minus" ]
for line in sys.stdin:
    if line.startswith( 'track' ):
        if fp != None:
            fp.close()
        fp = fps.pop(0)
        fp.write( "track type=wiggle_0 name=%s_%s\n" % (track_name, strands.pop(0) ) )
        continue
    fp.write( line )

fp.close()
