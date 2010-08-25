import sys
import rpy

def usage():
    print "./generate_binom_p_values.py num_bs_samples"

if __name__ == "__main__":
    if len( sys.argv ) != 2:
        usage()
        sys.exit( 1 )
        
    num_bs_samples = int( sys.argv[1] )
    for loop in reversed( xrange( num_bs_samples + 1 ) ):
        p_value = rpy.r.binom_test( loop, num_bs_samples, alternative='greater' )['p.value']
        print """    case %i:
            return %e;""" % ( loop, p_value )
