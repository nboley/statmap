import sys
import numpy
from numpy import array
import pylab

class ErrorDataRecord( object ):
    """Store an error data record.
    

    """
    def __init__( self, nreads, min_readkey, max_readkey,
                        max_readlen, min_qualscore, max_qualscore ):
        self.nreads = nreads
        
        self.min_readkey = min_readkey
        self.max_readkey = max_readkey
        
        self.max_readlen = max_readlen
        
        self.min_qualscore = min_qualscore
        self.max_qualscore = max_qualscore
        
        self.cnts_array_dim = (max_readlen, max_qualscore-min_qualscore+1)
        self.cnts = numpy.zeros( self.cnts_array_dim )
        self.mm_cnts = numpy.zeros( self.cnts_array_dim )
        
        return
    
    def marginal_cnts( self ):
        rv = { 'pos': {}, 'qual': {} }
        
        rv['pos']['mm_cnts'] = self.mm_cnts.sum( axis = 1  )
        rv['pos']['cnts'] = self.cnts.sum( axis = 1  )
        
        rv['qual']['mm_cnts'] = self.mm_cnts.sum( axis = 0  )
        rv['qual']['cnts'] = self.cnts.sum( axis = 0  )
        
        return rv
    
    def build_flat_arrays( self ):
        def flatten_array( values ):
            res = []
            for pos in xrange( 1, self.max_readlen+1 ):
                for qual_score in xrange(
                        self.min_qualscore, self.max_qualscore+1 ):
                    value = values[pos-1, qual_score-self.min_qualscore]
                    res.append( (pos, qual_score, value) )
            return zip( *res )
        
        return map( flatten_array, (self.mm_cnts, self.cnts) )
        
        
    def plot_marginals( self ):
        marginals = self.marginal_cnts()['pos']
        pylab.plot( range(1,self.max_readlen+1), 
                    marginals['mm_cnts']/marginals['cnts'] )
        pylab.show()
        return


def build_error_data_record_from_log_string( log_str, header = None ):
    data = map( int, log_str.split() )
    meta_data = data[:6]
    record = ErrorDataRecord( *meta_data )

    data = data[6:]
    assert len( data )%2 == 0
    mismatches = data[:-len(data)/2]
    cnts = data[len(data)/2:]
    assert len( mismatches ) == len( cnts )
    
    if header != None:
        header = header.split()
    
    i = 0
    for qualscore in xrange( record.min_qualscore, 
                             record.max_qualscore+1 ):
        for read_pos in xrange( 1, record.max_readlen+1 ):
            # get the cnts from the flat array
            mm_cnt = mismatches[i]
            cnt = cnts[i]
            i += 1

            # add them to the error data record
            pos = ( read_pos-1, qualscore-record.min_qualscore )
            record.mm_cnts[ pos ] = mm_cnt
            record.cnts[ pos ] = cnt
                
    return record
    
def load_error_data( fname ):
    header = None
    records = []
    with open( fname ) as fp:
        for line_num, line in enumerate( fp ):
            # skip the header
            if line_num == 0: 
                header = line
                continue
            
            # add the record to the records list
            rec = build_error_data_record_from_log_string(line, header)
            records.append( rec )

    return records

from numpy import matrix

def main():
    records = load_error_data( sys.argv[1] )
    for record in records:
        mm_cnts, cnts = record.build_flat_arrays()
        freqs = array(mm_cnts[2])/(array(cnts[2])+0.01)
        from scipy.interpolate import LSQBivariateSpline
        res = LSQBivariateSpline( cnts[0], cnts[1], freqs, range(1,50,5), range(63, 104,5) )
        z = matrix( res.ev(cnts[1], cnts[2]), nrows=50 )
        print z
        break
        #
        freqs = (record.mm_cnts/(record.cnts+1)).T
        x = array(range( 1, record.max_readlen + 1 ))
        y = array(range( record.min_qualscore, record.max_qualscore + 1 ))
        print type(x)
        print type(y)
        print type(freqs)
        pylab.contour( x, y, freqs )
        pylab.show()
        #record.plot_marginals()
        break
        """
        for key, val in record.marginal_cnts().iteritems():
            print key
            print val
        print
        print
        """

if __name__ == '__main__':
    main()
