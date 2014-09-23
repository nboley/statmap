import sys
import numpy
from numpy import array, ndarray
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
        
        self.cnts_array_dim = (max_readlen, 1)
        self.base_type_cnts = numpy.zeros( self.cnts_array_dim )
        self.base_type_mm_cnts = numpy.zeros( self.cnts_array_dim )
        
        return
    
    def marginal_cnts( self ):
        rv = { 'pos': {}, 'qual': {} }
        
        rv['pos']['mm_cnts'] = self.mm_cnts.sum( axis = 1  )
        rv['pos']['cnts'] = self.cnts.sum( axis = 1  )
        
        rv['qual']['mm_cnts'] = self.mm_cnts.sum( axis = 0  )
        rv['qual']['cnts'] = self.cnts.sum( axis = 0  )
        
        return rv
    
    def build_flat_arrays( self, skip_zeros=False ):
        res = []
        for pos in xrange( 1, self.max_readlen+1 ):
            for qual_score in xrange(
                    self.min_qualscore, self.max_qualscore+1 ):
                
                cnt = self.cnts[pos-1, qual_score-self.min_qualscore]
                if skip_zeros and cnt == 0:
                    continue
                
                mm_cnt = self.mm_cnts[pos-1, qual_score-self.min_qualscore]

                res.append( (pos, qual_score, mm_cnt, cnt) )
        
        return [ numpy.array(x, dtype=float) for x in zip( *res ) ]
        
    def get_knots_with_obs( self ):
        return zip(*self.cnts.nonzero())
    
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


def main():
    records = load_error_data( sys.argv[1] )
    for record in records:
        poss, quals, mm_cnts, cnts = record.build_flat_arrays(skip_zeros=True)
        freqs = (mm_cnts+0.01)/cnts
        logit_freqs = numpy.log(freqs) - numpy.log(1-freqs)
        
        ### Build the spline basis
        from scipy.interpolate import LSQBivariateSpline, SmoothBivariateSpline
        weights = 1/numpy.sqrt(freqs*(1-freqs))
        basis = SmoothBivariateSpline(poss, quals, logit_freqs )
        fit_vals = basis.ev(poss, quals)
        fit_freqs = 1/( 1.0 + numpy.exp( -fit_vals ) )

        #pylab.plot( logit_freqs, fit_vals )
        #pylab.show( )
        print logit_freqs
        print fit_vals
        break
    
        x = array(range( 1, record.max_readlen + 1 ))
        y = array(range( record.min_qualscore, record.max_qualscore + 1 ))
        #z = ndarray( record.cnts_array_dim, buffer=freqs )

        #pylab.contour( y, x, z )
        #pylab.show()
        
        print ndarray( record.cnts_array_dim, buffer=freqs )
        print
        print z
        break
        #
        freqs = (record.mm_cnts/(record.cnts+1)).T
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
