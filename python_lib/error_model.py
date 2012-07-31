import sys
import numpy

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
        rv = {}
        
        rv['pos_mm_cnts'] = self.mm_cnts.sum( axis = 1  )
        rv['pos_cnts'] = self.cnts.sum( axis = 1  )
        
        rv['qual_mm_cnts'] = self.mm_cnts.sum( axis = 0  )
        rv['qual_cnts'] = self.cnts.sum( axis = 0  )
        
        return rv
    
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
            if header != None:
                print i, read_pos, qualscore, header[i+6], header[i+len(cnts)+6]
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

if __name__ == '__main__':
    main()
