from tests import *

def test_error_rate_estimation( ):
    """Test to ensure that the error modelling code is working properly.
    
    We use the following error model:
    
    Error rate from position:
        0.01 + (pos/read_len)*0.04
        
    Error rate from char string:
        ?: +0.65
        a: +0.03
        d: +0.02
        f: +0.01
        h: +0.00
    """
    error_char_mappings = dict(zip('?adfh', [0.65, 0.03, 0.02, 0.01, 0.00] ))
    
    rl = read_len = 50
    indexed_seq_len = 20
    nsamples = 19900
    output_directory = "smo_test_error_rate_estimation"
    
    
    ###### Prepare the data for the test #######################################
    # build a random genome
    genome_of = open("tmp.genome.fa", "w")
    r_genome = build_random_genome( [2000,2000], ["1","2"] )
    write_genome_to_fasta( r_genome, genome_of, 1 )
    genome_of.close()
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in
    # the 5' direction
    fragments = sample_uniformly_from_genome( 
        r_genome, nsamples=nsamples, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )
    
    def build_error_str(read):
        # return a complete shit read 10% of the time
        if random.random() < 0.10:
            return "?"*len(read)
        
        return "".join( random.choice('?adfh') #'?adfh') 
                        for i in xrange(len(read)) )
        
    
    # mutate the reads according to some simple eror model
    def iter_mutated_reads( reads ):
        for read in reads:
            error_str = build_error_str(read)
            error_rates = [ 0.01 + 0.03*float(i)/len(read)
                            + error_char_mappings[char]
                            for i, char in enumerate( error_str ) ]
            
            mutated_read = []
            for base, error_prb in izip(read, error_rates):
                if random.random() < error_prb:
                    mutated_read.append( random.choice('acgtACGT') )
                else:
                    mutated_read.append( base )
            mutated_read = ''.join( mutated_read )
            yield ( mutated_read, error_str, read )

    mutated_reads = list( iter_mutated_reads( reads ) )

    # build and write the reads    
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_mutated_reads( mutated_reads, reads_of )
    reads_of.close()
    
    ###### Write out the test files, and run statmap ###########################
    
    ## Map the data
    read_fnames = [ "tmp.fastq", ]
    map_with_statmap( read_fnames, output_directory, indexed_seq_len, 
                      search_type='e', num_threads=8 )
    
    ###### Make sure that the error data looks correct #########################
    records = load_error_data(os.path.join(output_directory, "error_stats.log"))
    
    ###### Cleanup the created files ###########################################
    if CLEANUP:
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.fastq")
        shutil.rmtree(output_directory)


if __name__ == '__main__':
    try:
        test_error_rate_estimation()
    except Exception, inst:
        if not P_STATMAP_OUTPUT:
            stderr.seek(0)
            print stderr.read()
        
        raise
