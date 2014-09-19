from tests import *

error_chars = '?adfh'
error_char_rates = [0.65, 0.03, 0.02, 0.01, 0.00]

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
    error_char_mappings = dict(zip(error_chars, error_char_rates ))
    
    rl = read_len = 200
    indexed_seq_len = 20
    nsamples = 10000
    output_directory = "smo_test_error_rate_estimation"
    
    ###### Prepare the data for the test #######################################
    # build a random genome
    r_genome = build_random_genome( [2000,2000], ["1","2"] )
    genome_fnames = ( "tmp.genome.fa", )
    with open( genome_fnames[0], "w" ) as genome_of:
        write_genome_to_fasta( r_genome, genome_of, 1 )
    
    # sample uniformly from the genome. This gives us the sequences
    # that we need to map. Note that we dont RC them, so every read should be in
    # the 5' direction
    fragments = sample_uniformly_from_genome( 
        r_genome, nsamples=nsamples, frag_len=rl )
    reads = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=True, paired_end=False )

    fragments = sample_uniformly_from_genome( 
        r_genome, nsamples=nsamples, frag_len=rl )
    reads2 = build_reads_from_fragments( 
        r_genome, fragments, read_len=rl, rev_comp=False, paired_end=False )

    new_genome = build_random_genome( [2000,2000], ["1","2"] )
    fragments = sample_uniformly_from_genome( 
        new_genome, nsamples=nsamples, frag_len=rl )
    reads3 = build_reads_from_fragments( 
        new_genome, fragments, read_len=rl, rev_comp=True, paired_end=False )
    
    def build_error_str(read):
        # return a complete shit read 10% of the time
        if random.random() < 0.10:
            return "?"*len(read)
        
        return "".join( random.choice(error_chars) 
                        for i in xrange(len(read)) )
        
    
    # mutate the reads according to some simple eror model
    def iter_mutated_reads( reads, base_error_rate=0.10 ):
        for read in reads:
            error_str = build_error_str(read)
            error_rates = [ base_error_rate # + 0.15*float(i)/len(read)
                            #+ error_char_mappings[char]
                            for i, char in enumerate( error_str ) ]
            
            mutated_read = []
            for base, error_prb in izip(read, error_rates):
                if random.random() < error_prb:
                    mutated_read.append( random.choice('acgtACGT') )
                else:
                    mutated_read.append( base )
            mutated_read = ''.join( mutated_read )
            yield ( mutated_read, error_str, read )

    mutated_reads = list( iter_mutated_reads( reads, 0.10 ) )
    mutated_reads += list(iter_mutated_reads(reads2, 0.10))
    mutated_reads += list(iter_mutated_reads(reads3, 0.01))
    
    # build and write the reads    
    reads_of = open("tmp.fastq", "w")
    build_single_end_fastq_from_mutated_reads( mutated_reads, reads_of )
    reads_of.close()
    
    ###### Write out the test files, and run statmap ###########################
    
    ## Map the data
    read_fnames = [ "tmp.fastq", ]
    map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len, search_type='e', num_threads=-1 )
    
    ###### Make sure that the error data looks correct #########################
    # records = load_error_data(os.path.join(output_directory, "error_stats.log"))
    
    ###### Cleanup the created files ###########################################
    if CLEANUP:
        os.remove("./tmp.genome.fa")
        os.remove("./tmp.fastq")
        shutil.rmtree(output_directory)


if __name__ == '__main__':
    print "Starting test_error_rate_estimation.py ..."

    try:
        test_error_rate_estimation()
    except Exception, inst:
        if not P_STATMAP_OUTPUT:
            stderr.seek(0)
            print stderr.read()
        
        raise

    print "PASS:    test error rate estimation"
