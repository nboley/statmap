import sys
import random

import tests as sc # for genome building and sampling functions

def splice_introns_from_sequence( sequence,
                                  n_introns,
                                  min_intron_len,
                                  max_intron_len
                                ):
    introns = []
    intron_starts = sorted(random.sample( xrange(len(sequence)), n_introns ))
    for intron_start in intron_starts:
        # pick a random intron length bounded by min and max intron_len
        intron_len = random.randint(min_intron_len, max_intron_len)
        introns.append( ( intron_start, intron_start + intron_len ) )

    exons = []
    print "INTRONS:" # DEBUG
    for i, intron in enumerate(introns):
        # DEBUG
        print "(%d, %d)" % (intron[0], intron[1])
        if i == 0: # first intron
            exons.append( (0, intron[0]) )
        else:
            exons.append( ( introns[i-1][1], intron[0] ) )

    # add last exon
    exons.append( ( introns[-1][1], len(sequence) ) )

    exon_seqs = []
    for exon in exons:
        exon_seqs.append( sequence[exon[0]:exon[1]] )

    return ''.join(exon_seqs)

def splice_introns_from_genome( genome,
                                n_introns = 1,
                                min_intron_len = 50,
                                max_intron_len = 100
                              ):
    transcriptome = {}
    for chromosome, sequence in genome.items():
        spliced_sequence = splice_introns_from_sequence( sequence,
                n_introns, min_intron_len, max_intron_len )
        transcriptome[chromosome] = spliced_sequence

    return transcriptome

def main():

    output_directory = "smo_rnaseq_sim"

    read_len = 50
    n_reads  = 100
    frag_len = 200

    genome = sc.build_random_genome( [1000,], ["1",] )
    transcriptome = splice_introns_from_genome( genome )

    # sample fragments from the transcriptome
    fragments = sc.sample_uniformly_from_genome(
            transcriptome, nsamples=n_reads, frag_len=frag_len )
    reads = sc.build_reads_from_fragments(
            transcriptome, fragments,
            read_len=read_len,
            rev_comp=False,
            paired_end=False )

    # Write out the test files
    with open("tmp.genome.fa", "w") as genome_fp:
        sc.write_genome_to_fasta( genome, genome_fp, 1 )

    with open("tmp.fastq", "w") as reads_fp:
        sc.build_single_end_fastq_from_seqs( reads, reads_fp,
                untemplated_gs_perc=0 )

    # Map the data with Statmap
    read_fnames = [ "tmp.fastq", ]
    genome_fnames = [ "tmp.genome.fa" ]
    sc.map_with_statmap( read_fnames, output_directory,
                         # change this to something less than half to test the 
                         # moving intron code
                         indexed_seq_len=read_len/2,
                         assay='r',
                         genome_fnames=genome_fnames )

    # test the sam file to make sure that each of the reads appears
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    total_num_reads = sum( 1 for line in sam_fp )
    sam_fp.seek(0)

    if len(fragments) > total_num_reads:
        raise ValueError, \
                "Mapping returned the wrong number of reads ( %i vs expected %i )." \
                % ( total_num_reads, len(fragments) )

if __name__ == "__main__": main()
