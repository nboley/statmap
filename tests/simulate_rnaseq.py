import sys
import random
from itertools import izip
import re

import tests as sc # for genome building and sampling functions

MAX_ALLOWABLE_MISMATCHES = 3

def splice_introns_from_sequence( sequence,
                                  n_introns,
                                  min_intron_len,
                                  max_intron_len
                                ):
    # TODO revise intron sampling code - atm it's possible to pick introns
    # beyond the boundary of the sequence
    introns = []
    intron_starts = sorted(random.sample( xrange(len(sequence)), n_introns ))
    for intron_start in intron_starts:
        # pick a random intron length bounded by min and max intron_len
        intron_len = random.randint(min_intron_len, max_intron_len)
        introns.append( ( intron_start, intron_start + intron_len ) )

    # DEBUG
    #print "INTRONS:"

    exons = []
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

    return introns, ''.join(exon_seqs)

def splice_introns_from_genome( genome,
                                n_introns = 1,
                                min_intron_len = 10,
                                max_intron_len = 50
                              ):
    transcriptome = {}
    for chromosome, sequence in genome.items():
        introns, spliced_sequence = splice_introns_from_sequence( sequence,
                n_introns, min_intron_len, max_intron_len )
        transcriptome[chromosome] = spliced_sequence

    return introns, transcriptome

def num_mismatches( str1, str2 ):
    num_mismatches = 0
    for c1, c2 in izip( str1, str2 ):
        if c1 != c2:
            num_mismatches += 1

    return num_mismatches

def check_mapped_read( mapped_read, fragments, genome ):
    # find the corresponding fragment
    read_id = int(mapped_read[0])
    truth = fragments[read_id]

    # split the cigar string into entries
    cigar_re = r'([0-9]+)([MN])'
    cigar_entries = re.findall( cigar_re, mapped_read[5] )

    start_pos = int(mapped_read[3])
    contig = mapped_read[2]

    # Make sure the start and contig of the mapped read match the truth
    # Note - we can't compare the start locations, since we sampled fragments
    # from the transcriptome, so all of the start_bps after any introns are
    # offset. TODO - fix this? Or just use the genome comparison as the test?
    #if( truth[0] != contig or truth[1] != start_pos ):
    #    print "Mapped read does not match truth: found (%s, %i) for (%s, %i)" \
    #            % ( contig, start_pos, truth[0], truth[1] )
    #    print mapped_read
    #    print truth
    #    sys.exit(1)

    flag = int(mapped_read[1])
    rev_comp = False
    if flag & 0x10 > 0:
        rev_comp = True

    read_seq = mapped_read[9]

    # Keep track of the location in the genome to compare to matched sequences
    # (accounting for gaps)
    genome_pos = start_pos
    seq_pos = 0

    for entry in cigar_entries:
        entry_len = int( entry[0] )

        if entry[1] == 'M':
            # Check the matched sequence against the original genome
            # TODO rev comp if necessary
            read_segment = read_seq[seq_pos:seq_pos+entry_len]
            genome_segment = genome[truth[0]][genome_pos:genome_pos+entry_len]

            # Since the read might overlap the intron slightly and still have
            # been within Statmap's bounds, we say a match is good enough if it
            # has less than some number of mismatches
            mismatch_count = num_mismatches( read_segment.upper(),
                                             genome_segment.upper() )
            if( mismatch_count > MAX_ALLOWABLE_MISMATCHES ):
                print "Mapped read match segment does not match genome"
                print "read_id:", read_id
                print "had %i mismatches to the reference genome" % mismatch_count
                print read_segment
                print genome_segment
                sys.exit(1)

            # update the position in the underlying mapped read_len
            seq_pos += entry_len

        # Update the position in the genome to compare to
        if rev_comp:
            genome_pos -= entry_len
        else:
            genome_pos += entry_len

def main():

    output_directory = "smo_rnaseq_sim"

    read_len = 100
    n_reads  = 100
    frag_len = 200

    genome = sc.build_random_genome( [1000,], ["1",] )
    introns, transcriptome = splice_introns_from_genome( genome )

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
                         indexed_seq_len=20,
                         assay='r',
                         genome_fnames=genome_fnames )

    # test the sam file to make sure that each of the reads appears
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    total_num_reads = sum( 1 for mapped_read in sc.iter_sam_reads(sam_fp) )
    sam_fp.seek(0)

    # nonmapping read ids
    nonmapping_reads_fname = "./%s/mapped_reads.db.nonmapping" % output_directory
    nonmapping_read_ids = []
    with open( nonmapping_reads_fname ) as nonmapping_fp:
        for line in nonmapping_fp:
            read_id = int(line.strip())
            nonmapping_read_ids.append(read_id)

    if len(fragments) > total_num_reads:
        print "Mapped (%i/%i) reads" % ( total_num_reads, len(fragments) )

        print "Nonmapping reads (probably had index probes that overlapped the intron):"
        for read_id in nonmapping_read_ids:
            print "Read_id %i : " % read_id, fragments[read_id]

    for mapped_reads in sc.iter_sam_reads(sam_fp):
        for mapped_read in mapped_reads:
            check_mapped_read( mapped_read, fragments, genome )

if __name__ == "__main__": main()
