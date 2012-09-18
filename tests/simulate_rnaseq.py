import sys
import random
from itertools import izip
import re

import tests as sc # for genome building and sampling functions

def splice_introns_from_sequence( sequence,
                                  n_introns,
                                  min_intron_len,
                                  max_intron_len
                                ):
    introns = []
    # Subtract max_intron_len to account for the missing sample range at the end
    intron_starts = sorted(
            random.sample( xrange(len(sequence) - max_intron_len), n_introns ))
    for intron_start in intron_starts:
        intron_len = random.randint(min_intron_len, max_intron_len)
        introns.append( ( intron_start, intron_start + intron_len ) )

    # Now build the corresponding set of exons
    exons = []
    for i, intron in enumerate(introns):
        if i == 0: # first intron
            exons.append( (0, intron[0]) )
        else:
            exons.append( ( introns[i-1][1], intron[0] ) )

    # add last exon
    exons.append( ( introns[-1][1], len(sequence) ) )

    return introns, ''.join( [ sequence[exon[0]:exon[1]] for exon in exons ] )

def splice_introns_from_genome( genome,
                                n_introns = 1,
                                min_intron_len = 10,
                                max_intron_len = 50,
                              ):
    assert min_intron_len <= max_intron_len

    transcriptome = {}
    introns_dict = {}

    for chromosome, sequence in genome.items():
        introns, spliced_sequence = splice_introns_from_sequence( sequence,
                n_introns, min_intron_len, max_intron_len )
        transcriptome[chromosome] = spliced_sequence
        introns_dict[chromosome] = introns

    return introns_dict, transcriptome

def find_fragment_pos_in_genome( fragment, introns ):
    offset = 0
    for intron in introns[fragment[0]]: # sorted by start
        if intron[0] < fragment[1]:
            offset += ( intron[1] - intron[0] )

    return fragment[1] + offset, fragment[2] + offset

def categorize_fragments( genome, transcriptome, introns, fragments,
        read_len, indexed_seq_len ):

    categorized_fragments = {
            "no_overlap": [],
            "intron_overlap": [],
            "intron_overlap_probe": []
        }

    # A mapped fragment's read_id is the index of the fragment in fragments,
    # due to build_single_end_fastq_from_seqs
    for read_id, fragment in enumerate(fragments):
        # Fragment (start, stop) coordinates are in the transcriptome
        # Intron (start, stop) coordinates are in the genome

        # Convert fragment coordinates to coordinates in the genome
        frag_start, frag_stop = find_fragment_pos_in_genome( fragment, introns )
        #print "original: (%i, %i)" % ( fragment[1], fragment[2] )
        #print "genome  : (%i, %i)" % ( frag_start, frag_stop )

        # check against each intron
        for intron in introns[fragment[0]]:
            # a read's start cannot be inside an intron, since we sampled reads
            # from the transcriptome. a fragment overlaps the intron if the
            # start is < the intron's start and the end is >= the intron start
            if frag_start < intron[0] and frag_stop >= intron[0]:
                # if the range between the fragment start + indexed_seq_len or
                # the range from fragment end - indexed_seq_len is inside the
                # intron, then the index probe overlaps the intron
                if ( frag_start + indexed_seq_len > intron[0] or
                     frag_stop - indexed_seq_len < intron[1] ):
                    categorized_fragments["intron_overlap_probe"].append( read_id )
                else:
                    categorized_fragments["intron_overlap"].append( read_id )
            else:
                categorized_fragments["no_overlap"].append( read_id )

    return categorized_fragments

def num_mismatches( str1, str2 ):
    num_mismatches = 0
    for c1, c2 in izip( str1, str2 ):
        if c1 != c2:
            num_mismatches += 1

    return num_mismatches

def check_mapped_read_sequence( mapped_read, truth, genome,
        num_allowed_mismatches ):
    read_id = int(mapped_read[0])
    start_pos = int(mapped_read[3])

    # parse the cigar string into a list of entries
    cigar_re = r'([0-9]+)([MN])'
    cigar_entries = re.findall( cigar_re, mapped_read[5] )

    read_seq = mapped_read[9]

    flag = int(mapped_read[1])
    rev_comp = False
    if flag & 0x10 > 0:
        rev_comp = True

    # Keep track of the location in the genome to compare to matched sequences
    # (accounting for gaps)
    genome_pos = start_pos
    seq_pos = 0

    for entry in cigar_entries:
        entry_len = int( entry[0] )

        if entry[1] == 'M':
            # TODO rev comp if necessary
            read_segment = read_seq[seq_pos:seq_pos+entry_len]
            genome_segment = genome[truth[0]][genome_pos:genome_pos+entry_len]

            # Since the read might overlap the intron slightly and still have
            # been within Statmap's bounds, we say a match is good enough if it
            # has less than some number of mismatches
            mismatch_count = num_mismatches( read_segment.upper(),
                                             genome_segment.upper() )
            if( mismatch_count > num_allowed_mismatches ):
                print "Mapped read match segment does not match genome"
                print "read_id %i had %i mismatches to the reference genome" \
                        % (read_id, mismatch_count)
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

def check_mapped_read( mapped_read, fragments, categorized_fragments, introns,
        genome ):
    # find the source fragment
    read_id = int(mapped_read[0])
    truth = fragments[read_id]

    start_pos = int(mapped_read[3])
    contig = mapped_read[2]

    ### Check - mapped reads must have mapped to their known location in the genome

    # Make sure the start and contig of the mapped read match the truth
    # The fragment's start is in the transcriptome, so we need to find the
    # correct location for comparison to the mapped read's start position
    frag_start_in_genome, frag_stop_in_genome = find_fragment_pos_in_genome(
            truth, introns )
    if( truth[0] != contig or frag_start_in_genome != start_pos ):
        print "Mapped read does not match truth: found (%s, %i) for (%s, %i)" \
                % ( contig, start_pos, truth[0], frag_start_in_genome )
        print "Introns:", introns
        print "Mapped: ", mapped_read
        print "Truth:  ", truth
        sys.exit(1)

    if read_id in categorized_fragments['no_overlap']:
        # Check - if mapped read didn't overlap the intron, its sequence should
        # match the genome exactly
        check_mapped_read_sequence( mapped_read, truth, genome, 0 )
    elif read_id in categorized_fragments['intron_overlap']:
        # Check - if the mapped read overlapped the intron (but the index
        # probes did not), we should have mapped it perfectly (?)
        check_mapped_read_sequence( mapped_read, truth, genome, 3 )
    elif read_id in categorized_fragments['intron_overlap_probe']:
        # Check - if the mapped read had probes overlapping the intron (by just
        # a little bit), there may be some noise but Statmap would still have
        # mapped it. We say 3 mismatches are allowed, which empirically has been reasonable.
        check_mapped_read_sequence( mapped_read, truth, genome, 3 )
    else:
        print "ERROR       :  read id %i was not categorized. Fragment categorization code in the unit test failed." \
                % read_id
        sys.exit(1)

def check_statmap_output( output_directory, genome,
        transcriptome, introns, fragments, categorized_fragments ):

    # test the sam file to make sure that each of the reads appears
    sam_fp = open( "./%s/mapped_reads.sam" % output_directory )
    num_mapped_reads = sum( 1 for mapped_read in sc.iter_sam_reads(sam_fp) )
    sam_fp.seek(0)

    # nonmapping read ids
    nonmapping_reads_fname = "./%s/mapped_reads.db.nonmapping" % output_directory
    nonmapping_read_ids = []
    with open( nonmapping_reads_fname ) as nonmapping_fp:
        for line in nonmapping_fp:
            read_id = int(line.strip())
            nonmapping_read_ids.append(read_id)

    # Check that the nonmapping reads had index probes overlapping an intron.
    # Since the test reads are perfect samples from the genome, this is the
    # only way a read should be nonmapping.
    for nonmapping_read_id in nonmapping_read_ids:
        if nonmapping_read_id in categorized_fragments['intron_overlap_probe']:
            pass
        else:
            print "ERROR       :  Nonmapping read id did not have index probes overlapping the intron."
            sys.exit(1)

    for mapped_reads in sc.iter_sam_reads(sam_fp): # each read id
        for mapped_read in mapped_reads:           # each mapping for the readid
            check_mapped_read( mapped_read, fragments, categorized_fragments,
                    introns, genome )

def main():
    output_directory = "smo_rnaseq_sim"

    # Parameters
    # In order to test RNASeq, indexed_seq_len*2 < read_len - max_intron_len
    genome_len = 1000
    frag_len = 100
    read_len = 100
    indexed_seq_len = 20
    max_intron_len = 50
    n_reads  = 100

    genome = sc.build_random_genome( [genome_len,], ["1",] )
    introns, transcriptome = splice_introns_from_genome( genome,
            max_intron_len=max_intron_len )

    # sample fragments from the transcriptome
    fragments = sc.sample_uniformly_from_genome(
            transcriptome, nsamples=n_reads, frag_len=frag_len )
    reads = sc.build_reads_from_fragments(
            transcriptome, fragments,
            read_len=read_len,
            rev_comp=False,
            paired_end=False )

    categorized_fragments = categorize_fragments( genome,
            transcriptome, introns, fragments, read_len, indexed_seq_len )

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
                         indexed_seq_len=indexed_seq_len,
                         assay='r',
                         genome_fnames=genome_fnames )

    check_statmap_output( output_directory,
            genome, transcriptome, introns, fragments, categorized_fragments )

if __name__ == "__main__": main()
