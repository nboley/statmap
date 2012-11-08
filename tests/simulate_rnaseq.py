import sys
import random
from itertools import izip
import re
import gzip

import tests as sc # for genome building and sampling functions

def build_introns( sequence, n_introns, min_intron_len, max_intron_len ):
    # Generate random sets of intron starts until we have one that won't
    # produce overlapping introns
    intron_starts = None

    while(True):
        # Subtract max_intron_len to account for the missing sample range at
        # the end
        intron_starts = sorted( random.sample(
            xrange(len(sequence) - max_intron_len), n_introns ))

        # Since they're sorted, we can easily use a linear search to find the
        # minimum distance between intron starts
        min_distance = len(sequence)
        for i in xrange( len(intron_starts) - 1 ):
            # distance between the start of the current intron and the start of
            # the next intron
            intron_distance = intron_starts[i+1] - intron_starts[i]
            if intron_distance < min_distance:
                min_distance = intron_distance

        if min_distance > max_intron_len:
            # if there's no way we could build an intron that would overlap
            # with another, we're done generating possible intron starts
            break

    introns = []

    for intron_start in intron_starts:
        intron_len = random.randint(min_intron_len, max_intron_len)
        introns.append( ( intron_start, intron_start + intron_len ) )

    return introns

def splice_introns_from_sequence( sequence, n_introns, min_intron_len,
        max_intron_len ):
    introns = build_introns( sequence, n_introns, min_intron_len,
            max_intron_len )

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
    # Since reads are sampled from the transcriptome, their coordinates cannot
    # be compared to the genome without accounting the offset produced by any
    # introns.
    offset = 0
    for intron in introns[fragment[0]]: # sorted by start
        # if the intron comes before the fragment, add its length to the offset
        if intron[0] <= fragment[1]:
            offset += ( intron[1] - intron[0] )

    return fragment[1] + offset, fragment[2] + offset

def fragment_overlaps_intron( frag_start, frag_stop, intron_start ):
    # a read's start cannot be inside an intron, since we sampled reads
    # from the transcriptome. a fragment overlaps the intron if the
    # start is < the intron's start and the end is >= the intron start
    if frag_start < intron_start and frag_stop > intron_start:
        return True
    return False

def fragment_index_probe_overlaps_intron( frag_start, frag_stop, intron_start,
        intron_stop, probe_length ):
    # if the range between the fragment start + indexed_seq_len or
    # the range from fragment end - indexed_seq_len is inside the
    # intron, then the index probe overlaps the intron
    if ( frag_start + probe_length > intron_start or
         frag_stop - probe_length < intron_stop ):
        return True
    return False

def categorize_fragments( genome, transcriptome, introns, fragments,
        read_len, indexed_seq_len ):

    categorized_fragments = {
            "no_overlap": [],
            "overlaps_intron": [],
            "probe_overlaps_intron": []
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
            if fragment_overlaps_intron( frag_start, frag_stop, intron[0] ):
                if fragment_index_probe_overlaps_intron( frag_start, frag_stop,
                        intron[0], intron[1], indexed_seq_len ):
                    categorized_fragments["probe_overlaps_intron"].append( read_id )
                else:
                    categorized_fragments["overlaps_intron"].append( read_id )
            else:
                categorized_fragments["no_overlap"].append( read_id )

    return categorized_fragments

def count_mismatches( str1, str2 ):
    num_mismatches = 0
    for c1, c2 in izip( str1, str2 ):
        if c1 != c2:
            num_mismatches += 1

    return num_mismatches

def build_cigar_string_for_fragment( fragment, introns ):
    """Build a cigar string for the fragment given introns"""
    # for now, assume 1 intron
    assert len(introns) == 1

    f_start, f_stop = find_fragment_pos_in_genome( fragment, introns )
    #print "(%i, %i) -> (%i, %i)" % ( fragment[1], fragment[2], f_start, f_stop )

    cigar_entries = []
    # introns are pre-sorted by start
    for intron in introns[fragment[0]]:
        if fragment_overlaps_intron( f_start, f_stop, intron[0] ):
            cigar_entries.append( ( intron[0] - f_start, 'M' ) )
            cigar_entries.append( ( intron[1] - intron[0], 'N' ) )
            # TODO clean up
            cigar_entries.append( ( (f_stop - f_start) - (intron[0] - f_start), 'M' ) )
            
    if cigar_entries == []:
        # then the fragment did not overlap any introns - just has an M entry
        cigar_entries.append( ( fragment[2] - fragment[1], 'M' ) )

    return ''.join( [ "%i%c" % (t[0], t[1]) for t in cigar_entries ] )

def cigar_match( mapped_read, fragments, introns ):
    """Return True if the cigar string for this mapped read matches the known"""
    read_id = int(mapped_read[0])
    truth = fragments[read_id]

    true_cigar = build_cigar_string_for_fragment( truth, introns )

    # cigar string is 6th entry in SAM line
    if true_cigar == mapped_read[5]:
        return True

    # DEBUG (if match failed)
    #print "Mapped Cigar :", mapped_read[5] # cigar entry in SAM line
    #print "  True Cigar :", true_cigar

    return False

def check_mapped_read_sequence( mapped_read, fragments, categorized_fragments,
        introns, genome, num_mutations=0 ):
    # Extract useful information about the mapped read
    # find the source fragment
    read_id = int(mapped_read[0])
    truth = fragments[read_id]

    flag = int(mapped_read[1])
    rev_comp = False
    if flag & 0x10 > 0:
        rev_comp = True

    contig = mapped_read[2]
    start_pos = int(mapped_read[3])

    # parse the cigar string into a list of entries
    cigar_re = r'([0-9]+)([MNS])'
    cigar_string = mapped_read[5]
    cigar_entries = re.findall( cigar_re, cigar_string )

    read_seq = mapped_read[9]

    ## Check - start location should match the known fragment
    # Make sure the start and contig of the mapped read match the truth

    # The fragment's start is in the transcriptome, so we need to find the
    # correct location for comparison to the mapped read's start position
    frag_start_in_genome, frag_stop_in_genome = find_fragment_pos_in_genome(
            truth, introns )

    if( truth[0] != contig or frag_start_in_genome != start_pos ):
        if read_id in categorized_fragments['probe_overlaps_intron']:
            print "NOTICE      :  Mapped read %i (with probe overlapping intron) does not match truth: found (%s, %i) for (%s, %i)" \
                    % ( read_id, contig, start_pos, truth[0], frag_start_in_genome )
        else:
            print "ERROR       :  Mapped read %i does not match truth: found (%s, %i) for (%s, %i)" \
                    % ( read_id, contig, start_pos, truth[0], frag_start_in_genome )
            print "Introns:", introns
            print "Mapped: ", mapped_read
            print "Truth:  ", truth
            sys.exit(1)

    ## Check - all M sequence in the cigar string should match the genome
    # location they mapped to

    # If the read overlaps the intron, there may be multiple mappings found
    # with sequences that are close, but not exact matches with the original
    # genome. This is allowed.
    max_allowed_mismatches = num_mutations
    if ( read_id in categorized_fragments['overlaps_intron'] or
         read_id in categorized_fragments['probe_overlaps_intron'] ):
        max_allowed_mismatches += 3 # empirically chosen

    genome_pos = start_pos
    read_pos = 0

    for entry in cigar_entries:
        entry_len = int( entry[0] )
        entry_op = entry[1]

        # Note: the recheck algorithm in mapped_read.c needed to reverse
        # complement the genome sequence for comparison, since it was comparing
        # it to the original read's sequence. When we print a reverse
        # complement read out to SAM, we print the reverse complemented
        # sequence, so now we can just compare directly to the genome.

        if entry_op == 'M':
            genome_segment = genome[truth[0]][genome_pos:genome_pos+entry_len]
            read_segment = read_seq[read_pos:read_pos+entry_len]

            num_mismatches = count_mismatches( genome_segment.upper(),
                    read_segment.upper() )

            if( num_mismatches > max_allowed_mismatches ):
                print "ERROR       :  Mapped read %i match segment does not match genome (with %i mismatches)" \
                        % (read_id, num_mismatches)
                print "Truth:  ", truth
                print "Introns:", introns
                print "Mapped Read :", read_segment
                print "     Genome :", genome_segment
                print mapped_read

                if read_id in categorized_fragments['overlaps_intron']:
                    print "OVERLAPS INTRON"
                elif read_id in categorized_fragments['probe_overlaps_intron']:
                    print "PROBE OVERLAPS INTRON"
                else:
                    sys.exit(1)

            # Update pointers into sequence
            genome_pos += entry_len
            read_pos += entry_len

        elif entry_op == 'N':
            genome_pos += entry_len # skip the intron (in the genome sequence)
        elif entry_op == 'S':
            read_pos += entry_len # skip the soft clipped bases (in the read)
        else:
            print "ERROR : Invalid CIGAR entry op '%c'" % entry_op
            sys.exit(1)

    return

def check_statmap_output( output_directory, genome, transcriptome, introns,
        fragments, read_len, indexed_seq_len, num_mutations ):

    # Categorize the genome fragments by how they overlap the intron
    categorized_fragments = categorize_fragments( genome,
            transcriptome, introns, fragments, read_len, indexed_seq_len )

    # count mapping reads
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
        if nonmapping_read_id not in categorized_fragments['probe_overlaps_intron']:
            print "ERROR       :  Nonmapping read %i did not have index probes overlapping the intron." \
                    % nonmapping_read_id
            print "Introns :", introns
            print "Fragment :", fragments[nonmapping_read_id]
            sys.exit(1)

    for mapped_reads in sc.iter_sam_reads(sam_fp): # each read id
        read_id = int(mapped_reads[0][0]) # for error messages

        # Make sure at least one of the mappings for this read matches the
        # cigar string we expect from the known fragment and introns
        cigar_matched = False

        # each mapping for the readid
        for mapped_read in mapped_reads:
            check_mapped_read_sequence( mapped_read, fragments, categorized_fragments,
                    introns, genome, num_mutations )
            if cigar_match( mapped_read, fragments, introns ):
                cigar_matched = True

        if not cigar_matched:
            if read_id not in categorized_fragments['probe_overlaps_intron']:
                print "ERROR       :  Did not find a matching cigar string for read id %i (expected %s)" \
                        % (read_id, build_cigar_string_for_fragment(
                            fragments[read_id], introns ) )
                print "Introns :", introns
                print "  Truth :", fragments[read_id]
                print " Mapped :", mapped_read
                sys.exit(1)
            else:
                print "NOTICE      :  Did not find a matching cigar string for read id %i (expected %s)" \
                        % ( read_id, build_cigar_string_for_fragment(
                                fragments[read_id], introns ) )
                print "               Fragment had index probe overlapping an intron, so this was probably unavoidable."

def randomly_mutate_reads( reads, num_mutations ):
    mutated_reads = []

    for seq in reads:
        # pick num_mutations bp's to mutate
        mut_indices = random.sample( xrange(len(seq)), num_mutations )
        mut_seq = list(seq)

        for mut_index in mut_indices:
            orig_bp = mut_seq[mut_index]
            mut_seq[mut_index] = random.choice( list(
                    sc.bps_set.difference( set( [orig_bp, orig_bp.swapcase()] )
                )))

        mut_seq = ''.join(mut_seq)
        error_str = 'h'*len(mut_seq) # FOR NOW - doesn't matter to mismatch search type

        mutated_reads.append( sc.mut_read_t( mut_seq, error_str, seq ) )

    return mutated_reads

def main():
    print "Starting simulate_rnaseq.py ..."

    output_directory = "smo_rnaseq_sim"

    # Parameters
    # In order to test RNASeq, indexed_seq_len*2 < read_len - max_intron_len
    genome_len = 10000
    frag_len = 100
    read_len = 100
    indexed_seq_len = 20
    max_intron_len = 50
    n_reads  = 100000

    genome = sc.build_random_genome( [genome_len,], ["1",] )
    introns, transcriptome = splice_introns_from_genome( genome,
            max_intron_len=max_intron_len )

    # sample fragments from the transcriptome
    fragments = sc.sample_uniformly_from_genome(
            transcriptome, nsamples=n_reads, frag_len=frag_len )
    reads = sc.build_reads_from_fragments(
            transcriptome, fragments,
            read_len=read_len,
            rev_comp=True,
            paired_end=False )

    num_mutations=1
    mutated_reads = randomly_mutate_reads( reads, num_mutations=num_mutations )

    # Write out the test files
    genome_fnames = ( "tmp.genome.fa", )
    with open( genome_fnames[0], "w" ) as genome_fp:
        sc.write_genome_to_fasta( genome, genome_fp, 1 )

    read_fnames = ( "tmp.fastq", )
    with open( read_fnames[0], "w" ) as reads_fp:
        sc.build_single_end_fastq_from_mutated_reads( mutated_reads, reads_fp )

    # Map the data with Statmap
    sc.map_with_statmap( genome_fnames, read_fnames, output_directory,
            indexed_seq_len=indexed_seq_len, assay='r', search_type='e' )

    check_statmap_output( output_directory, genome, transcriptome, introns,
            fragments, read_len, indexed_seq_len, num_mutations )

    print "PASS:    simulate rnaseq"

if __name__ == "__main__": main()
