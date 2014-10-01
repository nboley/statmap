#!/usr/bin/python

### NOTE: Assumptions ###
# build_index has to parse arbitrary input for processing, and infers a lot.
# These are the assumptions it makes:
# 1. FASTA filenames have a .fa extension; diploid map files have .map
# 2. Each FASTA contains a single chromosome
#    XXX - add_chrs_from_fasta_file will process fasta's containing multiple chrs
# 3. Diploid chrs contain "paternal" or "maternal" in their identifier string

import sys
import os
import re

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )

from config_parsing import *
from diploid_map_data import *
from enums import *
from genome import *
from log import *
from mapped_read import *
from trace import *

def usage():
    print "Usage: ./build_index.py indexed_seq_len output_filename genome.fa(s) [diploid.map(s)]"
    sys.exit(1)

def is_diploid_dataset( files ):
    for f in files:
        if f.endswith(".map"):
            return True
    return False

def verify_groups( groups ):
    """
    Verify that the input file groups are valid.
    If any group is invalid, print error and exit.
    """
    # make sure each group is either a single .fa file (haploid group)
    # or two .fa's and a .map (diploid group)
    for prefix, group in groups.items():
        if len(group) == 1:
            # check fasta
            # for now, assume any non-.map file is a fasta
            if not group[0].endswith(".map"): pass
            else:
                print "FATAL    : Haploid file group contains non-FASTA file: %s" % group[0]
                sys.exit(1)
        elif len(group) == 3:
            # check 2 .fa, 1 .map
            num_fa, num_map = 0, 0
            for fn in group:
                # for now, assume any non-.map file is a fasta
                if not fn.endswith(".map"):
                    num_fa += 1
                elif fn.endswith(".map"):
                    num_map += 1
                else:
                    print "FATAL    : Diploid file group contains non-FASTA and non-MAP file: %s" % fn

            if num_fa == 2 and num_map == 1: pass
            else:
                print "FATAL    : Diploid file group is invalid; should contain 2 FASTA and 1 MAP file"
                print "Group    : ", group
                sys.exit(1)
        else:
            print "FATAL    : File group contains %i files" % len(group)
            print "There should be 1 (for a haploid chromosome) or 3 (for a diploid chromosome)"
            sys.exit(1)

def group_input_files( files ):
    groups = {}
    is_diploid = is_diploid_dataset( files )

    if is_diploid:
        # if diploid, group filenames by their prefix
        for f in files:
            # a prefix is delimited by _ or .
            prefix = re.split("[_\.]", f)[0]
            if prefix in groups.keys():
                groups[prefix].append(f)
            else:
                groups[prefix] = [f]
    else:
        # if haploid, treat each file as an independent chr
        # add each file to its own group, using the full filename as the "prefix"
        for f in files:
            groups[f] = [f]

    # verify groups - make sure input is logical
    verify_groups( groups )

    return groups, is_diploid

def get_chr_name( filename ):
    """
    Get the chromsome name from a FASTA file
    Removes newline and leading ">"
    """
    with open(filename) as fp:
        return fp.readline().strip(">\n")

def get_chr_source( filename ):
    """
    Get chromosome source from a FASTA file
    Assumes the FASTA is part of a diploid dataset, and the identifier string
    contains either "paternal" or "maternal"
    """
    with open(filename) as fp:
        identifier = fp.readline().strip()
        if "paternal" in identifier:
            return PATERNAL
        elif "maternal" in identifier:
            return MATERNAL
        else:
            return UNKNOWN

def is_diploid_group( group ):
    for fname in group:
        if fname.endswith(".map"):
            return True

    return False

def get_diploid_groups( groups ):
    diploid_groups = []
    for group in groups:
        if is_diploid_group( group ):
            diploid_groups.append( group )

    return diploid_groups

def build_diploid_map_data_from_group( group, map_data_t, genome ):
    # loop through files in group to get specific filenames
    for fn in group:
        if fn.endswith(".map"):
            map_fname = fn
        else:
            chr_source = get_chr_source( fn )
            if chr_source == PATERNAL:
                paternal_fname = fn
            elif chr_source == MATERNAL:
                maternal_fname = fn
            else:
                print >> sys.stderr, (
                    "Encounted chr of invalid origin %i "
                    "while parsing diploid genome" % chr_source )
                sys.exit(1)

    # get chr names from fasta files
    paternal_chr_name = get_chr_name( paternal_fname )
    maternal_chr_name = get_chr_name( maternal_fname )

    # get chr indexes from genome lookup
    paternal_chr_index = get_chr_index( genome, paternal_chr_name )
    maternal_chr_index = get_chr_index( genome, maternal_chr_name )

    # build diploid_map_data_t struct
    parse_map_file( map_fname, map_data_t, genome,
                    paternal_chr_index, maternal_chr_index )
    index_diploid_map_data( map_data_t )

def main():
    if len(sys.argv) < 4: usage()
    
    # parse arguments
    indexed_seq_len = int( sys.argv[1] )
    if indexed_seq_len%4 != 0:
        raise ValueError, "The indexed sequence length must be a multiple of 4"

    output_fname = sys.argv[2]
    index_fname = output_fname + ".index"

    # sort input files into groups
    groups, is_diploid = group_input_files( sys.argv[3:] )

    # init genome
    genome = init_genome()
    # add groups to genome
    for prefix, group in groups.items():
        if is_diploid_group( group ):
            for filename in group:
                # Add fasta files
                if not filename.endswith(".map"):
                    chr_source = get_chr_source( filename )
                    add_chrs_from_fasta_file( genome, filename, chr_source )
        else:
            for filename in group:
                add_chrs_from_fasta_file( genome, filename, REFERENCE )

    # get diploid chromosome file groups
    diploid_groups = get_diploid_groups( groups.values() )
    # initialize the index
    init_index( byref(genome.contents.index), indexed_seq_len,
                len(diploid_groups) )

    # TODO: build diploid map data for each diploid chr (file group)
    for i, group in enumerate(diploid_groups):
        build_diploid_map_data_from_group( group,
            byref(genome.contents.index.contents.diploid_maps.contents.maps[i]),
            genome )

    # index the genome
    index_genome( genome )
    
    # write the genome to file
    write_genome_to_disk( genome, output_fname )

    # write the index to disk
    build_ondisk_index( genome.contents.index, index_fname )

if __name__ == "__main__": main()
