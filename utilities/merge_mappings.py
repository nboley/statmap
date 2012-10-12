#!/usr/bin/env python

import sys
import os
import re

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from config_parsing import *
from mapped_read import *
from utils import *

def check_config_match( smo_dirs ):
    """ Make sure the input Statmap output directories were built from similar
    data - enough so that merging their mapped reads db's makes sense """

    configs = [ load_config_from_file( os.path.join( smo_dir, CONFIG_FNAME ) )
                for smo_dir in smo_dirs ]

    if len(configs) != len(smo_dirs):
        print "FATAL : Could not load config.dat from %i directories" \
                % ( len(smo_dirs) - len(configs) )
        sys.exit(1)

    ### Check the genome fnames

    # The reads must have been mapped to the same genome for it to make sense
    # for us to merge the mapped reads db's. Checking the genome filename is
    # not perfect (might have multiple genomes with the same filename, or
    # multiple copies of the same genome with different filenames), but this is
    # good enough for now.

    genome_fnames = [ config.contents.genome_fname for config in configs ]

    # Compare every genome_fname to the first
    gfn_match = [ gfn == genome_fnames[0] for gfn in genome_fnames ]
    if not all(gfn_match):
        print "FATAL : One or more genome filenames in the given Statmap output directories do not match: %s" \
                % genome_fnames
        sys.exit(1)

    return

def merge_mapped_reads_dbs( smo_dirs ):

    input_rdbs = [
            open_mapped_reads_db_for_reading(
                os.path.join( smo_dir, MAPPED_READS_DB_FNAME ))
            for smo_dir in smo_dirs ]

    # open the output db for the merged mapped reads
    output_rdb_fname = "mapped_reads.merged.db"
    output_rdb = open_mapped_reads_db_for_writing( output_rdb_fname )

    # pointer to mapped_read_t (pseudo structure, is just a pointer into
    # a block of mmapped memory)
    rd = c_void_p()
    total_mapped_reads = 0

    # iterate over the reads in each input rdb, copying to the merged rdb
    for input_rdb in input_rdbs:
        while get_next_read_from_mapped_reads_db( input_rdb, byref(rd) ) != EOF:
            add_read_to_mapped_reads_db( output_rdb, rd )
            total_mapped_reads += 1

    print "Merged %i mapped reads into %s" \
            % (total_mapped_reads, output_rdb_fname )

    # close the input mapped reads dbs
    for input_rdb in input_rdbs:
        close_mapped_reads_db( input_rdb )

    # close the output rdb to go back and write out the final size
    close_mapped_reads_db( output_rdb )

    # re-open the merged mapped reads db and check it
    output_rdb = open_mapped_reads_db_for_reading( output_rdb_fname )

    if output_rdb.contents.num_mapped_reads != total_mapped_reads:
        print "ERROR : Read %i mapped reads from input dbs, output merged db has %i reads" \
                % ( total_mapped_reads, output_rdb.contents.num_mapped_reads )

    close_mapped_reads_db( output_rdb )

def main():
    if len(sys.argv) < 3:
        print "Usage: ./merge_mappings.py statmap_output_dir_1 statmap_output_dir_2 [statmap_output_dir_3 ...]"
        sys.exit(1)

    smo_dirs = sys.argv[1:]

    check_config_match( smo_dirs )
    merge_mapped_reads_dbs( smo_dirs )

if __name__ == "__main__":
    main()
