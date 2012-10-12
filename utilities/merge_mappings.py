#!/usr/bin/env python

import sys
import os
import shutil
import fnmatch
from collections import defaultdict

# add python_lib to sys.path
sys.path.insert(0, os.path.normpath( sys.path[0] + "/../python_lib") )
from config_parsing import *
from mapped_read import *
from utils import *

def check_configs( smo_dirs ):
    """
    Make sure the Statmap output directories to merge were built from
    similar data - enough so that merging their mapped reads db's makes sense
    """

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

def copy_genome_symlinks( smo_dirs ):
    # Since check_configs made sure all of the mappings were built from the
    # same genome file, we can just copy the symlink from the first Statmap
    # output directory's config
    assert len(smo_dirs[0]) > 0 
    config = load_config_from_file( os.path.join( smo_dirs[0], CONFIG_FNAME ) )
    genome_abspath = config.contents.genome_fname

    os.symlink( genome_abspath, GENOME_FNAME )
    os.symlink( genome_abspath + ".index", GENOME_INDEX_FNAME )
    os.symlink( genome_abspath + ".index.pslocs",
            GENOME_INDEX_PSLOCS_FNAME )
    os.symlink( genome_abspath + ".index.dmap",
            GENOME_INDEX_DIPLOID_MAP_FNAME )

def make_iterative_mapping_directories():
    # TODO - it might be cleaner to put this code in the C in a function that
    # is always called to do iterative mapping. Otherwise we have to maintain
    # two codepaths that do the same thing.

    if SAVE_STARTING_SAMPLES:
        os.mkdir( STARTING_SAMPLES_PATH )

    if SAVE_SAMPLES:
        os.mkdir( RELAXED_SAMPLES_PATH )

    if CALL_PEAKS:
        os.mkdir( CALLED_PEAKS_OUTPUT_DIRECTORY )

    if SAVE_BOOTSTRAP_SAMPLES:
        os.mkdir( BOOTSTRAP_SAMPLES_PATH )
        os.mkdir( BOOTSTRAP_SAMPLES_MAX_PATH )
        os.mkdir( BOOTSTRAP_SAMPLES_MIN_PATH )
        os.mkdir( BOOTSTRAP_SAMPLES_ALL_PATH )

def merge_reads( smo_dirs ):
    # There are several possible combinations of reads files in a Statmap
    # output directory, depending on whether the input was single or paired
    # end, had a negative control, etc.

    # Assuming each reads file begins with read.*, we build a list of read
    # files to cat together based on filename. This will be a dictionary of the
    # read filename to a list of paths to this read filename in each of the
    # Statmap output directories

    read_fname_dict = defaultdict(list)

    for smo_dir in smo_dirs:
        read_fnames = fnmatch.filter(os.listdir(smo_dir), "reads.*")
        for read_fname in read_fnames:
            read_fname_dict[read_fname].append(
                    os.path.join(smo_dir, read_fname) )

    # Two ways to actually cat the files. One is to call `cat` (only works on
    # Linux, but might be faster), and the other is to use shutil.copyfileobj

    # `cat` version
    #for read_fname, path_list in read_fname_dict.items():
        #cat_cmd = [ "cat" ] + path_list + [ ">", read_fname ]
        #os.system( " ".join(cat_cmd) )

    # shutil.copyfileobj version
    for read_fname, path_list in read_fname_dict.items():
        with open( read_fname, "w" ) as read_fp:
            for path in path_list:
                shutil.copyfileobj( open( path ), read_fp )

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

def merge_nonmapping_reads( out_fp, in_fname, offset ):
    with open( in_fname ) as in_fp:
        for line in in_fp:
            updated_read_id = int(line.strip()) + offset
            out_fp.write( str(updated_read_id) + "\n" )

def merge_nonmapped_reads_dbs( smo_dirs ):
    # Note that merge_mapped_reads_dbs will have already created these files
    # through the call to open_mapped_reads_db_for_writing. However, the files
    # will be empty and we can just overwrite them.

    # It is not impossible to use the Statmap interface
    # (add_nonmapping_read_to_mapped_reads_db,
    # add_unmappable_read_to_mapped_reads_db) to do this, perhaps integrating
    # with merge_mapped_reads_dbs. However, it is not worth it atm.

    # the nonmapping mapped reads db files are lists of read_ids delimited by
    # newlines. all we have to do is keep track of the total number of reads
    # merged "so far" so we can correctly offset the nonmapping read ids
    reads_merged = 0

    merged_nonmapping_fp = open( NONMAPPING_READS_DB_FNAME, "w" )
    merged_unmappable_fp = open( UNMAPPABLE_READS_DB_FNAME, "w" )

    for smo_dir in smo_dirs:
        merge_nonmapping_reads( merged_nonmapping_fp, os.path.join(smo_dir,
            NONMAPPING_READS_DB_FNAME), reads_merged )
        merge_nonmapping_reads( merged_unmappable_fp, os.path.join(smo_dir,
            UNMAPPABLE_READS_DB_FNAME), reads_merged )

        # find the total number of input reads for this set of mappings - this
        # gives the offset for the next set of read_ids
        full_reads_fp = None
        if os.path.isfile( os.path.join( smo_dir, "reads.unpaired") ):
            # Single end reads
            full_reads_fp = open( os.path.join( smo_dir, "reads.unpaired" ) )
        elif os.path.isfile( os.path.join( smo_dir, "reads.pair1" ) ):
            # Paired end (both files are the same length)
            full_reads_fp = open( os.path.join( smo_dir, "reads.pair1" ) )
        else:
            print "FATAL : Could not find input read files for single or paired end"
            raise

        # update the number of reads merged by the total number of input reads
        # for this set of mappings (we assume this is FASTQ, so there are
        # 4 lines for each read)
        reads_merged += sum( 1 for line in full_reads_fp ) / 4

        full_reads_fp.close()

    merged_nonmapping_fp.close()
    merged_unmappable_fp.close()

def merge_mappings( output_dir, smo_dirs ):
    # Use the absolute paths of Statmap output directories. Once we've changed
    # into the merged output directory, the paths to Statmap directories given
    # on the command line will be unresolvable unless they are absolute
    for i, smo_dir in enumerate(smo_dirs):
        smo_dirs[i] = os.path.abspath(smo_dir)

    # Create the output directory and change into it
    try:
        os.mkdir( output_dir )
        os.chdir( output_dir )
    except IOError as e:
        print \
"ERROR : Failed to create or change to given output directory %s" % output_dir
        raise e

    copy_genome_symlinks( smo_dirs )
    make_iterative_mapping_directories()

    merge_reads( smo_dirs )
    merge_mapped_reads_dbs( smo_dirs )
    merge_nonmapped_reads_dbs( smo_dirs )

def main():
    if len(sys.argv) < 4:
        print "Usage: ./merge_mappings.py output_dir statmap_dir_1 statmap_dir_2 [statmap_dir_3 ...]"
        sys.exit(1)

    output_dir = sys.argv[1]
    smo_dirs = sys.argv[2:]

    # Make sure it is valid to join the given set of mappings
    check_configs( smo_dirs )
    merge_mappings( output_dir, smo_dirs )

if __name__ == "__main__":
    main()
