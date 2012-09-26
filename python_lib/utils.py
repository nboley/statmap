#!/usr/bin/python
'''
Python interface functions for Statmap utilities
'''

import sys
import os

from config_parsing import *
from genome import *
from rawread import *
from mapped_read import *
from trace import *

### FILE* <-> Python File Object ###
# http://svn.python.org/projects/ctypes/trunk/ctypeslib/ctypeslib/contrib/pythonhdr.py
import ctypes
try:
    class FILE(ctypes.Structure):
        pass
    FILE_ptr = ctypes.POINTER(FILE)

    PyFile_FromFile = ctypes.pythonapi.PyFile_FromFile
    PyFile_FromFile.restype = ctypes.py_object
    PyFile_FromFile.argtypes = [FILE_ptr,
                                ctypes.c_char_p,
                                ctypes.c_char_p,
                                ctypes.CFUNCTYPE(ctypes.c_int, FILE_ptr)]

    PyFile_AsFile = ctypes.pythonapi.PyFile_AsFile
    PyFile_AsFile.restype = FILE_ptr
    PyFile_AsFile.argtypes = [ctypes.py_object]
except AttributeError:
    del FILE_ptr

### StatmapOutput ###

CONFIG_FNAME                = "config.dat"
GENOME_FNAME                = "genome.bin"
MAPPED_READS_DB_FNAME       = "mapped_reads.db"
MAPPED_NC_READS_DB_FNAME    = "mapped_NC_reads.db"
FL_DIST_FNAME               = "estimated_fl_dist.txt"

RAWREADS_FNAMES             = [
                                "reads.unpaired",
                                "reads.pair1",
                                "reads.pair2",
                              ]

NC_RAWREADS_FNAMES          = [
                                "reads.NC.unpaired",
                                "reads.NC.pair1",
                                "reads.NC.pair2",
                              ]

class StatmapOutput:
    '''
    Loads all configuration information from a Statmap output directory so it
    can be post-processed by utilities

    Loads data structures into memory lazily. Defaults to only loading the
    genome - use load_ flags to load additional data structures.
    '''
    def _set_num_threads(self, num_threads):
        if num_threads < 1: # num_threads not set by caller, uses default 0
            try:
                # try to determine number of available threads from os
                import multiprocessing
                num_threads = multiprocessing.cpu_count()
                # never set the number of threads to more than 8
                if num_threads > 8:
                    num_threads = 8
            except NotImplementedError:
                # if we can't determine the number of threads, set it to 1
                num_threads = 1
            
        # Set the number of threads used by the shared library code - allows us
        # to specify the number of threads for C functions called by a utility
        # to use
        statmap_o.set_num_threads( num_threads )

    def _set_trace_update_fn(self):
        # set trace update function to use ( depends on assay type )
        if self.config.contents.assay_type == CAGE:
            self.trace_update_fn = \
                statmap_o.update_CAGE_trace_expectation_from_location
        elif self.config.contents.assay_type == CHIP_SEQ:
            self.trace_update_fn = \
                statmap_o.update_chipseq_trace_expectation_from_location

    def _load_mapped_reads( self, load_nc ):
        # load default mapped reads db
        self.mpd_rdb = open_mapped_reads_db_for_reading(
                MAPPED_READS_DB_FNAME )

        # if the mapped reads db is empty, we can't do anything with it
        if self.mpd_rdb.contents.num_mapped_reads == 0:
            print "Mapped reads db is empty (num_mapped_reads=%i)." \
                    % ( self.mpd_rdb.contents.num_mapped_reads )
            sys.exit(-1)

        # load fl dist (if needed - depends on assay type)
        if self.config.contents.assay_type == CHIP_SEQ:
            build_fl_dist_from_filename( self.mpd_rdb, FL_DIST_FNAME )
            build_chipseq_bs_density( self.mpd_rdb.contents.fl_dist )

        # load cond probs db
        self.cond_prbs_db = init_cond_prbs_db_from_mpd_rdb( self.mpd_rdb )

        if load_nc:
            # load negative control reads
            self.NC_mpd_rdb = open_mapped_reads_db_for_reading(
                    MAPPED_NC_READS_DB_FNAME )

            # if the mapped reads db is empty, we can't do anything with it
            if self.NC_mpd_rdb.contents.num_mapped_reads == 0:
                print "Mapped reads db is empty (num_mapped_reads=%i)." \
                        % ( self.NC_mpd_rdb.contents.num_mapped_reads )
                sys.exit(-1)

            # load the fragment length distribution estimate
            build_fl_dist_from_filename( self.NC_mpd_rdb, FL_DIST_FNAME )
            build_chipseq_bs_density( self.NC_mpd_rdb.contents.fl_dist )

            self.NC_cond_prbs_db = init_cond_prbs_db_from_mpd_rdb( self.NC_mpd_rdb )

    def _load_rawreads(self, load_nc ):
        # Note: if loading rawread db does not work, None is returned (NULL ptr)
        self.rawread_db = populate_rawread_db( *RAWREADS_FNAMES )
        if load_nc:
            self.NC_rawread_db = populate_rawread_db( *NC_RAWREADS_FNAMES )

    def __init__(   self,
                    output_directory,
                    num_threads=0,
                    load_mapped_reads=False,
                    load_raw_reads=False,
                    load_nc=False, ):

        # change working directory to Statmap output directory
        try:
            os.chdir( output_directory )
        except IOError as e:
            print "Could not change to directory : %s" % output_directory; raise

        # load statmap run's configuration
        self.config = load_config_from_file( CONFIG_FNAME )

        # These values are saved in the configuration and should not be changed
        statmap_o.set_min_num_hq_bps(
            self.config.contents.min_num_hq_bps )
        statmap_o.set_max_reference_insert_len(
            self.config.contents.max_reference_insert_len )

        self._set_num_threads( num_threads )
        
        # load genome
        self.genome = load_genome_from_disk( GENOME_FNAME )

        # load mapped reads
        if load_mapped_reads:
            self._load_mapped_reads( load_nc )

        # load the rawread db (if requested)
        if load_raw_reads:
            self._load_rawreads( load_nc )

def test():
    """test code in this file"""
    if len(sys.argv) != 2:
        print "Usage: ./utils.py statmap_output_directory/"; sys.exit(1)

    smo = StatmapOutput( sys.argv[1],
                         load_mapped_reads=True,
                         load_raw_reads=True, )

    # see if everything loaded properly

    # did config load properly?
    print "==== CONFIG ===="
    print "genome_fname:", smo.config.contents.genome_fname
    print "min_match_penalty:", smo.config.contents.min_match_penalty
    print "max_penalty_spread:", smo.config.contents.max_penalty_spread

    # did global variables get set properly?
    print "==== GLOBALS ===="
    print "num_threads:", statmap_o.get_num_threads()
    print "min_num_hq_bps:", statmap_o.get_min_num_hq_bps()
    print "max_reference_insert_len:", statmap_o.get_max_reference_insert_len()

    # did genome load properly?
    print "==== GENOME ===="
    print "num_chrs:", smo.genome.contents.num_chrs

    # did mapped reads db load properly?
    print "==== MAPPED_READS_DB ===="
    print "mmapped_data_size:", smo.mpd_rdb.contents.mmapped_data_size
    print "num_mmapped_reads:", smo.mpd_rdb.contents.num_mmapped_reads

    # if necessary, did the frequency length distribution load correctly?
    if smo.config.contents.assay_type == CHIP_SEQ:
        print "==== FL_DIST ===="
        print "min_fl:", smo.mpd_rdb.contents.fl_dist.contents.min_fl
        print "max_fl:", smo.mpd_rdb.contents.fl_dist.contents.max_fl

    # did cond_prbs_db load properly?
    print "==== COND_PRBS_DB ===="
    print "max_rd_id:", smo.cond_prbs_db.contents.max_rd_id

    # did rawread db load properly?
    print "==== RAWREAD_DB ===="
    print "readkey:", smo.rawread_db.contents.readkey
    print "file_type:", smo.rawread_db.contents.file_type

if __name__ == "__main__": test()
