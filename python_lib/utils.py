#!/usr/bin/python
'''
Python interface functions for Statmap utilities
'''

import sys
import os

from config_parsing import *
from genome import *
from mapped_read import *
from trace import *

CONFIG_FNAME            = "config.dat"
GENOME_FNAME            = "genome.bin"
MAPPED_READS_DB_FNAME   = "mapped_reads.db"
FL_DIST_FNAME           = "estimated_fl_dist.txt"

RAWREADS_FNAMES         = [
                            "reads.unpaired",
                            "reads.pair1",
                            "reads.pair2",
                          ]

NC_RAWREADS_FNAMES      = [
                            "reads.NC.unpaired",
                            "reads.NC.pair1",
                            "reads.NC.pair2",
                          ]

class StatmapOutput:
    '''
    Loads all configuration information from a Statmap output directory so it
    can be post-processed by utilities
    '''
    def __init__(   self,
                    output_directory,
                    num_threads=1,
                    load_rawread_db=False, ):

        # change working directory to statmap output directory
        try:
            os.chdir( output_directory )
        except IOError as e:
            print "Could not change to directory : %s" % output_directory; raise

        self.config = load_config_from_file( CONFIG_FNAME )

        # set global variables in shared library
        statmap_o.set_num_threads( num_threads )
        statmap_o.set_min_num_hq_bps(
            self.config.contents.min_num_hq_bps) # stored in config

        # load genome and mapped reads db
        self.genome = load_genome_from_disk( GENOME_FNAME )
        self.mpd_rdb = open_mapped_reads_db( MAPPED_READS_DB_FNAME )

        # load fl dist if necessary
        if self.config.contents.assay_type == CHIP_SEQ:
            build_fl_dist_from_filename( self.mpd_rdb, FL_DIST_FNAME )
            # TODO: make so parameter is fl_dist_t
            build_chipseq_bs_density( self.mpd_rdb )

        # load cond probs db
        self.cond_prbs_db = init_cond_prbs_db_from_mpd_rdb( self.mpd_rdb )

        # prepare mapped reads db
        mmap_mapped_reads_db( self.mpd_rdb )
        index_mapped_reads_db( self.mpd_rdb )

        # load the rawread db (only if we need it)
        self.rawread_db, self.NC_rawread_db = None, None
        # if loading rawread db does not work, will return None (NULL ptr)
        if load_rawread_db:
            self.rawread_db = populate_rawread_db( *RAWREADS_FNAMES )
            self.NC_rawread_db = populate_rawread_db( *NC_RAWREADS_FNAMES )

def test():
    """test code in this file"""
    if len(sys.argv) != 2:
        print "Usage: ./utils.py statmap_output_directory/"; sys.exit(1)

    smo = StatmapOutput( sys.argv[1] )

    # see if everything loaded properly

    # did config loaded properly?
    print "==== CONFIG ===="
    print "genome_fname:", smo.config.contents.genome_fname
    print "min_match_penalty:", smo.config.contents.min_match_penalty
    print "max_penalty_spread:", smo.config.contents.max_penalty_spread

    # did global variables get set properly?
    print "==== GLOBALS ===="
    print "num_threads:", statmap_o.get_num_threads()
    print "min_num_hq_bps:", statmap_o.get_min_num_hq_bps()

    # did genome load properly?
    print "==== GENOME ===="
    print "num_chrs:", smo.genome.contents.num_chrs

    # did mapped reads db load properly?
    print "==== MAPPED_READS_DB ===="
    print "mmapped_data_size:", smo.mpd_rdb.contents.mmapped_data_size
    print "num_mmapped_reads:", smo.mpd_rdb.contents.num_mmapped_reads

    # if necessary, did the frequency length distribution load correctly?
    # TODO: crux of double pointer indirection issue. wtf
    '''
    if smo.config.contents.assay_type == CHIP_SEQ:
        print "==== FL_DIST ===="
        print "min_fl:", smo.mpd_rdb.contents.fl_dist
        print "max_fl:", smo.mpd_rdb.contents.fl_dist
    '''

    # did cond_prbs_db load properly?
    print "==== COND_PRBS_DB ===="
    print "max_rd_id:", smo.cond_prbs_db.contents.max_rd_id

if __name__ == "__main__": test()
