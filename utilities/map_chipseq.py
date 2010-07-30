#
#
#
#


if __name__ == '__main__':
    # map the real data fully
    # DONE - just call statmap in the normal way

    # marginally map the negative control 
    # DONE - just call statmap in the normal way

    # for each bootstrap samples:
        # update the NC control mappings from the updated real density
        # TODO - really this just needs to parse the marginal density wiggle 
        # into a track, and then call update_chip_seq mappings from trace, and 
        # update trace from chipseq mappings

        # now, we have two mappings - the NC and the real
        # compare them basepair by basepair. If they are different, call the BP a peak
        # finally, compare the positive strand and negative strand on the real and, 
        #   if they are the same, call it a true peak
    print "Hello World."
