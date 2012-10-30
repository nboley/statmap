#!/usr/bin/env python
"""
Produce CSV with info about each read mapped by statmap

Assumes that reads are logged by each thread like this:

read_id %f
...info (format "key value")...
end read_id %f
"""

import sys
import re
from types import StringType
import csv

from collections import namedtuple, defaultdict, OrderedDict

# for future reference - (basic) log line regex
_log_line_re = r'\[(?P<date>[\d-]+) (?P<time>[\d:]+)\] \[(?P<tid>[\d]+)\] (?P<level>[A-Z]+)\s+: (?P<msg>.+)'

def parse_value(val):
    """
    Given a string, try to parse it as an int, a float, or a string
    """
    assert type(val) is StringType # make sure we got passed a string to parse
    try:
        return int(val)
    except ValueError:
        try:
            return float(val)
        except ValueError:
            return val

def parse_log(log_fname):
    """
    Returns a dictionary of logged information about each read processed.
    Handles output from single or multithreaded Statmap runs
    """
    # Sort lines into bins based on their thread id
    threads = defaultdict(list)
    with open(log_fname) as log_fp:
        # Skip log output from bootstrap
        for line in log_fp:
            if "Finding candidate mappings" in line:
                break

        #import pdb; pdb.set_trace();

        for line in log_fp:
            line_data = re.match(_log_line_re, line)
            if line_data:
                threads[int(line_data.group('tid'))].append(line)
            else:
                print "Error parsing regex for line: %s" % line
                continue


    #import pdb; pdb.set_trace()
    reads = []
    for thread_id, thread_log in threads.iteritems():
        log_lines_iterator = iter(thread_log)
        for line in log_lines_iterator:
            if "read_id" in line:
                # get all the log lines about this read
                read_lines = []
                while "end read_id" not in line:
                    read_lines.append(line)
                    line = log_lines_iterator.next()

                # build read dictionary
                read_dict = OrderedDict()
                for read_line in read_lines:
                    line_data = re.match(_log_line_re, read_line)
                    log_msg = line_data.group("msg")
                    k, v = log_msg.split()
                    read_dict[k] = parse_value(v)

                reads.append(read_dict)
            else:
                continue

    return reads

def main():
    if len(sys.argv) != 2:
        print "Usage: ./plot_statmap_log.py statmap.log"
        sys.exit(1)

    log_fname = sys.argv[1]
    # build dictionary of logged data for each read ID
    reads = parse_log(log_fname)

    # sort by read_id
    reads.sort(key=lambda read: read['read_id'])

    #import pdb; pdb.set_trace()

    # write out to csv
    csv_fname = "statmap.csv"
    csv_fp = open(csv_fname, "w")
    try:
        # TODO in a particular order?
        # TODO use OrderedDict so it has the same order as the log file output?
        fieldnames = reads[0].keys()
        writer = csv.DictWriter(csv_fp, fieldnames=fieldnames)
        headers = dict( (n,n) for n in fieldnames )
        writer.writerow(headers)
        for read in reads:
            writer.writerow(read)
    finally:
        csv_fp.close()

if __name__ == "__main__":
    main()
