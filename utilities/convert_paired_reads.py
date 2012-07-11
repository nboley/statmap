#!/usr/bin/python
'''
Converts paired reads from the "old format" to the "new format"

In the old format, paired end reads were identified by '/1' or '/2' suffixes.
It is unclear whether these reads were always in 2 separate files or if they
were sometimes all grouped together into one.

In the new format, paired end reads are always distributed as two files and
reads correspond across the files. These files are specified by the -1 and
-2 flags in Statmap. The read names may or may not have the '/1' and '/2'
suffixes.

So, this conversion is fairly easy. We will use two flags to indicate the names
of the new format output files: -1 and -2, like Statmap. After that, we will
take any number of input filenames on the command line. 
'''

from itertools import islice

def parse_args():
    import argparse; parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument("--output-file-1", "-1",
            dest="of1",
            help="Output filename for first reads"
        )
    parser.add_argument("--output-file-2", "-2",
            dest="of2",
            help="Output filename for second reads"
        )
    parser.add_argument("input_files",
            nargs="+", # one or more optional arguments
            help="Any number of input files in old format FASTQ"
        )

    return parser.parse_args()

def convert_old_fastq( old, new1, new2 ):
    '''
    Convert old fastq into new fastq
    '''
    # open the input file
    with open(old) as infp:
        while True:
            # each read is 4 lines in the FASTQ
            next_read = list(islice(infp, 4))
            if not next_read:
                break

            # process next_read
            if next_read[0].strip().endswith("/1"):
                # if read name ends with /1, append it to new1
                new1.writelines(next_read)
            elif next_read[0].strip().endswith("/2"):
                # if it ends with /2, append it to new2
                new2.writelines(next_read)
            else:
                # else, warning
                print ( "WARNING     :  Expected read name in old style FASTQ "
                        "to end in '/1' or '/2'. Got: " +  next_read[0]
                    )

def main():
    args = parse_args()

    # open file pointers to the output files
    print args.of1
    print args.of2
    ofp1 = open( args.of1, "w" )
    ofp2 = open( args.of2, "w" )

    for input_file in args.input_files:
        convert_old_fastq( input_file, ofp1, ofp2 )

    # close output files
    ofp1.close()
    ofp2.close()

if __name__ == "__main__": main()
