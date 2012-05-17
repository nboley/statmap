import psycopg2
import gzip
import subprocess
import sys
import os
from time import sleep

def get_filenames_by_sample():
    conn = psycopg2.connect("host=eel dbname=labtrack user=nboley")
    cur = conn.cursor()
    query = """
    SELECT sample_name, originating_lab_id, array_agg(file) as files 
    FROM ( 
             SELECT sample_name, originating_lab_id, unnest(file) as file 
               FROM solexa_lane_view 
              WHERE sample_name ~ 'CAGE' 
                AND sample_name ~ 'Dpse'
    ) as foo 
    WHERE file IS NOT NULL 
    GROUP BY sample_name, originating_lab_id;
    """
    
    cur.execute( query )
    # remove spaces from the sample names, and return  the results
    #return [ ( sn.replace(" ", "_"), lab_id.split(";")[0], fns )
    return [ ( sn.replace(" ", "_"), lab_id, fns )
             for sn, lab_id, fns in cur.fetchall() ]

def write_samples_file( samples ):
    '''
    Given a list of tuples of (sn, lab_id, fns) from get_filenames_by_sample,
    generate a tab-delimited file of samples to process

    Each line (tab-delimited):
    sample_type    bio_sample_id    comma sep filenames
    '''
    sys.stdout.write( '\n'.join(
        [
            '\t'.join(( sn, lab_id, ','.join( fns ) ))
            for sn, lab_id, fns in samples
        ]
    ))

def main():
    # Writes to stdout - redirect to save to a file
    write_samples_file( get_filenames_by_sample() )

if __name__ == "__main__": main()
