"""
Automate the process of creating standard star field observing scripts
for the g and r band only.
"""

import glob
import os.path
import json

import maglites.utils.fileio

############################################################

INPUT_STANDARDS_DIR = '/usr/remote/user/sispi/decam/ExposureScripts/DES/standards/'
INPUT_FILTER_SETS = ['grizY', 'Yzirg']
OUTPUT_STANDARDS_DIR = '/usr/remote/user/DECamObserver/kadrlica/maglites/ExposureScripts/standards_test/'
OUTPUT_FILTER_SETS = ['gr', 'rg']

############################################################

def convert(infile, outfile, output_filter_set):
    """
    Convert a single standard star field
    """
    input_data = maglites.utils.fileio.read_json(infile)

    ra, dec = None, None
    for ii in range(0, len(input_data)):
        if 'RA' in input_data[ii].keys() and 'dec' in input_data[ii].keys():
            ra = input_data[ii]['RA']
            dec = input_data[ii]['dec']
            break

    if ra is None or dec is None:
        print 'Coordinates not properly set'

    output_data = []
    for filter in output_filter_set:
        for ii in range(0, len(input_data)):
            if input_data[ii]['filter'] == filter:
                output_data.append(input_data[ii])
                output_data[-1]['RA'] = ra
                output_data[-1]['dec'] = dec
                break

    if len(output_data) != len(output_filter_set):
        print 'Something terrible happened'

    maglites.utils.fileio.write_json(outfile, output_data)

############################################################

def run():
    for input_filter_set, output_filter_set in zip(INPUT_FILTER_SETS, OUTPUT_FILTER_SETS):
        infiles = glob.glob('%s/%s/*'%(INPUT_STANDARDS_DIR, input_filter_set))
        for infile in infiles:
            outfile = '%s/%s/%s'%(OUTPUT_STANDARDS_DIR, 
                                  output_filter_set, 
                                  os.path.basename(infile).replace(input_filter_set, output_filter_set))
            print 'Infile:', infile
            print 'Outfile:', outfile, '\n'
            convert(infile, outfile, output_filter_set)
            #break
        #break

############################################################

if __name__ == "__main__":
        run()

############################################################

