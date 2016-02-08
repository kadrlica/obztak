"""
Decide which fields to observe and time windows to observe.
"""

import os,sys
import numpy as np
import ephem
import logging
import pylab as plt
import numpy.lib.recfunctions as recfunc

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.ortho
from maglites.field import FieldArray
from maglites.utils.constants import BANDS
from maglites.utils import fileio

############################################################

def prepareObservationWindows(nights, horizon=-14., outfile=None):
    """
    Use -14 deg twilight as default for start and end of observation windows.
    """

    observatory = ephem.Observer()
    observatory.lon = maglites.utils.constants.LON_CTIO
    observatory.lat = maglites.utils.constants.LAT_CTIO
    observatory.elevation = maglites.utils.constants.ELEVATION_CTIO
    observatory.horizon = '%.2f'%(horizon)

    sun = ephem.Sun()

    observation_windows = []
    for date, mode in nights:
        observatory.date = '%s 03:00'%(date)
        observatory.date = observatory.date + 24. * ephem.hour 

        time_setting = ephem.Date(observatory.previous_setting(sun))
        time_rising = ephem.Date(observatory.next_rising(sun))
        time_midpoint = ephem.Date(0.5 * (time_setting + time_rising))

        if mode == 'full':
            observation_windows.append([time_setting, time_rising])
        elif mode == 'first':
            observation_windows.append([time_setting, time_midpoint])
        elif mode == 'second':
            observation_windows.append([time_midpoint, time_rising])

    dtype=[('UTC_START','S20'),('UTC_END','S20')]
    observation_windows = np.rec.fromrecords(observation_windows,dtype=dtype)

    if outfile:
        fileio.rec2csv(outfile,observation_windows)

    return observation_windows

############################################################

def prepareTargetList(infile=None, outfile=None, plot=True):

    DECAM_DITHERS = [[0,0],[8/3.,11/3.],[-8/3.,-11/3.],[8/3.,0]]

    # Import the dither function here...
    def dither(ra,dec):
        return ra,dec

    if infile is None:
        infile = os.path.expandvars('$MAGLITESDIR/maglites/data/smash_fields_alltiles.txt')
    #data = np.recfromtxt('smash_fields_alltiles.txt', names=['RA', 'DEC'])
    data = np.recfromtxt(infile, names=True)
    

    sel = maglites.utils.projector.footprint(data['RA'],data['DEC'])

    # The selection could be done after the dither...
    smash_id = data['ID'][sel]
    ra       = data['RA'][sel]
    dec      = data['DEC'][sel]

    nhexes = sel.sum()
    ntilings = len(DECAM_DITHERS)
    nbands = len(BANDS)
    nfields = nhexes*nbands*ntilings

    logging.info("Number of hexes: %d"%nhexes)
    logging.info("Number of tilings: %d"%ntilings)
    logging.info("Number of filters: %d"%nbands)
    
    fields = FieldArray(nfields)
    fields['ID'] = np.arange(1,nfields+1)
    fields['SMASH_ID'] = np.tile(np.repeat(smash_id,nbands),ntilings)
    fields['PRIORITY'].fill(1)
    fields['TILING'] = np.repeat(np.arange(1,ntilings+1),nhexes*nbands)
    fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)

    for i in range(ntilings):
        idx0 = i*nhexes*nbands
        idx1 = idx0+nhexes*nbands
        ra_dither,dec_dither = dither(ra,dec)
        fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
        fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

    logging.info("Number of target fields: %d"%len(fields))

    if plot:
        fig, basemap = maglites.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

        proj = maglites.utils.ortho.safeProj(basemap,fields['RA'],fields['DEC'])
        basemap.scatter(*proj, color='orange', edgecolor='none', s=50)
        if outfile:
            outfig = os.path.splitext(outfile)[0]+'.png'
            fig.savefig(outfig,bbox_inches='tight')
        if not sys.flags.interactive:
            plt.show(block=True)

    if outfile: fields.write(outfile)

    return fields

def dither_fields(ra,dec,offset=(0.4482,0.5975)):
    """
    ra     : RA of nominal field center
    dec    : Dec of nominal field center
    offset : Units of CCD dimension
    """
    pixel_scale = 0.2626 # arcsec

    # note that North is in the -x direction in FITS coordintes
    nx,ny = 4096,2048 

    DECAM = (7,12) # CCD dimensions
    ccdsize = (nx*pixel_scale,ny*pixel_scale)


############################################################


def main():
    """
    # One season prediction
    nights = [['2016/2/10', 'second'],
              ['2016/2/11', 'second'],
              ['2016/2/12', 'second'],
              ['2016/2/13', 'second'],
              ['2016/2/14', 'second'],
              ['2016/2/15', 'second'],
              ['2016/6/27', 'full'],
              ['2016/6/28', 'full'],
              ['2016/6/29', 'full']]
    """

    # Two seasons prediction
    nights = [['2016/2/10', 'second'],
              ['2016/2/11', 'second'],
              ['2016/2/12', 'second'],
              ['2016/2/13', 'second'],
              ['2016/2/14', 'second'],
              ['2016/2/15', 'second'],
              ['2016/6/27', 'full'],
              ['2016/6/28', 'full'],
              ['2016/6/29', 'full'],
              ['2017/2/18', 'second'],
              ['2017/2/19', 'second'],
              ['2017/2/20', 'second'],
              ['2017/2/21', 'second'],
              ['2017/2/22', 'second'],
              ['2017/2/23', 'second'],
              ['2017/6/27', 'full'],
              ['2017/6/28', 'full'],
              ['2017/6/29', 'full']]

    args = parser().parse_args()

    observation_windows = prepareObservationWindows(nights, outfile=args.windows)

    #data, data2 = prepareTargetList('smash_fields_alltiles.txt', outfile='list.txt')
    prepareTargetList('%s/maglites/data/smash_fields_alltiles.txt'%(os.environ['MAGLITESDIR']), outfile=args.fields,plot=args.plot)
    

def parser():
    import argparse
    description = __doc__
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=formatter)
    parser.add_argument('-p','--plot',action='store_true',
                        help='Plot output.')
    parser.add_argument('-f','--fields',default='target_fields.csv',
                        help='List of all target fields.')
    parser.add_argument('-w','--windows',default='observation_windows.csv',
                        help='List of observation windows.')
    return parser

if __name__ == '__main__':
    main()

############################################################
