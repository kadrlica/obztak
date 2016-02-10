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

from maglites.utils.ortho import datestring
from maglites.field import FieldArray
from maglites.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS
from maglites.utils import fileio

plt.ion()

############################################################

def prepareObservationWindows(nights, horizon=-14., standards=True, outfile=None):
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
            if standards: 
                time_midpoint = datestring(time_midpoint - STANDARDS)
            observation_windows.append([time_setting, time_midpoint])
        elif mode == 'second':
            if standards: 
                time_midpoint = datestring(time_midpoint + STANDARDS)
            observation_windows.append([time_midpoint, time_rising])

    dtype=[('UTC_START','S20'),('UTC_END','S20')]
    observation_windows = np.rec.fromrecords(observation_windows,dtype=dtype)

    if outfile:
        fileio.rec2csv(outfile,observation_windows)

    return observation_windows

############################################################

def prepareTargetList(infile=None, outfile=None, mode='smash_dither', plot=True):
    # Import the dither function here...
    #def dither(ra,dec,dx,dy):
    #    return ra,dec

    if mode is None or mode.lower() == 'none':
        def dither(ra,dec,dx,dy):
            return ra,dec
        TILINGS = [(0,0),(0,0),(0,0),(0,0)]
    elif mode.lower() == 'smash_dither':
        TILINGS = [(0,0), (1.0,0.0), (-1.0,0.0), (0.0,-0.75)]
        dither = smash_dither
    elif mode.lower() == 'smash_rotate':
        TILINGS = [(0,0), (0.75,0.75), (-0.75,0.75), (0.0,-0.75)]
        dither = smash_rotate
    elif mode.lower() == 'decam_dither':
        TILINGS = [(0., 0.),(8/3.*CCD_X, -11/3.*CCD_Y),
                   (8/3.*CCD_X, 8/3.*CCD_Y),(-8/3.*CCD_X, 0.)]
    
    if infile is None:
        infile = os.path.expandvars('$MAGLITESDIR/maglites/data/smash_fields_alltiles.txt')
    data = np.recfromtxt(infile, names=True)
    
    # Apply footprint selection after tiling/dither
    #sel = maglites.utils.projector.footprint(data['RA'],data['DEC'])

    # This is currently a non-op
    smash_id = data['ID']
    ra       = data['RA']
    dec      = data['DEC']

    nhexes = len(data)
    #ntilings = len(DECAM_DITHERS)
    ntilings = len(TILINGS)
    nbands = len(BANDS)
    nfields = nhexes*nbands*ntilings

    logging.info("Number of hexes: %d"%nhexes)
    logging.info("Number of tilings: %d"%ntilings)
    logging.info("Number of filters: %d"%nbands)
    
    fields = FieldArray(nfields)
    fields['HEX'] = np.tile(np.repeat(smash_id,nbands),ntilings)
    fields['PRIORITY'].fill(1)
    fields['TILING'] = np.repeat(np.arange(1,ntilings+1),nhexes*nbands)
    fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)

    #for i in range(ntilings):
    for i,tiling in enumerate(TILINGS):
        idx0 = i*nhexes*nbands
        idx1 = idx0+nhexes*nbands
        ra_dither,dec_dither = dither(ra,dec,tiling[0],tiling[1])
        fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
        fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

    # Apply footprint selection after tiling/dither
    sel = maglites.utils.projector.footprint(fields['RA'],fields['DEC'])
    sel = sel & (fields['DEC'] > maglites.utils.constants.SOUTHERN_REACH)
    fields = fields[sel]

    logging.info("Number of target fields: %d"%len(fields))

    if plot:
        fig, basemap = maglites.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

        proj = maglites.utils.ortho.safeProj(basemap,fields['RA'],fields['DEC'])
        basemap.scatter(*proj, c=fields['TILING'], edgecolor='none', s=50, cmap='Spectral',vmin=0,vmax=len(TILINGS))
        colorbar = plt.colorbar(label='Tiling')
        
        if outfile:
            outfig = os.path.splitext(outfile)[0]+'.png'
            fig.savefig(outfig,bbox_inches='tight')
        if not sys.flags.interactive:
            plt.show(block=True)

    if outfile: fields.write(outfile)

    return fields

############################################################
#
#def dither_fields(ra,dec,offset=(0.4482,0.5975)):
#    """
#    ra     : RA of nominal field center
#    dec    : Dec of nominal field center
#    offset : Units of CCD dimension
#    """
#    pixel_scale = 0.2626 # arcsec
#
#    # note that North is in the -x direction in FITS coordintes
#    nx,ny = 4096,2048 
#
#    DECAM = (7,12) # CCD dimensions
#    ccdsize = (nx*pixel_scale,ny*pixel_scale)
#
############################################################

def smash_dither(ra,dec,dx,dy):
    """
    Convert to SMASH coordinates, then dither, and convert back.
    dx, dy specify offset in decimal degrees.
    """
    ra0,dec0 = SMASH_POLE
    # Rotate the pole at SMASH_POLE to (0,-90)
    R1 = maglites.utils.projector.SphericalRotator(ra0,90+dec0); 
    ra1,dec1 = R1.rotate(ra,dec)
    # Dither the offset 
    ra2 = ra1+dx/np.cos(np.radians(dec1))
    dec2 = dec1+dy
    # Rotate back to the original frame (keeping R2)
    return R1.rotate(ra2,dec2,invert=True)

############################################################

def decam_dither(ra,dec,dx,dy):
    """
    Conventional dither in DECam coordinates.
    dx, dy specify offset in decimal degrees.
    Consider alternative DECam dither in "image" coordinates (see scratch/dither.py).
    """
    out = []
    for _ra,_dec in zip(ra,dec):
        out.append(maglites.utils.projector.SphericalRotator(_ra,_dec).rotate(dx,dy,invert=True))
    return np.array(out).T

############################################################

def smash_rotate(ra,dec,dx,dy):
    """
    Rotate to SMASH coordinates, then rotate by desired offset, and convert back.
    dx, dy specify offset in decimal degrees.
    """
    ra0,dec0 = SMASH_POLE
    # Rotate the pole at SMASH_POLE to (0,-90)
    R1 = maglites.utils.projector.SphericalRotator(ra0,90+dec0); 
    ra1,dec1 = R1.rotate(ra,dec)
    # Rotate around the SMASH pole
    R2 = maglites.utils.projector.SphericalRotator(dx,dy)
    ra2,dec2 = R2.rotate(ra1,dec1)
    # Rotate back to the original frame (keeping the R2 shift)
    return R1.rotate(ra2,dec2,invert=True)

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

    observation_windows = prepareObservationWindows(nights, outfile=args.windows, standards = args.standards)

    #data, data2 = prepareTargetList('smash_fields_alltiles.txt', outfile='list.txt')
    prepareTargetList('%s/maglites/data/smash_fields_alltiles.txt'%(os.environ['MAGLITESDIR']), outfile=args.fields,plot=args.plot)
    
############################################################

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
    parser.add_argument('-d','--dither',default='smash_dither',
                        help='Dithering scheme.')
    parser.add_argument('--no-standards',action='store_false',dest='standards', help = "Don't include time for standard star observations.")
    return parser

############################################################

if __name__ == '__main__':
    main()

############################################################
