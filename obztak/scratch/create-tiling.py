#!/usr/bin/env python
"""
Create various tiling schemes.
"""
__author__ = "Alex Drlica-Wagner"
import fitsio
import numpy as np

from obztak.utils import fileio
from obztak.survey import Survey
from obztak.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS,COLORS,DECAM

from obztak.ddgs import DdgsSurvey
from obztak.maglites import MaglitesSurvey
from obztak.maglites2 import Maglites2Survey

DECALS = fitsio.read(fileio.get_datafile('decam-tiles_obstatus.fits'))
SMASH  = np.genfromtxt(fileio.get_datafile('smash_fields_alltiles.txt'),
                       dtype=[('TILEID','>i4'),('RA','>f8'),('DEC','>f8')])
N = len(DECALS)/3
DTYPE = [
    ('TILEID','>i4'),('PASS','>i2'),('RA','>f8'),('DEC','>f8'),
    ('IN_DES',bool),('IN_DECALS',bool),('IN_SMASH',bool),('IN_MAGLITES',bool),
]

def create_fields(base, offsets, dither):
    base = np.copy(base)
    fields =[]
    for i,off in enumerate(offsets):
        data = np.zeros(len(base),dtype=DTYPE)
        data['TILEID'] = base['TILEID']
        data['PASS'] = i+1
        ra_dither,dec_dither = dither(base['RA'],base['DEC'],off[0],off[1])
        data['RA']  = ra_dither
        data['DEC'] = dec_dither
        fields.append(data)
    return np.concatenate(fields)

def create_smash_fields():
    offsets = [(0.,0.), (1.0,0.0), (-1.0,0.0), (0.0,-0.75)]
    base = SMASH
    dither = Survey.smash_dither
    return create_fields(base, offsets, dither)

def create_bliss_fields():
    offsets = [(0., 0.), (8/3.*CCD_X, -11/3.*CCD_Y),
               (8/3.*CCD_X, 8/3.*CCD_Y), (-8/3.*CCD_X, 0.)]
    base = DECALS[:N]
    dither = Survey.decals_rotate
    return create_fields(base, offsets, dither)

def create_decals_fields():
    """ Replicate first 3 DECaLS tilings from
    decam-tiles_obstatus.fits and add a 4th tiling.
    """
    offsets = [(0.,0.),(-0.2917, 0.0833),
               (-0.5861, 0.1333),
               (-0.8805,0.1833)
    ]
    base = DECALS[:N]
    dither = Survey.decals_rotate
    return create_fields(base, offsets, dither)

def find_decals_rotate(n=1):
    """ Find the DECaLS rotation scheme. """
    if n > 3: raise ValueError("Only 3 chunks of DECaLS exposures.")
    from scipy.optimize import root, minimize
    ra0,dec0 = DECALS[:N]['RA'],DECALS[:N]['DEC']
    ra1,dec1 = DECALS[n*N:(n+1)*N]['RA'],DECALS[n*N:(n+1)*N]['DEC']
    def fn(p):
        dx,dy=p
        ra2,dec2 = Survey.decals_rotate(ra0,dec0,dx,dy)
        return np.sum( (ra1-ra2)**2 + (dec1-dec2)**2)

    result = minimize(fn,x0=(-0.3,0.1), bounds=[[-1.0,1.0],[-1.0,1.0]],
                       method='BFGS',tol=1e-5)
    print result['x']
    return result

def set_footprints(data):
    """ Set survey footprint information """
    data['IN_DES'] = Survey.footprintDES(data['RA'],data['DEC'])
    data['IN_DECALS'] = Survey.footprintDECALS(data['RA'],data['DEC'])
    data['IN_SMASH'] = Survey.footprintSMASH(data['RA'],data['DEC'],DECAM)

    data['IN_MAGLITES'] = MaglitesSurvey.footprint(data['RA'],data['DEC'])
    data['IN_MAGLITES'] |= Maglites2Survey.footprint(data['RA'],data['DEC'])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    smash = create_smash_fields()
    set_footprints(smash)
    outfile = 'decam-tiles-smash-v1.fits.gz'
    fitsio.write(outfile,smash,clobber=True)

    decals = create_decals_fields()
    set_footprints(decals)
    outfile = 'decam-tiles-decals-v1.fits.gz'
    fitsio.write(outfile,decals,clobber=True)

    # v0 : had the 0-rotation being executed
    # v1 : should be functionally identical to v0
    bliss = create_bliss_fields()
    set_footprints(bliss)
    outfile = 'decam-tiles-bliss-v2.fits.gz'
    fitsio.write(outfile,bliss,clobber=True)
