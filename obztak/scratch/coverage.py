#!/usr/bin/env python
"""
Generic python script.
"""
import sys 
from collections import OrderedDict as odict

import healpy as hp
import numpy as np
import pylab as plt

from obztak.utils.database import Database
from obztak.utils.ortho import DECamBasemap, DECamMcBride

db = Database()
db.connect()

query ="""
SELECT id as expnum, telra as ra, teldec as dec, expTime, filter, 
COALESCE(qc_teff,'NaN') as teff
FROM exposure where exptime >= 30 and discard = False and delivered = True 
and flavor = 'object' and telra between 0 and 360 and teldec between -90 and 90
and filter in ('u','g','r','i','z','Y') and propid NOT LIKE '%-9999'
ORDER BY id;
"""

data = db.query2recarray(query)

exposures = odict([
        ('u',data[data['filter'] =='u']),
        ('g',data[data['filter'] =='g']),
        ('r',data[data['filter'] =='r']),
        ('i',data[data['filter'] =='i']),
        ('z',data[data['filter'] =='z']),
        ('Y',data[data['filter'] =='Y']),
        ])

for b,exp in exposures.items():
    nan = np.isnan(exp['teff'])
    median = np.median(exp[~nan]['teff'])
    print "Median teff for %s-band: %s"%(b,median)

    # Set the teff to a random value from the distribution
    idx = np.random.randint(len(nan)-nan.sum(),size=nan.sum())
    teff = exp['teff'][np.where(~nan)[0][idx]]
    exp['teff'][nan] = teff
                  
def phi2lon(phi): return np.degrees(phi)
def lon2phi(lon): return np.radians(lon)

def theta2lat(theta): return 90. - np.degrees(theta)
def lat2theta(lat): return np.radians(90. - lat)

def ang2vec(lon, lat):
    theta = lat2theta(lat)
    phi = lon2phi(lon)
    vec = hp.ang2vec(theta, phi)
    return vec

def ang2disc(nside, lon, lat, radius, inclusive=False, fact=4, nest=False):
    """
    Wrap `query_disc` to use lon, lat, and radius in degrees.
    """
    vec = ang2vec(lon,lat)
    return hp.query_disc(nside,vec,radius,inclusive,fact,nest)

NSIDE = 1024 # resolution
DECAM = 1.1  # DECam radius

sum_skymaps = odict([
        ('u',np.zeros(hp.nside2npix(NSIDE))),
        ('g',np.zeros(hp.nside2npix(NSIDE))),
        ('r',np.zeros(hp.nside2npix(NSIDE))),
        ('i',np.zeros(hp.nside2npix(NSIDE))),
        ('z',np.zeros(hp.nside2npix(NSIDE))),
        ('Y',np.zeros(hp.nside2npix(NSIDE))),
        ])

max_skymaps = odict([
        ('u',np.zeros(hp.nside2npix(NSIDE))),
        ('g',np.zeros(hp.nside2npix(NSIDE))),
        ('r',np.zeros(hp.nside2npix(NSIDE))),
        ('i',np.zeros(hp.nside2npix(NSIDE))),
        ('z',np.zeros(hp.nside2npix(NSIDE))),
        ('Y',np.zeros(hp.nside2npix(NSIDE))),
        ])

for (band,_max),(band,_sum) in zip(max_skymaps.items(),sum_skymaps.items()):
    print band
    d2 = exposures[band]
    vec = ang2vec(d2['ra'],d2['dec'])
    rad = np.radians(DECAM)
    for i,(v,d) in enumerate(zip(vec,d2)):
        print '\r%s/%s'%(i+1,len(vec)),
        sys.stdout.flush()
        pix = hp.query_disc(NSIDE,v,rad,inclusive=False,fact=4,nest=False)
        _max[pix] = np.clip(_max[pix],d['teff']*d['exptime'],None)
        _sum[pix] += d['teff']*d['exptime']
    print

outbase = "decam_sum_expmap_%s_n%s"
fig = plt.figure(1)
bmap = DECamMcBride()
bmap.draw_parallels(); bmap.draw_meridians()

for band,sky in sum_skymaps.items():
    outfile = outbase%(band,NSIDE)+'.fits.gz'
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky)

    outfile = outbase%(band,NSIDE)+'_mbt.png'
    print "Writing %s..."%outfile
    bmap.draw_hpxmap(np.log10(sky))
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_car.png'
    print "Writing %s..."%outfile
    hp.cartview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(TEFF*EXPTIME)',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_hist.png'
    plt.hist(sky,bins=np.linspace(1,1e3,50),color=COLORS['band'])
    plt.title('%s-band'%band); plt.xlabel('sum(TEFF * EXPTIME)')
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

fig = plt.figure(1)
bmap = DECamMcBride()
bmap.draw_parallels(); bmap.draw_meridians()

outbase = "decam_max_expmap_%s_n%s"
for band,sky in max_skymaps.items():
    outfile = outbase%(band,NSIDE)+'.fits.gz'
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky)

    outfile = outbase%(band,NSIDE)+'_mbt.png'
    print "Writing %s..."%outfile
    bmap.draw_hpxmap(np.log10(sky),vmin=np.log10(30),vmax=np.log10(300))
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_car.png'
    print "Writing %s..."%outfile
    hp.cartview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(max(EXPTIME))',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_hist.png'
    plt.hist(sky,bins=np.linspace(1,1e3,50),color=COLORS['band'])
    plt.title('%s-band'%band); plt.xlabel('max(TEFF * EXPTIME)')
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

fig = plt.figure(1)
outbase = "decam_sum_90s_%s_n%s"
for band,sky in sum_skymaps.items():
    sky = (sky > 90)
    outfile = outbase%(band,NSIDE)+'_mol.png'
    print "Writing %s..."%outfile
    hp.mollview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(sum(TEFF*EXPTIME))',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_car.png'
    print "Writing %s..."%outfile
    hp.cartview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(TEFF*EXPTIME) > 90s',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'.fits.gz'
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky,dtype=bool)

fig = plt.figure(1)
outbase = "decam_max_60s_%s_n%s"
for band,sky in max_skymaps.items():
    sky = (sky > 60)
    outfile = outbase%(band,NSIDE)+'_mol.png'
    print "Writing %s..."%outfile
    hp.mollview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(max(TEFF*EXPTIME)) > 60s',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'_car.png'
    print "Writing %s..."%outfile
    hp.cartview(np.log10(sky),title='DECam Coverage (%s-band)'%band,
                unit='log10(TEFF*EXPTIME)',fig=1)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    outfile = outbase%(band,NSIDE)+'.fits.gz'
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky,dtype=bool)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
