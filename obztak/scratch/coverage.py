#!/usr/bin/env python
"""
Generic python script.
"""
import os,sys
from collections import OrderedDict as odict
import warnings
import datetime

import numpy as np
import pylab as plt
import pandas as pd
import healpy as hp
import fitsio

from obztak.utils.database import Database
from obztak.utils.constants import COLORS
from obztak.utils.ortho import DECamBasemap, DECamMcBride
from obztak.utils import fileio
from skymap.survey import MaglitesSkymap

warnings.simplefilter('ignore', UserWarning)

date = datetime.datetime.now().strftime('%Y%m%d')
outfile = 'decam-exposures-%s.fits.gz'%date

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('outfile',default=outfile)
parser.add_argument('-m','--maps',action='store_true')
parser.add_argument('-p','--plot',action='store_true')
args = parser.parse_args()


NSIDE = 1024 # resolution
DECAM = 1.1  # DECam radius
#BANDS = ['u','g','r','i','z','Y','VR']
BANDS = ['u','g','r','i','z','Y']
#BANDS = ['VR']

SkymapCls,suffix = DECamMcBride,'_mbt'
#SkymapCls,suffix = MaglitesSkymap,'_ort'

# COALESCE(qc_teff,'NaN')
QUERY ="""
SELECT id as expnum, telra as ra, teldec as dec, exptime, filter, propid,
(CASE WHEN qc_teff is NULL THEN 'NaN' WHEN qc_teff=-1 THEN 'NAN' ELSE qc_teff END) as teff,
(CASE WHEN qc_fwhm is NULL THEN 'NaN' WHEN qc_fwhm=-1 THEN 'NAN' ELSE qc_fwhm END) as fwhm
FROM exposure WHERE
aborted=False and exposed=True and digitized=True and built=True and delivered=True and discard=False
and flavor = 'object' and telra between 0 and 360 and teldec between -90 and 90
and exptime >= 30
and filter in (%s) and propid NOT LIKE '%%-9999'
-- and propid = '2014B-0404' -- DECaLS
-- and (date < current_date - interval '18 months' or propid='2019A-0305')
ORDER BY id;
"""%(",".join(["'%s'"%b for b in BANDS]))

shephard = """
SELECT ra,teldec,filter,exptime,propid,
to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE
FROM exposure WHERE
(propid = '2016B-0288' or propid = '2017A-0367')
and flavor = 'object' order by date
"""

print(QUERY)

def merge_qa(data,filename):
    """
    Update qa properties based on other values.
    """
    df = fileio.read_csv(filename)
    print("Loaded %i QA values..."%(len(df)))

    x = pd.DataFrame(data).merge(df,left_on='expnum',right_on='expnum',how='left')
    sel = np.isnan(data['teff']) & ~np.isnan(x['teff_y'])

    print("Updating %i QA values..."%(sel.sum()))
    data['teff'][sel] = x[sel]['teff_y']

    if 'fwhm' in data.dtype.names and 'fwhm_y' in x.columns:
        data['fwhm'][sel] = x[sel]['fwhm_y']


args.db = True
if args.db:
    db = Database()
    db.connect()
    data = db.query2recarray(QUERY)

    args.delve = True
    if args.delve:
        print("Reading DELVE QA values...")
        merge_qa(data,'data/delve-qa-20210127.csv.gz')

    if os.path.exists(args.outfile): os.remove(args.outfile)
    print("Writing %s..."%args.outfile)
    fitsio.write(args.outfile,data)

# Do we want to make maps?
if not args.maps: sys.exit()

data = fitsio.read(args.outfile)

exposures = odict([(b,data[data['filter'] ==b]) for b in BANDS])
sum_skymaps = odict([(b,np.zeros(hp.nside2npix(NSIDE))) for b in BANDS])
max_skymaps = odict([(b,np.zeros(hp.nside2npix(NSIDE))) for b in BANDS])

for band,exp in exposures.items():
    print "%s-band..."%band
    nan = np.isnan(exp['teff'])
    median = np.median(exp[~nan]['teff'])
    print "  Median teff: %.2f"%(median)
    print "  Fraction without teff: %i%%"%(100*nan.sum()/float(len(nan)))
    if not nan.sum(): continue

    ## DEPRECATED: Now using DELVE QA info.
    ## Set the teff to a random value from the distribution
    #idx = np.random.randint(len(nan)-nan.sum(),size=nan.sum())
    #teff = exp['teff'][np.where(~nan)[0][idx]]
    #exp['teff'][nan] = teff

for (band,_max),(band,_sum) in zip(max_skymaps.items(),sum_skymaps.items()):
    print band
    d2 = exposures[band]
    #d2 = d2[ (d2['teff'] >= 0.1) ]
    d2 = d2[ (d2['teff'] >= 0.2) ]
    #d2 = d2[ (d2['teff'] >= 0.3) ]
    d2 = d2[ (d2['teff']*d2['exptime'] > 18) ] # teff*Texp > 18s (0.2 * 90s)
    d2 = d2[ (d2['fwhm'] < 2.0) ]
    vec = hp.ang2vec(d2['ra'],d2['dec'],lonlat=True)
    rad = np.radians(DECAM)
    for i,(v,d) in enumerate(zip(vec,d2)):
        print '\r%s/%s'%(i+1,len(vec)),
        sys.stdout.flush()
        pix = hp.query_disc(NSIDE,v,rad,inclusive=False,fact=4,nest=False)
        _max[pix] = np.clip(_max[pix],d['teff']*d['exptime'],None)
        _sum[pix] += d['teff']*d['exptime']

    hdr = dict(DATE=str(datetime.date.today()))
    outfile = "decam_sum_expmap_%s_n%s.fits.gz"%(band,NSIDE)
    print "Writing %s..."%outfile
    hp.write_map(outfile,_sum,extra_header=hdr)

    outfile = "decam_max_expmap_%s_n%s.fits.gz"%(band,NSIDE)
    print "Writing %s..."%outfile
    hp.write_map(outfile,_max,extra_header=hdr)

    print

# Do we want to plot?
if not args.plot: sys.exit(0)

cbar_kwargs = dict(orientation='horizontal',aspect=40)
fig = plt.figure(1); plt.clf()
outbase = "decam_sum_expmap_%s_n%s"
label = r'$\log_{10} \sum(t_{\rm eff} t_{\rm exp})$'
for band,sky in sum_skymaps.items():
    if not sky.sum():
        print ("No exposures found in %s band; skipping..."%band)
        continue

    #outfile = outbase%(band,NSIDE)+'.fits.gz'
    #print "Writing %s..."%outfile
    #hp.write_map(outfile,sky)

    title = '%s-band'%band
    outfile = outbase%(band,NSIDE)+suffix+'.png'
    print "Writing %s..."%outfile
    bmap = SkymapCls(); bmap.draw_des()
    bmap.draw_hpxmap(np.log10(sky));
    plt.colorbar(label=label,**cbar_kwargs)
    plt.title(title)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    #outfile = outbase%(band,NSIDE)+'_cyl.png'
    #print "Writing %s..."%outfile
    #bmap = DECamBasemap(projection='cyl',celestial=True); bmap.draw_des()
    #bmap.draw_hpxmap(np.log10(sky));
    #plt.colorbar(label=label,**cbar_kwargs)
    #plt.title(title)
    #plt.savefig(outfile,bbox_inches='tight')
    #plt.clf()

    outfile = outbase%(band,NSIDE)+'_hist.png'
    print "Writing %s..."%outfile
    plt.hist(sky,bins=np.linspace(1,1e3,50),color=COLORS[band])
    plt.title('%s-band'%band); plt.xlabel('sum(TEFF * EXPTIME)')
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

fig = plt.figure(1); plt.clf()
outbase = "decam_max_expmap_%s_n%s"
label = r'$\log_{10} (\max(t_{\rm eff} t_{\rm exp})$'
for band,sky in max_skymaps.items():
    if not sky.sum():
        print ("No exposures found in %s band; skipping..."%band)
        continue

    #outfile = outbase%(band,NSIDE)+'.fits.gz'
    #print "Writing %s..."%outfile
    #hp.write_map(outfile,sky)

    outfile = outbase%(band,NSIDE)+suffix+'.png'
    print "Writing %s..."%outfile
    title = '%s-band'%band
    bmap = SkymapCls(); bmap.draw_des()
    bmap.draw_hpxmap(np.log10(sky));
    plt.colorbar(label=label,**cbar_kwargs)
    plt.title(title);
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    #outfile = outbase%(band,NSIDE)+'_cyl.png'
    #print "Writing %s..."%outfile
    #bmap = DECamBasemap(projection='cyl',celestial=True); bmap.draw_des()
    #bmap.draw_hpxmap(np.log10(sky));
    #plt.colorbar(label=label,**cbar_kwargs)
    #plt.title(title)
    #plt.savefig(outfile,bbox_inches='tight')
    #plt.clf()

    outfile = outbase%(band,NSIDE)+'_hist.png'
    print "Writing %s..."%outfile
    plt.hist(sky,bins=np.linspace(1,5e2,50),color=COLORS[band])
    plt.title('%s-band'%band); plt.xlabel('max(TEFF * EXPTIME)')
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

fig = plt.figure(1); plt.clf()
outbase = "decam_sum_90s_%s_n%s"
label = r'$\sum(t_{\rm eff} t_{\rm exp}) > 90$'
for band,sky in sum_skymaps.items():
    if not sky.sum():
        print ("No exposures found in %s band; skipping..."%band)
        continue

    sky = (sky > 90)
    outfile = outbase%(band,NSIDE)+'.fits.gz'
    title = '%s-band'%band
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky,dtype=bool)

    outfile = outbase%(band,NSIDE)+suffix+'.png'
    print "Writing %s..."%outfile
    bmap = SkymapCls(); bmap.draw_des()
    bmap.draw_hpxmap(np.log10(sky));
    plt.title(title)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    #outfile = outbase%(band,NSIDE)+'_cyl.png'
    #print "Writing %s..."%outfile
    #bmap = DECamBasemap(projection='cyl',celestial=True); bmap.draw_des()
    #bmap.draw_hpxmap(np.log10(sky));
    #plt.title(title)
    #plt.savefig(outfile,bbox_inches='tight')
    #plt.clf()


fig = plt.figure(1); plt.clf()
outbase = "decam_max_60s_%s_n%s"
label = r'$\max(t_{\rm eff} t_{\rm exp}) > 60$'
for band,sky in max_skymaps.items():
    if not sky.sum():
        print ("No exposures found in %s band; skipping..."%band)
        continue

    sky = (sky > 60)
    outfile = outbase%(band,NSIDE)+'.fits.gz'
    title = '%s-band'%band
    print "Writing %s..."%outfile
    hp.write_map(outfile,sky,dtype=bool)

    outfile = outbase%(band,NSIDE)+suffix+'.png'
    print "Writing %s..."%outfile
    bmap = SkymapCls(); bmap.draw_des()
    bmap.draw_hpxmap(np.log10(sky));
    plt.title(title)
    plt.savefig(outfile,bbox_inches='tight')
    plt.clf()

    #outfile = outbase%(band,NSIDE)+'_cyl.png'
    #print "Writing %s..."%outfile
    #bmap = DECamBasemap(projection='cyl',celestial=True); bmap.draw_des()
    #bmap.draw_hpxmap(np.log10(sky));
    #plt.title(title)
    #plt.savefig(outfile,bbox_inches='tight')
    #plt.clf()

