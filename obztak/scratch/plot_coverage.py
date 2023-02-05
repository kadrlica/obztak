#!/usr/bin/env python
"""
Plot healpix map coverage.
"""
__author__ = "Alex Drlica-Wagner"
import fitsio
import pylab as plt
import healpy as hp
import matplotlib.cm

import skymap
from skymap.utils import setdefaults

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-o','--outfile',default=None)
parser.add_argument('-n','--nite',type=int,default=20210815)
parser.add_argument('--vmax',default=0.3*3*90,type=float)
args = parser.parse_args()

basename = 'maps/%s/decam_sum_expmap_%s_n1024.fits.gz'
#kwargs = dict(vmin=0, vmax=3*0.3*90)
#kwargs = dict(vmin=0, vmax=3*90)
#kwargs = dict(vmin=0, vmax=900)
kwargs = dict(vmin=0, vmax=args.vmax)
epsilon = 1.0e-7

fig,ax = plt.subplots(2,2,figsize=(16,8))
plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.015,right=0.99,
                    bottom=0.01,top=0.99)

bands = ['g','r','i','z']
for i,b in enumerate(bands):
    plt.sca(ax.flat[i])
    filename = basename%(args.nite,b)
    print("Reading %s..."%filename)
    hpxmap = hp.read_map(filename)

    smap = skymap.SurveyMcBryde(parallels=False,meridians=False)
    para = smap.draw_parallels(fontsize=12)
    para[0][1][0].set_text('')
    meri = smap.draw_meridians(fontsize=12)

    cmap = matplotlib.cm.get_cmap(smap.CMAPS[b])
    cmap.set_under('w')
    kwargs['cmap'] = cmap
    smap.draw_hpxmap(hpxmap - epsilon,**kwargs)
    smap.draw_des(color='0.7')
    smap.draw_milky_way(10)
    plt.gca().set_title('%s band'%b,fontsize=14)

plt.suptitle('DECam Coverage (%s)'%args.nite,fontsize=18)

if not args.outfile:
    args.outfile = 'decam_coverage_%s.png'%args.nite

plt.savefig(args.outfile)
