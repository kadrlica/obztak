#!/usr/bin/env python
"""
Create movie of fields
"""
__author__ = "Alex Drlica-Wagner"
import os
import shutil
import glob
import subprocess
import tempfile

import pylab as plt
import pandas as pd
import numpy as np
import numpy as np

import fitsio
import skymap
from skymap.utils import setdefaults

def hollywood(infiles,outfile=None,delay=10):
    print("Lights, Camera, Action...")
    infiles = np.atleast_1d(infiles)
    if not len(infiles):
        msg = "No input files found"
        raise ValueError(msg)

    infiles = ' '.join(infiles)
    if not outfile: outfile = infiles[0].replace('.png','.gif')
    #cmd='convert -delay %i -quality 100 %s %s'%(delay,infiles,outfile)
    cmd='convert -delay %i -loop 1 -quality 100 %s %s'%(delay,infiles,outfile)
    #print(cmd)
    subprocess.check_call(cmd,shell=True)

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('filename')
parser.add_argument('-o','--outfile')
parser.add_argument('-k','--chunk',default=20,type=int,
                    help='number of nights per frame')
parser.add_argument('-d','--delay',default=10,type=int,
                    help='delay between frames')
args = parser.parse_args()

print("Reading %s..."%args.filename)
if args.filename.endswith(('.fits','.fits.gz','.fz')):
    fields = fitsio.read(args.filename)
else:
    fields = pd.read_csv(args.filename,comment='#').to_records(index=False)

names = fields.dtype.names
names = [str(n.lower()) for n in names]
names = ['ra' if n=='tradeg' else n for n in names]
names = ['dec' if n=='tdecdeg' else n for n in names]
names = ['filter' if n=='band' else n for n in names]
fields.dtype.names = names

# Datetime selection
sel  = ~(fields['datetime'] == 'None')
sel &= ~(np.char.startswith(fields['datetime'],'19')) # Filter bad times
# DES Selections
if True:
    sel &= fields['propid' ] == '2012B-0001'
    sel &= fields['exptime' ] > 30
    sel &= ((fields['ra' ] > 295) | (fields['ra'] < 97))
    sel &= (fields['dec' ] < 10)
fields = fields[sel]

datetime = np.char.replace(fields['datetime'],' ','T').astype(np.datetime64)
timedelta = np.timedelta64(args.chunk, 'D')
tmin = np.datetime64('2012-10-01')
#tmin = np.datetime64('2013-08-01')
#tmin = np.datetime64('2018-08-01')
tmax = tmin + timedelta

outdir = tempfile.mkdtemp()
shutil.rmtree(outdir,ignore_errors=True)
os.makedirs(outdir)

fig,ax = plt.subplots(2,2,figsize=(16,9))
plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,
                    bottom=0.01,top=0.99)
kwargs = dict(edgecolor='none', alpha=0.3, vmin=-1, vmax=2, s=12)
bands = ['g','r','i','z']

fig = plt.figure(figsize=(24,9))
gs=plt.GridSpec(2, 6)
gs.update(left=0.02,right=0.98,bottom=0.02,top=0.97,wspace=0.1,hspace=0.1)
ax = []
ax.append(fig.add_subplot(gs[0,0:2]))
ax.append(fig.add_subplot(gs[0,2:4]))
ax.append(fig.add_subplot(gs[0,4:6]))
ax.append(fig.add_subplot(gs[1,1:3]))
ax.append(fig.add_subplot(gs[1,3:5]))
ax = np.array(ax)
kwargs = dict(edgecolor='none', alpha=0.1, vmin=-1, vmax=2, s=12)
bands = ['g','r','i','z','Y']


filenames = []
while True:

    sel = (datetime >= tmin) & (datetime <= tmax)
    for i,b in enumerate(bands):

        plt.sca(ax.flat[i])
        plt.cla()

        f = fields[sel & (fields['filter']==b)]
        smap = skymap.SurveyMcBryde()
        smap.draw_des(lw=1,color='k')
        smap.draw_milky_way(10,color='k')
        smap.draw_fields(f,**kwargs)
        plt.gca().set_title('%s-band'%b)

    plt.suptitle('%s'%tmax,fontsize=16)

    filename = os.path.join(outdir,'decam_%s.png'%tmax)
    print("Writing %s..."%filename)
    plt.savefig(filename,dpi=100)
    filenames.append(filename)

    if tmax > datetime.max(): break
    else: tmax += timedelta


if not args.outfile:
    basename = os.path.basename(args.filename)
    args.outfile = basename.split('.')[0]+'.gif'

print("Writing %s..."%args.outfile)
hollywood(filenames,args.outfile,args.delay)
shutil.rmtree(outdir,ignore_errors=True)
