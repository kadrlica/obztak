#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import argparse
import pylab as plt
import pandas as pd
import numpy as np

import fitsio
import skymap
from skymap.utils import setdefaults

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('filename')
parser.add_argument('-o','--outfile')
parser.add_argument('-p','--program',default=None)
parser.add_argument('--propid',default=None)
parser.add_argument('--priority',nargs='+',action='append',default=None,type=int)
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

if args.program:
    fields = fields[fields['program'] == args.program]
if args.propid:
    fields = fields[fields['propid'] == args.propid]
if args.priority:
    fields = fields[np.in1d(fields['priority'],args.priority)]

#fields = fields[fields['expnum'] <= 283862] # DES Y1A1
#fields = fields[fields['expnum'] <= 516819] # DES Y3A2

fig,ax = plt.subplots(2,2,figsize=(16,9))
plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,
                    bottom=0.01,top=0.99)

kwargs = dict(edgecolor='none', alpha=0.2, vmin=-1, vmax=2, s=12)

bands = ['g','r','i','z']
for i,b in enumerate(bands):
    plt.sca(ax.flat[i])

    f = fields[fields['filter'] == b]

    smap = skymap.SurveyMcBryde()
    smap.draw_des()
    smap.draw_milky_way(10)
    smap.draw_fields(f,**kwargs)
    plt.gca().set_title('%s-band'%b)

plt.suptitle('DECam Fields',fontsize=16)

if not args.outfile:
    basename = os.path.basename(args.filename)
    args.outfile = basename.split('.')[0]+'.png'

print("Writing %s..."%args.outfile)
plt.savefig(args.outfile)
plt.ion()
