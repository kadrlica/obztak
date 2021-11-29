#!/usr/bin/env python
"""
Update survey fields
"""
__author__ = "Alex Drlica-Wagner"
import copy

import fitsio
import numpy as np
import pylab as plt

import skymap

from obztak.utils import fileio
import obztak.delve
from obztak.delve import DelveFieldArray
BANDS = ['g','r','i','z']

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('old')
parser.add_argument('new')
parser.add_argument('-o','--outfile',default='update_target_fields.csv')
args = parser.parse_args()

# Old field array
old = DelveFieldArray.load(args.old)
# New field array
new = DelveFieldArray.load(args.new)
# Completed DELVE observations
db  = DelveFieldArray.load_database()
# New completed fields from DELVE observations
delve = np.in1d(new.unique_id,db.unique_id)

print("Comparing to old and new target fields...")
print(args.old)
print(args.new)
if len(old) != len(new): print("Different number of fields")

### New fields (not done before, now done)
#done = (new['PRIORITY'] < 0) & (old['PRIORITY'] >= 0)
done = (new['PRIORITY'] < 0) & np.in1d(new.unique_id, old.unique_id[old['PRIORITY'] >= 0])

print("%i exposures newly considered done"%done.sum())
for b in BANDS:
    print("  %s: %i "%(b,(done & (new['FILTER']==b)).sum()))

plt.figure()
smap = skymap.SurveyMcBryde()
smap.draw_fields(new[done])
smap.draw_des()
plt.title('Exposures newly considered done')
plt.show()

### New fields that were previously considered done and are not now
undone = (new['PRIORITY'] >= 0) & np.in1d(new.unique_id, old.unique_id[old['PRIORITY']<0])

print("%i exposures no longer considered done"%undone.sum())
for b in BANDS:
    print("  %s: %i "%(b,(undone & (new['FILTER']==b)).sum()))

plt.figure()
smap = skymap.SurveyMcBryde()
smap.draw_fields(new[undone])
smap.draw_des()
plt.title('Exposures no longer considered done')
plt.show()

# This block just updates based on previous field list
if False:
    # Write here
    out = DelveFieldArray.load(args.old)

    ### There are two ways of doing this that should give the same answers...
    print("Running DelveSurvey.update_covered_fields...")
    update = obztak.delve.DelveSurvey.update_covered_fields(old)
    # newly done
    done = (update['PRIORITY'] < 0) & (old['PRIORITY'] >= 0)
    # observed by delve
    obs = np.in1d(update.unique_id,db.unique_id)

    plt.figure()
    smap = skymap.SurveyMcBryde()
    smap.draw_fields(update[done & ~obs])
    smap.draw_des()
    #smap.draw_fields(update[done])
    plt.title('Exposures updated in output')

    print("Writing %s..."%args.outfile)
    update.write(args.outfile)

    # double check
    assert len(fileio.read_csv(args.old)) == len(fileio.read_csv(args.outfile))

print("REMINDER: gzip the output file and move to data directory.")

# Check with DEROSITAS
if False:
    deros_g = DelveFieldArray.load('derositas/data/test_g_all_20210127.csv')
    deros_i = DelveFieldArray.load('derositas/data/test_i_all_20210127.csv')
    plt.figure()
    smap = skymap.SurveyMcBryde()
    smap.draw_fields(deros_g)
    smap.draw_fields(deros_i)
    smap.draw_des()
    #smap.draw_fields(update[done])
    plt.title('DeROSITAS TODO')

# Set missing DeROSITAS field priority
if True:
    import pandas as pd
    uid = pd.read_csv('derositas/derositas-uid-20210615.csv')
    sel = np.in1d(new.unique_id,uid)
    new['PRIORITY'][sel] = 1
    new.write(args.outfile)
