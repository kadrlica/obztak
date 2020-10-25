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

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('old')
parser.add_argument('new')
parser.add_argument('-o','--outfile',default='update_target_fields.csv')
args = parser.parse_args()

db = DelveFieldArray.load_database()
old = DelveFieldArray.load(args.old)
new = DelveFieldArray.load(args.new)

print("Running comparing to new fields...")

if len(old) != len(new):
    print("Different number of fields")

delve = np.in1d(new.unique_id,db.unique_id)

#done = (new['PRIORITY'] < 0) & (old['PRIORITY'] >= 0)
done = (new['PRIORITY'] < 0) & np.in1d(new.unique_id, old.unique_id[old['PRIORITY'] >= 0])

plt.figure()
smap = skymap.SurveyMcBryde()
smap.draw_fields(new[done & ~delve])
smap.draw_des()
plt.title('New')
plt.show()

# Write here
out = DelveFieldArray.load(args.old)

### There are two ways of doing this that should give the same answers...
print("Running DelveSurvey.update_covered_fields...")
update = obztak.delve.DelveSurvey.update_covered_fields(old)
done = (update['PRIORITY'] < 0) & (old['PRIORITY'] >= 0)
delve = np.in1d(update.unique_id,db.unique_id)

plt.figure()
smap = skymap.SurveyMcBryde()
smap.draw_fields(update[done & ~delve])
#smap.draw_fields(update[done])
plt.title('Update')

print("Writing %s..."%args.outfile)
update.write(args.outfile)

# double check
assert len(fileio.read_csv(args.old)) == len(fileio.read_csv(args.outfile))

print("REMINDER: gzip the output file and move to data directory.")
