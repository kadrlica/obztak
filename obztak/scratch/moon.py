#!/usr/bin/env python
"""
Output exposure metadata with moon info.
"""
import sys
from collections import OrderedDict as odict

import healpy as hp
import numpy as np
import pylab as plt
import numpy.lib.recfunctions as recfns
import fitsio

import obztak.tactician
from obztak.utils.database import Database

#BANDS = ['u','g','r','i','z','Y','VR']
BANDS = ['g','r','i']

# COALESCE(qc_teff,'NaN')
query ="""
SELECT id as expnum, date, telra as ra, teldec as dec, exptime, filter,
(CASE WHEN moonangl is NULL THEN 'NaN' ELSE moonangl END) as moonangl,
(CASE WHEN qc_teff is NULL THEN 'NaN' WHEN qc_teff=-1 THEN 'NAN' ELSE qc_teff END) as teff,
qc_fwhm::FLOAT, qc_cloud::FLOAT, qc_sky::FLOAT
FROM exposure where
aborted=False and exposed=True and digitized=True and built=True and delivered=True and discard=False
and flavor = 'object' and telra between 0 and 360 and teldec between -90 and 90
and exptime >= 30
and filter in (%s) and propid NOT LIKE '%%-9999'
-- and date < current_date - interval '18 months'
and qc_teff >= 0 and qc_cloud is not NULL
ORDER BY id;
"""%(",".join(["'%s'"%b for b in BANDS]))
#and filter in (%s) and propid = '2012B-0001'

print(query)

db = Database()
db.connect()
data = db.query2recarray(query)

tac = obztak.tactician.Tactician()

print("Calculating moon phase and altitude...")
alt,phase = [],[]
for d in data:
    tac.set_date(d['date'])
    phase.append(tac.moon.phase)
    alt.append(tac.moon.alt)
alt = np.degrees(alt)
data = recfns.rec_append_fields(data,['moonphase','moonalt'],[phase,alt])

# Convert to types that fits understands
dtype = data.dtype.descr
dtype[data.dtype.names.index('date')] = ('date','|S32')

out = np.recarray(len(data),dtype=dtype)
out[:] = data[:]
outfile = 'decam-exposures-moon.fits.gz'
fitsio.write(outfile,out,clobber=True)
