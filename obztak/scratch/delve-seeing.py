#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import pandas as pd
import numpy as np
import pylab as plt
from scipy.interpolate import interp1d
from collections import OrderedDict as odict

from obztak.utils import fileio
from obztak.delve import DelveTactician, DelveFieldArray
SEXB_RA, SEXB_DEC = 150.00, 5.33
WAVETRANS = odict([
    ( 'g' , 0.9067 ),
    ( 'r' , 0.9609 ),
    ( 'i' ,    1.0 ),
    ( 'z' ,  1.036 ),
    ( 'Y' , 1.0523 ),
])
DECAM2DIMM = 1.0916
DECAMINST = 0.5 # instrumental PSF

def dimm2decam(dimm_fwhm, airmass=1.0, band='i', fwhm_inst = 0.5):
    dimm_fwhm = np.atleast_1d(dimm_fwhm)
    airmass = np.atleast_1d(airmass)
    const = 1/WAVETRANS[band] * 1/DECAM2DIMM
    if dimm_fwhm.ndim == airmass.ndim:
        fwhm = const * (dimm_fwhm*(airmass**(0.6)))
    else:
        fwhm = const * (dimm_fwhm[:,np.newaxis]*(airmass**(0.6)))
    return np.hypot( fwhm, fwhm_inst)

filename = fileio.get_datafile("delve-windows.csv")
windows = pd.read_csv(filename,comment='#',parse_dates=['UTC_START','UTC_END'])

filename = fileio.get_datafile('delve-target-fields.csv')
fields = DelveFieldArray.read(filename)
sel = ((fields['PROGRAM'] == 'delve-deep') & (fields['DEC'] > -5))
sel = [38428]
sel = (fields['TILING'] == 3)
f = fields[sel]
# Sex B
sexb = f[(f['HEX'] == 15854)&(f['FILTER']=='i')]
idx = [9409]
f = f[idx]

tac = DelveTactician()
tac.set_target_fields(f)

times,airmass,moon_phase,moon_alt = [],[],[],[]
delta = pd.Timedelta('5 minutes')
for idx,(start,end) in windows.iterrows():
    t = start
    while t < end:
        tac.set_date(t)
        times.append(t)
        airmass.append(tac.airmass)
        moon_phase.append(tac.moon.phase)
        moon_alt.append(tac.moon.alt)
        t += delta

times = pd.DatetimeIndex(times)
airmass = np.array(airmass)
moon_phase = np.array(moon_phase)
moon_alt = np.array(moon_alt)

# Some configurables
i = 0 # for single field
amax = 1.3
seli = (airmass[:,i] < amax)
selg = (airmass[:,i] < amax) & ( (moon_phase < 10) | (moon_alt < -0.04) )
sel = seli

yval_g = 35.0 * 1.3
yval_i = 23.3 * 1.3
yval = yval_g #+ yval_i

samples = [
    ('dimm_fwhm','k','DIMM (zenith)'),
    ('decam_fwhm_i','blue','DECam (i,zenith)'),
    ('field_fwhm_i','m','DECam (i,field)'),
    #('field_fwhm_g','g','DECam (g,field)'),
]

bins = np.linspace(0,2,50)
hist = dict(
    bins = bins,
    centers = (bins[:-1]+bins[1:])/2.,
    dimm_fwhm = [],
    decam_fwhm_i = [],
    field_fwhm_i = [],
    field_fwhm_g = [],
)

for n in range(1,11):
    dirname = '/Users/kadrlica/delve/observing/data'
    filename =os.path.join(dirname,'delve_sim_%02d.csv'%n)
    print("Reading %s..."%filename)
    df = pd.read_csv(filename,names=['date','fwhm'],
                     parse_dates=['date'],index_col='date')

    interp = interp1d(df.index.to_julian_date(),df['fwhm'],
                      bounds_error=False,fill_value='extrapolate',
                      assume_sorted=True)

    dimm_fwhm = interp(times.to_julian_date())
    print("Calculating DECam FWHM i-band zenith...")
    decam_fwhm_i = dimm2decam(dimm_fwhm,airmass=1.0,band='i')
    print("Calculating field g-band FWHM ...")
    field_fwhm_g = dimm2decam(dimm_fwhm,airmass=airmass,band='g')
    #print("Calculating field r-band FWHM ...")
    #field_fwhm_r = dimm2decam(dimm_fwhm,airmass=airmass,band='r')
    print("Calculating field i-band FWHM ...")
    field_fwhm_i = dimm2decam(dimm_fwhm,airmass=airmass,band='i')

    kwargs = dict(bins=bins)
    hist['dimm_fwhm'].append(np.histogram(dimm_fwhm[sel],**kwargs)[0])
    hist['decam_fwhm_i'].append(np.histogram(decam_fwhm_i[sel],**kwargs)[0])
    hist['field_fwhm_i'].append(np.histogram(field_fwhm_i[:,i][sel],**kwargs)[0])
    hist['field_fwhm_g'].append(np.histogram(field_fwhm_g[:,i][sel],**kwargs)[0])

tstep = delta.seconds/3600. # Convert to hours
hist['dimm_fwhm']    = np.cumsum(np.array(hist['dimm_fwhm']),axis=1)*tstep
hist['decam_fwhm_i'] = np.cumsum(np.array(hist['decam_fwhm_i']),axis=1)*tstep
hist['field_fwhm_i'] = np.cumsum(np.array(hist['field_fwhm_i']),axis=1)*tstep
hist['field_fwhm_g'] = np.cumsum(np.array(hist['field_fwhm_g']),axis=1)*tstep



plt.figure()
plt.axhline(yval,ls='--',c='gray',lw=2)
for i,(n,c,l) in enumerate(samples):
    plt.fill_between(hist['centers'],color=c,alpha=0.2,
                     y1=hist[n].min(axis=0),
                     y2=hist[n].max(axis=0),
    )
    plt.plot(hist['centers'],np.median(hist[n],axis=0),color=c,label=l)

    xval = interp1d(np.median(hist[n],axis=0),hist['centers'])(yval)
    plt.axvline(xval,ls='--',c=c,lw=1)
    plt.annotate('fwhm=%.2f'%xval,(0.6,0.05+0.05*i), xycoords='axes fraction',
                 ha='left',va='top',color=c)

plt.xlim(0.0,2.0);
plt.legend(loc='upper left')
plt.xlabel('FWHM (arcsec)')
plt.ylabel('Available Time (hours)')

"""
kwargs = dict(cumulative=True,histtype='step',lw=2)
#plot airmass
plt.figure()
data = np.nan_to_num(airmass[:,idx])
weights = delta.seconds/3600. * np.ones_like(data) # hours
plt.hist(data,bins=np.linspace(0,2,50), weights=weights,
         color='blue',**kwargs);
plt.xlim(0.0,2.0);
plt.xlabel('Airmass')
plt.ylabel('Available Time (hours)')

#plot airmass
plt.figure()
yval = 75.8 #58.3
sel = (airmass[:,idx] < 1.4)
data = np.nan_to_num(dimm_fwhm[sel])
weights = delta.seconds/3600. * np.ones_like(data) # hours
n,b,p = plt.hist(data,bins=np.linspace(0,2,50), weights=weights, label='DIMM (zenith)',
         color='blue',**kwargs);
xval1 = interp1d(n,(b[1:]+b[:-1])/2.)(yval)

data = np.nan_to_num(decam_fwhm_i[sel])
weights = delta.seconds/3600. * np.ones_like(data) # hours
n,b,p = plt.hist(data,bins=np.linspace(0,2,50), weights=weights, label='DECam (i,zenith)',
         color='red',**kwargs);
xval2 = interp1d(n,(b[1:]+b[:-1])/2.)(yval)

data = np.nan_to_num(field_fwhm_i[:,idx][sel])
weights = delta.seconds/3600. * np.ones_like(data) # hours
n,b,p = plt.hist(data,bins=np.linspace(0,2,50), weights=weights, label='DECam (i,field)',
         color='green',**kwargs);
xval3 = interp1d(n,(b[1:]+b[:-1])/2.)(yval)

plt.axhline(yval,ls='--',c='gray',lw=2)
plt.axvline(xval1,ls='--',c='blue',lw=2)
plt.axvline(xval2,ls='--',c='red',lw=2)
plt.axvline(xval3,ls='--',c='green',lw=2)
plt.annotate('fwhm=%.2f'%xval1,(1.1*xval3,0.9*yval),ha='left',va='top',color='b')
plt.annotate('fwhm=%.2f'%xval2,(1.1*xval3,0.8*yval),ha='left',va='top',color='r')
plt.annotate('fwhm=%.2f'%xval3,(1.1*xval3,0.7*yval),ha='left',va='top',color='g')

plt.xlim(0.0,2.0);
plt.legend(loc='bottom right')
plt.xlabel('FWHM (arcsec)')
plt.ylabel('Available Time (hours)')


plt.fill_between(hist['centers'],color='g',alpha=0.2,
                 y1=hist['decam_fwhm_i'].min(axis=0),
                 y2=hist['decam_fwhm_i'].max(axis=0),
             )
plt.plot(hist['centers'],np.median(hist['decam_fwhm_i'],axis=0),color='g',
         label='DECam (i,zenith)')

plt.fill_between(hist['centers'],color='r',alpha=0.2,
                 y1=hist['field_fwhm_i'].min(axis=0),
                 y2=hist['field_fwhm_i'].max(axis=0),
             )
plt.plot(hist['centers'],np.median(hist['field_fwhm_i'],axis=0),color='r',
         label='DECam (i,field)')


"""
