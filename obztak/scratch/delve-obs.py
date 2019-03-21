#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import numpy as np
import pylab as plt
import datetime as dt
import dateutil.parser
from dateutil.parser import parse as dparse
import fitsio
import scipy.stats
from collections import OrderedDict as odict

import obztak.delve
from obztak.utils import fileio
import obztak.tactician
import obztak.survey
from obztak.utils import fileio

# CTIO midpoint varies from 04:40 to 05:00 UTC over the course of the year.

SEMESTERS = odict([
    ('2019A',(dparse('2019-02-02T04:50UTC'),dparse('2019-08-01T04:50UTC'),20)),
    ('2019B',(dparse('2019-08-02T04:50UTC'),dparse('2020-02-01T04:50UTC'),21)),
    ('2020A',(dparse('2020-02-02T04:50UTC'),dparse('2020-08-01T04:50UTC'),22)),
    ('2020B',(dparse('2020-08-02T04:50UTC'),dparse('2021-02-01T04:50UTC'),21)),
    ('2021A',(dparse('2021-02-02T04:50UTC'),dparse('2021-08-01T04:50UTC'),22)),
    ('2021B',(dparse('2021-08-02T04:50UTC'),dparse('2022-02-01T04:50UTC'),20)),
])

def get_semester(date):
    for key,(start,stop,nights) in SEMESTERS.items():
        if (date >= start) and (date <= stop):
            return key
    raise ValueError(str(date))

def choose_2019A(data):
    nights = [
        [dparse('2019-02-07T04:50UTC'), 'second'],
        [dparse('2019-02-08T04:50UTC'), 'second'],
        [dparse('2019-02-09T04:50UTC'), 'second'],
        [dparse('2019-02-12T04:50UTC'), 'full  '],
        [dparse('2019-02-13T04:50UTC'), 'full  '],
        [dparse('2019-02-14T04:50UTC'), 'second'],
        [dparse('2019-02-15T04:50UTC'), 'full  '],
        [dparse('2019-02-24T04:50UTC'), 'second'],
        [dparse('2019-02-25T04:50UTC'), 'second'],
        [dparse('2019-02-26T04:50UTC'), 'second'],
        [dparse('2019-02-27T04:50UTC'), 'second'],
        [dparse('2019-02-28T04:50UTC'), 'second'],
        [dparse('2019-03-01T04:50UTC'), 'second'],
        [dparse('2019-05-12T04:50UTC'), 'full  '],
        [dparse('2019-05-13T04:50UTC'), 'full  '],
        [dparse('2019-05-28T04:50UTC'), 'second'],
        [dparse('2019-05-29T04:50UTC'), 'second'],
        [dparse('2019-05-30T04:50UTC'), 'second'],
        [dparse('2019-05-31T04:50UTC'), 'second'],
        [dparse('2019-06-01T04:50UTC'), 'second'],
        [dparse('2019-06-02T04:50UTC'), 'second'],
        [dparse('2019-06-03T04:50UTC'), 'second'],
        [dparse('2019-06-04T04:50UTC'), 'second'],
        [dparse('2019-06-05T04:50UTC'), 'full  '],
        [dparse('2019-06-06T04:50UTC'), 'full  '],
        [dparse('2019-06-07T04:50UTC'), 'full  '],
        [dparse('2019-06-08T04:50UTC'), 'full  '],
        [dparse('2019-06-09T04:50UTC'), 'full  '],
        [dparse('2019-06-23T04:50UTC'), 'second'],
        [dparse('2019-06-24T04:50UTC'), 'second'],
        [dparse('2019-06-25T04:50UTC'), 'second'],
        [dparse('2019-06-26T04:50UTC'), 'second'],
        [dparse('2019-06-27T04:50UTC'), 'second'],
        [dparse('2019-06-28T04:50UTC'), 'second'],
    ]

    sel = np.in1d(data['date'],[n[0] for n in nights])
    choice = data[sel]
    for n in nights:
        choice['half'][choice['date'] == n[0]] = n[1]

    return choice

def plot_objects():
    yval = plt.ylim()[0] + 0.75*(plt.ylim()[1] - plt.ylim()[0])
    kwargs = dict(rotation=90.,ha='center',va='center')
    plt.axvline(80.8939,ls='--',color='gray')
    plt.annotate('LMC', (75,yval), color='gray', **kwargs)
    plt.axvline(13.1867,ls='--',color='gray')
    plt.annotate('SMC', (8,yval), color='gray', **kwargs)
    plt.axvline(150,ls='--',color='gray')
    plt.annotate('Sex B', (140,yval), color='gray', **kwargs)
    plt.axvline(3.79,ls='--',color='gray')
    plt.annotate('NGC 55', (-2,yval), color='gray', **kwargs)
    plt.axvline(13.72,ls='--',color='gray')
    plt.annotate('NGC 300', (19,yval), color='gray', **kwargs)
    plt.axvline(228.56,ls='--',color='gray')
    plt.annotate('ESO274-001', (223,yval), color='gray', **kwargs)
    plt.axvline(330.67,ls='--',color='gray')
    plt.annotate('IC 5152', (325,yval), color='gray', **kwargs)


begin = SEMESTERS['2019A'][0]
end   = SEMESTERS['2021B'][1]
size = (end - begin).days

dtype = [('semester','S5'),('date',object),('sel',bool),('moon',float),
         ('zenith',float),('weight',float),('half','S6')]
data = np.recarray(size,dtype=dtype)

# Create the array of possible dates
data['date'] = np.array([begin + dt.timedelta(days=i) for i in xrange(size)])
data['half'] = 'full  '
tac = obztak.tactician.Tactician()

print("Calculating moon phase and altitude...")
for i,d in enumerate(data):
    tac.set_date(str(d['date']).replace('T',' ')+'UTC')
    data[i]['moon'] = tac.moon.phase
    data[i]['zenith'] = tac.zenith_angle[0]
    data[i]['semester'] = get_semester(d['date'])

#data = recfns.rec_append_fields(data,['moonphase','moonalt'],[phase,alt])

fields = obztak.delve.DelveFieldArray.read(fileio.get_datafile('delve-target-fields.csv'))
sel  = fields['PRIORITY'] > 0
sel &= ~((fields['TILING'] > 3) & (fields['PROGRAM'] == 'delve-wide'))
sel &= ~((fields['TILING'] > 3) & (fields['PROGRAM'] == 'delve-mc'))
fields = fields[sel]

print("Calculating KDE...")
pts = np.repeat(fields['RA'],np.round(fields['EXPTIME']/90.).astype(int))
pts = np.concatenate([pts, pts[pts < 5.]+360, pts[pts > 355.] - 360])
kde = scipy.stats.gaussian_kde(pts,bw_method=7.0/360.)
#kde = scipy.stats.gaussian_kde(pts)
data['weight'] = kde(data['zenith'])
# Remove exposure close to full moon
data['weight'][data['moon'] > 80.] = 0

print("Making choices...")
np.random.seed(0)
schedule = []
for key,(start,stop,size) in SEMESTERS.items():
    if key == '2019A':
        choice = choose_2019A(data)
    else:
        d = data[data['semester'] == key]
        w = d['weight']/d['weight'].sum()
        idx = np.random.choice(np.arange(len(d)),size,p=w,replace=False)
        choice = d[np.sort(idx)]
    schedule.append(choice)

schedule = np.concatenate(schedule)
filename = 'delve-schedule.csv'
print("Writing schedule %s..."%filename)
fileio.rec2csv(filename,schedule)

print("Creating windows...")
survey = obztak.survey.Survey()
nights = [[d['date'].strftime('%Y/%m/%d'),d['half']] for d in schedule]
outfile = 'delve-windows-predicted.csv'
windows = survey.prepare_windows(nights,outfile=outfile)


################
### Plotting ###
################

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
plt.sca(ax1)

x = np.linspace(0,360.,3600)
kwargs = dict(histtype='step',lw=1.5,bins=np.linspace(-10,370,381),normed=True)
plt.hist(fields['RA'],weights=fields['EXPTIME'],color='k',label='raw',**kwargs)
#plt.hist(pts,color='r',label='points',**kwargs)
plt.plot(x,kde(x),lw=2,c='b',label='kde')
w = 0.5 + 0.5*(schedule['half']=='full  ')
plt.hist(schedule['zenith'],color='g',label='schedule',weights=w, **kwargs)
#plot_objects()

plt.legend()
plt.ylabel("Normalized Observation Time")
plt.xlabel("RA (deg)")
ax1.set_xlim(-30,380)
ax2.set_xlim(ax1.get_xlim())

tickdates = ['2019-01-01T04:50UTC','2019-02-01T04:50UTC','2019-03-01T04:50UTC',
             '2019-04-01T04:50UTC','2019-05-01T04:50UTC','2019-06-01T04:50UTC',
             '2019-07-01T04:50UTC','2019-08-01T04:50UTC','2019-09-01T04:50UTC',
             '2019-10-01T04:50UTC','2019-11-01T04:50UTC','2019-12-01T04:50UTC',
]

ticks = []
ticklabels = []
for d in tickdates:
    date = dparse(d)
    tac.set_date(date)
    ticks.append(tac.zenith_angle[0])
    ticklabels.append(date.strftime('%b %d'))
ax2.set_xticks(ticks)
ax2.set_xticklabels(ticklabels,rotation=20)
ax2.set_xlabel("Midnight Zenith")


"""
g = fields[fields['FILTER'] == 'g']
r = fields[fields['FILTER'] == 'r']
i = fields[fields['FILTER'] == 'i']

kwargs = dict(histtype='step',lw=1.5,bins=np.linspace(0,360,37))
plt.figure()
plt.hist(fields['RA'],weights=fields['EXPTIME']/3600.,color='k',label='all',**kwargs)
plt.hist(g['RA'],weights=g['EXPTIME']/3600.,color='g',label='g',**kwargs)
plt.hist(r['RA'],weights=r['EXPTIME']/3600.,color='r',label='r',**kwargs)
plt.hist(i['RA'],weights=i['EXPTIME']/3600.,color='purple',label='i',**kwargs)
plt.legend()
plt.ylabel("Observation time (hours)")
plt.xlabel("RA (deg)")

kwargs = dict(histtype='step',lw=1.5,bins=np.linspace(0,360,37))
plt.figure()
plt.hist(fields['RA'],weights=fields['EXPTIME']/3600.,color='k',label='all',**kwargs)

for (program,c) in [ ('wide','g'),('mc','r'),('deep','b')]:
    f = fields[np.char.find(fields['PROGRAM'],program) > 0]
    plt.hist(f['RA'],weights=f['EXPTIME']/3600.,label=program,color=c,**kwargs)

plt.legend()
plt.ylabel("Observation time (hours)")
plt.xlabel("RA (deg)")
"""
