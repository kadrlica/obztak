#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import time

import pylab as plt
import maglites.utils.projector
import numpy as np

from maglites.utils.projector import SphericalRotator
from maglites.utils.projector import angsep, match, footprint
from maglites.utils.ortho import makePlot,safeProj
from maglites.utils.constants import SMASH_POLE
import maglites.utils.fileio as fileio

datadir = fileio.get_datadir()
filename  = os.path.join(datadir,'smash_fields_alltiles.txt')
data = np.recfromtxt(filename,names=True)
d = data

def smash_dither(dx,dy,ra=d['RA'],dec=d['DEC']):
    ra0,dec0 = SMASH_POLE
    # Rotate the pole at SMASH_POLE to (0,-90)
    R1 = SphericalRotator(ra0,90+dec0); 
    ra1,dec1 = R1.rotate(ra,dec)
    # Dither the offset 
    ra2 = ra1+dx/np.cos(np.radians(dec1))
    dec2 = dec1+dy
    # Rotate back to the original frame (keeping R2)
    return R1.rotate(ra2,dec2,invert=True)

def decam_dither(dx,dy,ra=d['RA'],dec=d['DEC']):
    out = []
    for _ra,_dec in zip(ra,dec):
        out.append(SphericalRotator(_ra,_dec).rotate(dx,dy,invert=True))
    return np.array(out).T

def smash_rotate(dx,dy,ra=d['RA'],dec=d['DEC']):
    ra0,dec0 = SMASH_POLE
    # Rotate the pole at SMASH_POLE to (0,-90)
    R1 = SphericalRotator(ra0,90+dec0); 
    ra1,dec1 = R1.rotate(ra,dec)
    # Rotate around the SMASH pole
    R2 = SphericalRotator(dx,dy)
    ra2,dec2 = R2.rotate(ra1,dec1)
    # Rotate back to the original frame (keeping the R2 shift)
    return R1.rotate(ra2,dec2,invert=True)


def plot(ra1,dec1,ra2=d['RA'],dec2=d['DEC'],sel=None):
    if sel is None: sel = slice(None,None)
    sep = angsep(ra1,dec1,ra2,dec2)
    plt.clf(); 
    plt.scatter(ra1[sel],dec1[sel],c=sep[sel],s=5,edgecolor='none'); 
    plt.colorbar()
    plt.draw()
    time.sleep(0.01)
    return sep

def zoom(ax=None):
    if ax is None: ax = plt.gca()
    ax.set_xlim(2.7e6,9.2e6)
    ax.set_ylim(4.1e6,9.1e6)

    #ax.set_xlim(5.1e6,7.6e6)
    #ax.set_ylim(5.1e6,7.6e6)

plt.ion()

#tilings = [(0,0),(0.75,0.75),(-0.75,0.75),(0,0.75)]
#tilings = [(0,0),(0.62,0.62),(-.62,-0.62)]#, (0.7,-0.2)]
tilings = [(0,0),(1.0,0.0),(-1.0,0.0),(0.0,0.75)]
colors = ['k','r','b','g']
kwargs = dict(date='2016/02/11 03:00:00',center=(0,-90),airmass=False,moon=False,smash=False)
sc_kwargs = dict(edgecolor='none',s=50,vmin=0.3,vmax=1.6)
sc_kwargs = dict(edgecolor='none',s=50)

"""
# Make some plots for the original tilings
fig,basemap=makePlot(**kwargs)
ra,dec =d['RA'],d['DEC']

x,y = safeProj(basemap,ra,dec)
basemap.scatter(x,y,facecolor='none',edgecolor='k', s=50)
plt.savefig('smash_hexs.png',bbox_inches='tight')

fig,basemap=makePlot(**kwargs)
R1 = SphericalRotator(SMASH_POLE[0],90+SMASH_POLE[-1]); 
ra1,dec1 = R1.rotate(ra,dec)
basemap.scatter(x,y,c=dec1,edgecolor='none',s=50)
plt.colorbar(label ='SMASH Declination (deg)')
tropic = (dec1 > -30.5) & (dec1 < -29.5)
basemap.scatter(x[tropic],y[tropic],c='k',edgecolor='none',s=50)
plt.savefig('smash_dec.png',bbox_inches='tight')
zoom()
plt.savefig('smash_dec_zoom.png',bbox_inches='tight')

fig,basemap=makePlot(**kwargs)
idx1,idx2,sep = match(ra,dec,ra,dec,nnearest=2)
basemap.scatter(x,y,c=sep,**sc_kwargs)
plt.colorbar(label = 'Minimum Separation (deg)')
plt.savefig('smash_sep.png',bbox_inches='tight')

plt.figure()
plt.hist(sep,bins=np.linspace(0,2,50),lw='2',histtype='step',color='k')
plt.xlabel('Minimum Separation (deg)')
plt.ylabel('Number of Fields')
plt.savefig('smash_hist.png',bbox_inches='tight')
"""

# Make the tilings
ra_list,dec_list = [],[]
for i,tiling in enumerate(tilings):
    _ra,_dec=smash_dither(*tiling,ra=d['RA'],dec=d['DEC'])
    sel = footprint(_ra,_dec)
    ra_list.append(_ra[sel]); dec_list.append(_dec[sel])

"""
# Plot the tilings
ra,dec = [],[]
fig_hist = plt.figure('hist')
for i,(_ra,_dec) in enumerate(zip(ra_list,dec_list)):
    num = i+1
    ra,dec = np.hstack(zip(ra_list[:i+1],dec_list[:i+1]))

    sep = match(ra,dec,ra,dec,nnearest=2)[-1]

    fig,basemap = makePlot(**kwargs)
    x,y = safeProj(basemap,ra,dec)
    basemap.scatter(x, y, c=sep, **sc_kwargs)
    plt.colorbar(label='Minimum Separation (deg)')
    title = 'Tiling %s'%(num)
    plt.title(title)
    plt.savefig('separation_tiling_%i.png'%num,bbox_inches='tight')
    #zoom()
    #plt.savefig('separation_tiling_%i_zoom.png'%num,bbox_inches='tight')
    
    plt.figure('hist')
    plt.hist(sep,bins=np.linspace(0,2,100),normed=True,lw='2',histtype='step',color=colors[i],label='Tiling %s'%num)

plt.figure('hist')
plt.legend(loc='upper right',fontsize=12)
plt.xlabel('Minimum Separation (deg)')
plt.ylabel('Normalized Number of Fields')
plt.savefig('maglites_hist.png',bbox_inches='tight')
"""

fig,basemap = makePlot(**kwargs)
for i,(_ra,_dec) in enumerate(zip(ra_list,dec_list)):
    num = i+1
    _x,_y = safeProj(basemap,_ra,_dec)
    basemap.scatter(_x,_y,facecolor='none',edgecolor=colors[i],s=50,label='Tiling %i'%num)

    plt.legend(loc='lower left',fontsize=10)
    plt.title('Tiling %s'%(i+1))
    plt.savefig('tiling_position_%i.png'%num,bbox_inches='tight')
    zoom()
    plt.savefig('tiling_position_%i_zoom.png'%num,bbox_inches='tight')
    plt.gca().relim()
    plt.gca().autoscale()


fig,basemap = makePlot(**kwargs)
PIXSCALE=0.2626
CCD = [4096*PIXSCALE/3600.,2048*PIXSCALE/3600.]

tilings2 = [(0,0),(8/3.*CCD[0],-11/3.*CCD[1]),(-8/3.*CCD[0],+11/3.*CCD[1])]
tilings2 = [(0,0),(8/3.*CCD[0],-11/3.*CCD[1]),(+8/3.*CCD[0],+11/3.*CCD[1])]
tilings2 = [(0,0),(8/3.*CCD[0],-11/3.*CCD[1]),(+5/3.*CCD[0],+8/3.*CCD[1])]

sel = footprint(d['RA'],d['DEC'])
ra,dec = d['RA'][sel],d['DEC'][sel]

ra_list2,dec_list2 = [],[]
for i,tiling in enumerate(tilings2):
    _ra,_dec=decam_dither(*tiling,ra=ra,dec=dec)
    ra_list2.append(_ra); dec_list2.append(_dec)


for i,(_ra,_dec) in enumerate(zip(ra_list2,dec_list2)):
    num = i+1
    _x,_y = safeProj(basemap,_ra,_dec)
    basemap.scatter(_x,_y,facecolor='none',edgecolor=colors[i],s=50,label='Tiling %i'%num)

    plt.legend(loc='lower left',fontsize=10)
    plt.title('Tiling %s'%(i+1))
    plt.savefig('decam_tiling_position_%i.png'%num,bbox_inches='tight')
    zoom()
    plt.savefig('decam_tiling_position_%i_zoom.png'%num,bbox_inches='tight')
    plt.gca().relim()
    plt.gca().autoscale()




"""
for x2 in np.arange(-90,90,30):
    for y2 in np.arange(-90,90,30):
        for x1 in np.arange(-180,180,30):
            for y1 in np.arange(-90,90,30):
                plot(x1,y1,x2,y2)
"""


#for x1 in [-190]:
#    for y1 in [-90]:
#        for x2 in [-45]:
#            for y2 in [-25]:
#                plot(x1,y1,x2,y2)

