import os
from mpl_toolkits.basemap import Basemap
import numpy as np
import ephem
import matplotlib.pyplot as plt
import matplotlib
import time
import logging

import maglites.utils.projector
import constants

plt.ion()

############################################################

params = {
    #'backend': 'eps',
    'axes.labelsize': 16,
    #'text.fontsize': 12,           
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    #'text.usetex': True,       # ADW: Slow and no reason for tex right now
    #'font.family':'serif',
    #'font.serif':'Computer Modern Roman',
    #'figure.figsize': fig_size,
    'font.size': 12
    }
matplotlib.rcParams.update(params)

############################################################

#FIGSIZE = (10.5,8.5)
#FIGSIZE = (10.5 / 2., 8.5 / 2.)
#SCALE = np.sqrt((8.0*6.0)/(FIGSIZE[0]*FIGSIZE[1]))
#DPI = 80

#LMC_RA = 80.8939   
#LMC_DEC = -69.7561

############################################################

def safeProj(proj, lon, lat):
    """ Remove points outside of projection """
    x, y = proj(np.asarray(lon),np.asarray(lat))
    x[x > 1e29] = None
    y[y > 1e29] = None
    #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
    return x, y

############################################################

def drawDES(basemap):
    infile = '%s/maglites/data/round13-poly.txt'%(os.environ['MAGLITESDIR'])
    reader_poly = open(infile)
    lines_poly = reader_poly.readlines()
    reader_poly.close()

    ra_poly = []
    dec_poly = []
    for line in lines_poly:
        if '#' in line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        ra_poly.append(float(parts[0]))
        dec_poly.append(float(parts[1]))

    l_poly, b_poly = maglites.utils.projector.celToGal(ra_poly, dec_poly)

    proj = safeProj(basemap, ra_poly, dec_poly)
    basemap.plot(*proj, color='red', lw=2)

############################################################

def drawSMASH(basemap, s=50):
    # SMASH fields
    infile = '%s/maglites/data/smash_fields_final.txt'%(os.environ['MAGLITESDIR'])
    reader = open(infile)
    lines = reader.readlines()
    reader.close()

    ra_smash = []
    dec_smash = []
    for ii in range(0, len(lines)):
        if '#' in lines[ii]:
            continue
        lines[ii] = lines[ii].replace('\xc2\xa0', '')
        parts = lines[ii].split()
        ra_smash.append(float(parts[4]))
        dec_smash.append(float(parts[5]))

    ra_smash = np.array(ra_smash)
    dec_smash = np.array(dec_smash)

    proj = safeProj(basemap, ra_smash, dec_smash)
    basemap.scatter(*proj, edgecolor='black', color='none', marker='h', s=s)

    #basemap.scatter(ra_smash, dec_smash, latlon=True, edgecolor='black', color='none', marker='h', s=50)


############################################################

def drawMAGLITES(basemap):
    infile = '%s/maglites/data/maglites-poly.txt'%(os.environ['MAGLITESDIR'])
    reader_poly = open(infile)
    lines_poly = reader_poly.readlines()
    reader_poly.close()

    ra_poly = []
    dec_poly = []
    for line in lines_poly:
        if '#' in line:
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        ra_poly.append(float(parts[0]))
        dec_poly.append(float(parts[1]))

    l_poly, b_poly = maglites.utils.projector.celToGal(ra_poly, dec_poly)

    proj = safeProj(basemap, ra_poly, dec_poly)
    basemap.plot(*proj, color='blue', lw=2)

############################################################

def drawAirmassContour(basemap, observatory, airmass, n=360, s=50):
    #airmass = 1. / cos(90. - altitude)
    #90 - alt = arccos(1. / airmass)
    altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)

    ra_contour = np.zeros(n)
    dec_contour = np.zeros(n)
    for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, n)):
        ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
        ra_contour[ii] = np.degrees(ra_radians)
        dec_contour[ii] = np.degrees(dec_radians)
    proj = safeProj(basemap, ra_contour, dec_contour)
    basemap.plot(*proj, color='green', lw=2)

    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)
    proj = safeProj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))
    basemap.scatter(*proj, color='green', edgecolor='none', s=s)

############################################################

def drawMoon(basemap, date):
    moon = ephem.Moon()
    moon.compute(date)
    ra_moon = np.degrees(moon.ra)
    dec_moon = np.degrees(moon.dec)

    proj = safeProj(basemap, np.array([ra_moon]), np.array([dec_moon]))

    if np.isnan(proj[0]).all() or np.isnan(proj[1]).all(): return

    basemap.scatter(*proj, color='%.2f'%(0.01 * moon.phase), edgecolor='black', s=500)
    color = 'black' if moon.phase > 50. else 'white'
    plt.text(proj[0], proj[1], '%.2f'%(0.01 * moon.phase), 
             fontsize=10, ha='center', va='center', color=color)

############################################################

def datestring(date): 
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],4)
    minutes = datetuple[-2]
    minutes += seconds//60
    seconds = seconds%60.
    return '%s/%s/%s %02i:%02i:%07.4f'%(datetuple[:-2]+(minutes,seconds))

############################################################

def makePlot(date=None, name=None, figsize=(10.5,8.5), dpi=80, s=50, center=None, airmass=True, moon=True, des=True, smash=True, maglites=True):
    """
    Create map in orthographic projection
    """
    if date is None: date = ephem.now()
    if type(date) != ephem.Date:
        date = ephem.Date(date)

    observatory = ephem.Observer()
    observatory.lon = constants.LON_CTIO
    observatory.lat = constants.LAT_CTIO
    observatory.elevation = constants.ELEVATION_CTIO
    observatory.date = date
    
    #fig, ax = plt.subplots(fig='ortho', figsize=FIGSIZE, dpi=DPI)
    #fig = plt.figure('ortho')
    #ax = plt.subplots(figure=fig, figsize=FIGSIZE, dpi=DPI)
    fig = plt.figure(name, figsize=figsize, dpi=dpi)

    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)

    # Zenith position
    #lon_zen = LMC_RA; lat_zen = LMC_DEC
    lon_zen = ra_zenith; lat_zen = dec_zenith

    # Create the basemap
    proj_kwargs = dict(projection='ortho', celestial=True)
    if center is None:
        lon_0, lat_0 = -lon_zen, lat_zen # Center position
    else:
        lon_0, lat_0 = center[0], center[1]

    proj_kwargs.update(lon_0=lon_0, lat_0=lat_0)
    basemap = Basemap(**proj_kwargs)

    parallels = np.arange(-90.,120.,30.)
    basemap.drawparallels(parallels)
    meridians = np.arange(0.,420.,60.)
    basemap.drawmeridians(meridians)

    if des:   drawDES(basemap)
    if smash: drawSMASH(basemap, s=s)
    if maglites: drawMAGLITES(basemap)
    if airmass: drawAirmassContour(basemap, observatory, 2., s=s)
    if moon: drawMoon(basemap, date)
    plt.title('%s UTC'%(datestring(date)))

    #return fig, ax, basemap
    return fig, basemap

def plotFields(fields, step = 1):
    # ADW: Need to be careful about the size of the marker. It
    # does not change with the size of the frame so it is
    # really safest to scale to the size of the zenith circle
    # (see PlotPointings). That said, s=50 is probably roughly ok.

    ra,dec = fields['RA'],fields['DEC']

    completed = []

    kwargs = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    for i,f in enumerate(fields):
        msg = "%s: time=%s, "%(fields['ID'][i],f['DATE'].split(' ')[-1])
        msg +="ra=%(RA)-8.4f, dec=%(DEC)-8.4f, exptime=%(EXPTIME).0f"%f
        logging.info(msg)
        if plt.get_fignums(): plt.cla()
        fig, basemap = maglites.utils.ortho.makePlot(f['DATE'],name='ortho')
        
        proj = maglites.utils.ortho.safeProj(basemap, ra, dec)
        basemap.scatter(*proj, c=np.zeros(len(ra)), **kwargs)
    
        proj = maglites.utils.ortho.safeProj(basemap, ra[:i+1], dec[:i+1])
        basemap.scatter(*proj, c=fields['TILING'][:i+1], **kwargs)

        # Draw colorbar in existing axis
        if len(fig.axes) == 2:
            colorbar = plt.colorbar(cax=fig.axes[-1])
            colorbar
        else:
            colorbar = plt.colorbar()
        colorbar.set_label('Tiling')

        # Show the selected field
        proj = maglites.utils.ortho.safeProj(basemap, ra[i:i+1], dec[i:i+1])
        basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)
            
        plt.draw()
        time.sleep(0.01)
    

############################################################

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

