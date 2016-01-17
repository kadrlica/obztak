import os
from mpl_toolkits.basemap import Basemap
import numpy as np
import ephem
import matplotlib.pyplot as plt
import matplotlib

import utils
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
    'text.usetex': True,
    #'figure.figsize': fig_size,
    'font.family':'serif',
    'font.serif':'Computer Modern Roman',
    'font.size': 12
    }

matplotlib.rcParams.update(params)

############################################################

FIGSIZE = (10.5,8.5)
SCALE = np.sqrt((8.0*6.0)/(FIGSIZE[0]*FIGSIZE[1]))
DPI = 80

LMC_RA = 80.8939   
LMC_DEC = -69.7561

############################################################

def safe_proj(proj, lon, lat):
    """ Remove points outside of projection """
    x, y = proj(np.asarray(lon),np.asarray(lat))
    x[x > 1e29] = None
    y[y > 1e29] = None
    return x, y

############################################################

def drawDES(basemap):
    infile = '%s/data/round13-poly.txt'%(os.environ['MAGLITESDIR'])
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

    l_poly, b_poly = utils.celToGal(ra_poly, dec_poly)

    proj = safe_proj(basemap, ra_poly, dec_poly)
    basemap.plot(*proj, color='red', lw=2)

############################################################

def drawSMASH(basemap):
    # SMASH fields
    infile = 'smash_fields_final.txt'
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

    proj = safe_proj(basemap, ra_smash, dec_smash)
    basemap.scatter(*proj, edgecolor='black', color='none', marker='h', s=50)

    #basemap.scatter(ra_smash, dec_smash, latlon=True, edgecolor='black', color='none', marker='h', s=50)

############################################################

def drawAirmassContour(basemap, observatory, airmass, n=360):
    #airmass = 1. / cos(90. - altitude)
    #90 - alt = arccos(1. / airmass)
    altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)

    ra_contour = np.zeros(n)
    dec_contour = np.zeros(n)
    for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, n)):
        ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
        ra_contour[ii] = np.degrees(ra_radians)
        dec_contour[ii] = np.degrees(dec_radians)
    proj = safe_proj(basemap, ra_contour, dec_contour)
    basemap.plot(*proj, color='green', lw=2)

    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)
    proj = safe_proj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))
    basemap.scatter(*proj, color='green', edgecolor='none', s=50)

############################################################

def datestring(date):
    if type(date) != ephem.Date:
        date = ephem.Date(date) 
    return '%s/%s/%s %02i:%02i:%02i'%date.tuple()

############################################################

def makePlot(date):
    
    observatory = ephem.Observer()
    observatory.lon = constants.LON_CTIO
    observatory.lat = constants.LAT_CTIO
    observatory.elevation = constants.ELEVATION_CTIO
    observatory.date = date
    
    #fig, ax = plt.subplots(fig='ortho', figsize=FIGSIZE, dpi=DPI)
    #fig = plt.figure('ortho')
    #ax = plt.subplots(figure=fig, figsize=FIGSIZE, dpi=DPI)
    fig = plt.figure('ortho', figsize=FIGSIZE, dpi=DPI)

    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)

    # Zenith position
    #lon_zen = LMC_RA; lat_zen = LMC_DEC
    lon_zen = ra_zenith; lat_zen = dec_zenith

    # Create the basemap
    proj_kwargs = dict(projection='ortho', celestial=True)
    lon_0, lat_0 = -lon_zen, lat_zen # Center position
    proj_kwargs.update(lon_0=lon_0, lat_0=lat_0)
    basemap = Basemap(**proj_kwargs)
    parallels = np.arange(-90.,120.,30.)
    basemap.drawparallels(parallels)
    meridians = np.arange(0.,420.,60.)
    basemap.drawmeridians(meridians)

    drawDES(basemap)
    drawSMASH(basemap)
    drawAirmassContour(basemap, observatory, 2.)

    plt.title(datestring(date))

    #return fig, ax, basemap
    return fig, basemap

############################################################
