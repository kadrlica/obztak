import os
from mpl_toolkits.basemap import Basemap
import numpy as np
import ephem
import matplotlib.pyplot as plt
import matplotlib
import time
import logging

import maglites.utils.projector
import maglites.utils.constants as constants

from maglites.field import FieldArray

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

def plotField(field, target_fields=None, completed_fields=None, **kwargs):
    """
    Plot a specific target field.
    """
    defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    for k,v in defaults.items():
        kwargs.setdefault(k,v)

    if isinstance(field,np.core.records.record):
        tmp = FieldArray(1)
        tmp[0] = field
        field = tmp

    msg = "  Plotting -- "
    msg += "%s (time=%.8s, "%(field['ID'][0],field['DATE'][0].split(' ')[-1])
    msg +="ra=%(RA)-8.4f, dec=%(DEC)-8.4f)"%field[0]
    logging.info(msg)

    if plt.get_fignums(): plt.cla()

    fig, basemap = maglites.utils.ortho.makePlot(field['DATE'][0],name='ortho')
        
    # Plot target fields
    if target_fields is not None:
        proj = maglites.utils.ortho.safeProj(basemap, target_fields['RA'], target_fields['DEC'])
        basemap.scatter(*proj, c=np.zeros(len(target_fields)), **kwargs)

    # Plot completed fields        
    if completed_fields is not None:
        proj = maglites.utils.ortho.safeProj(basemap,completed_fields['RA'],completed_fields['DEC'])
        basemap.scatter(*proj, c=completed_fields['TILING'], **kwargs)

    # Draw colorbar in existing axis
    if len(fig.axes) == 2:
        colorbar = plt.colorbar(cax=fig.axes[-1])
    else:
        colorbar = plt.colorbar()
    colorbar.set_label('Tiling')

    # Show the selected field
    proj = maglites.utils.ortho.safeProj(basemap, field['RA'], field['DEC'])
    basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)
            
    plt.draw()


def plotFields(fields=None,target_fields=None,completed_fields=None, **kwargs):
    # ADW: Need to be careful about the size of the marker. It
    # does not change with the size of the frame so it is
    # really safest to scale to the size of the zenith circle
    # (see PlotPointings). That said, s=50 is probably roughly ok.
    if fields is None:
        fields = completed_fields[-1]

    if isinstance(fields,np.core.records.record):
        tmp = FieldArray(1)
        tmp[0] = fields
        fields = tmp

    for i,f in enumerate(fields):
        plotField(fields[i],target_fields,completed_fields,**kwargs)

        if completed_fields is None: completed_fields = FieldArray(0)
        completed_fields = completed_fields + fields[[i]]

        time.sleep(0.01)
    
def plotWeight(field, target_fields, weight, **kwargs):
    if isinstance(field,FieldArray):
        field = field[-1]

    date = ephem.Date(field['DATE'])

    if plt.get_fignums(): plt.cla()
    fig, basemap = maglites.utils.ortho.makePlot(date,name='weight')
    
    index_sort = np.argsort(weight)[::-1]
    proj = maglites.utils.ortho.safeProj(basemap, target_fields['RA'][index_sort], target_fields['DEC'][index_sort])
    weight_min = np.min(weight)
    basemap.scatter(*proj, c=weight[index_sort], edgecolor='none', s=50, vmin=weight_min, vmax=weight_min + 300., cmap='Spectral')

    #cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
    #proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
    #basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
    
    """
    cut_accomplished = np.in1d(self.target_fields['ID'],self.accomplished_fields['ID'])
    proj = maglites.utils.ortho.safeProj(basemap, 
                                         self.target_fields['RA'][~cut_accomplished], 
                                         self.target_fields['DEC'][~cut_accomplished])
    basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    
    proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
    basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    """

    # Draw colorbar in existing axis
    if len(fig.axes) == 2:
        colorbar = plt.colorbar(cax=fig.axes[-1])
    else:
        colorbar = plt.colorbar()
    colorbar.set_label('Weight')
    
    # Show the selected field
    proj = maglites.utils.ortho.safeProj(basemap, [field['RA']], [field['DEC']])
    basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

    plt.draw()


############################################################

def datestring(date,precision=4): 
    """
    Convert an ephem.Date object to a string with increased precision
    
    Parameters:
    -----------
    date      : ephem.Date object
    precision : Output precision 
    
    Returns:
    --------
    datestr   : String representation of the date
    """
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],precision)
    minutes = datetuple[-2]
    minutes += seconds//60
    seconds = seconds%60.
    datestr = '%d/%02d/%d %02i:%02i:%07.4f'%(datetuple[:-2]+(minutes,seconds))
    return datestr

############################################################

def nite2utc(nite, observer=None):
    import dateutil.parser
    import dateutil.tz
    import datetime
    sun = ephem.Sun()

    if observer is None:
        observer = ephem.Observer()
        observer.lon = constants.LON_CTIO
        observer.lat = constants.LAT_CTIO
        observer.elevation = constants.ELEVATION_CTIO
    
    if not isinstance(nite,datetime.datetime):
        nite = dateutil.parser.parse(nite)
    nite = nite.replace(hour=12,tzinfo=dateutil.tz.tzlocal())
    utc = ephem.Date(nite - nite.utcoffset())
    observer.date = utc   
   
    return observer.next_antitransit(sun)

def utc2nite(utc, observer=None):
    sun = ephem.Sun()

    if observer is None:
        observer = ephem.Observer()
        observer.lon = constants.LON_CTIO
        observer.lat = constants.LAT_CTIO
        observer.elevation = constants.ELEVATION_CTIO

    observer.date = utc

    if observer.previous_setting(sun) > observer.previous_rising(sun):
        # It's night time, use the date of the previous setting
        nite = ephem.localtime(observer.previous_setting(sun))
    else:
        # It's daytime, use the next setting
        nite = ephem.localtime(observer.next_setting(sun))

    return ephem.Date(nite)

def get_nite(utc=None):
    """Convert from a date and time to the 'nite'. 

    A 'nite' is defined by the day (UTC) at noon local time in Chile
    before observing started. This follows the usual convention of
    naming the nite after the day (local time) on which it starts.

    Parameters:
    -----------
    utc : The date/time (UTC) to calculate the nite from.
    
    Returns:
    --------
    nite : An ephem.Date object containing the nite (at sunset)

    """
    if not utc: utc = ephem.now()
    return utc2nite(utc)


############################################################

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

