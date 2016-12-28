import os
from os.path import expandvars
import shutil
from mpl_toolkits.basemap import Basemap
import numpy as np
import ephem
import matplotlib.pyplot as plt
import matplotlib
import time
import logging
import tempfile
import subprocess

import maglites.utils.projector
import maglites.utils.constants as constants
from maglites.utils import fileio
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

class DECamBasemap(Basemap):

    def __init__(self, *args, **kwargs):
        super(DECamBasemap,self).__init__(self,*args,**kwargs)

    def proj(self,lon,lat):
        """ Remove points outside of projection """
        x, y = self(np.atleast_1d(lon),np.atleast_1d(lat))
        x[x > 1e29] = None
        y[y > 1e29] = None
        #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
        return x, y

    def draw_polygon(self,filename,**kwargs):
        """ Draw a polygon footprint on this Basemap instance.
        """
        defaults=dict(color='k', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        perim = np.loadtxt(filename,dtype=[('ra',float),('dec',float)])
        xy = self.proj(perim['ra'],perim['dec'])
        self.plot(*xy,**kwargs)

    def draw_maglites(self,**kwargs):
        defaults=dict(color='blue', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        filename = os.path.join(fileio.get_datadir(),'maglites-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_des(self,**kwargs):
        """ Draw the DES footprint on this Basemap instance.
        """
        defaults=dict(color='red', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        filename = os.path.join(fileio.get_datadir(),'round13-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_smash(self,**kwargs):
        defaults=dict(facecolor='none',color='k')
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        filename = os.path.join(fileio.get_datadir(),'smash_fields_final.txt')
        smash=np.genfromtxt(filename,dtype=[('ra',float),('dec',float)],usecols=[4,5])
        xy = self.proj(smash['ra'],smash['dec'])
        self.scatter(*xy,**kwargs)

    def draw_airmass(self, observatory, airmass, npts=360, **kwargs):
        defaults = dict(color='green', lw=2)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)
        ra_contour = np.zeros(npts)
        dec_contour = np.zeros(npts)
        for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, npts)):
            ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
            ra_contour[ii] = np.degrees(ra_radians)
            dec_contour[ii] = np.degrees(dec_radians)
        xy = self.proj(ra_contour, dec_contour)
        self.plot(*xy, **kwargs)
         
        self.drawZenith(observatory)

    def draw_zenith(self, observatory):
        """
        Plot a to-scale representation of the focal plane size at the zenith.
        """
        defaults = dict(color='green',alpha=0.75,lw=1.5)
        for k,v in defaults.items():
            kwargs.setdefault(k,v)

        # RA and Dec of zenith
        ra_zenith, dec_zenith = np.degrees(observatory.radec_of(0, '90'))
        xy = self.proj(ra_zenith, dec_zenith)
         
        self.plot(*xy,marker='+',ms=10,mew=1.5, **kwargs)
        self.tissot(ra_zenith, dec_zenith, constants.DECAM, 100, fc='none',**kwargs)


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

def safeProj(proj, lon, lat):
    """ Remove points outside of projection """
    x, y = proj(np.atleast_1d(lon),np.atleast_1d(lat))
    x[x > 1e29] = None
    y[y > 1e29] = None
    #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
    return x, y

############################################################

def drawDES(basemap, color='red'):
    infile = os.path.join(fileio.get_datadir(),'round13-poly.txt')
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
    basemap.plot(*proj, color=color, lw=3)

############################################################

def drawSMASH(basemap, color='none', edgecolor='black', marker='h', s=50):
    # SMASH fields
    infile = os.path.join(fileio.get_datadir(),'smash_fields_final.txt')
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
    basemap.scatter(*proj, edgecolor=edgecolor, color=color, marker=marker, s=s)

    #basemap.scatter(ra_smash, dec_smash, latlon=True, edgecolor='black', color='none', marker='h', s=50)


############################################################

def drawMAGLITES(basemap, color='blue'):
    infile = os.path.join(fileio.get_datadir(),'maglites-poly.txt')
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
    basemap.plot(*proj, color=color, lw=3)

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

    drawZenith(basemap, observatory)
    #ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    #ra_zenith = np.degrees(ra_zenith)
    #dec_zenith = np.degrees(dec_zenith)
    #proj = safeProj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))
    #basemap.scatter(*proj, color='green', edgecolor='none', s=s)

def drawZenith(basemap, observatory):
    """
    Plot a to-scale representation of the focal plane size at the zenith.
    """
    ra_zenith, dec_zenith = observatory.radec_of(0, '90') # RA and Dec of zenith
    ra_zenith = np.degrees(ra_zenith)
    dec_zenith = np.degrees(dec_zenith)
    proj = safeProj(basemap, np.array([ra_zenith]), np.array([dec_zenith]))

    zen_kwargs = dict(color='green',alpha=0.75,lw=1.5,zorder=1000)
    basemap.plot(*proj,marker='+',ms=10,mew=1.5, **zen_kwargs)
    basemap.tissot(ra_zenith, dec_zenith, constants.DECAM, 100, fc='none',**zen_kwargs)


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
    plt.cla()

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
    basemap = DECamBasemap(**proj_kwargs)

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

def plotField(field, target_fields=None, completed_fields=None, options_basemap={}, **kwargs):
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
    msg +="ra=%(RA)-6.2f, dec=%(DEC)-6.2f, secz=%(AIRMASS)-4.2f)"%field[0]
    logging.info(msg)

    #if plt.get_fignums(): plt.cla()

    defaults = dict(date=field['DATE'][0], name='ortho')
    options_basemap = dict(options_basemap)
    for k,v in defaults.items():
        options_basemap.setdefault(k,v)
    #fig, basemap = maglites.utils.ortho.makePlot(field['DATE'][0],name='ortho',**options_basemap)
    fig, basemap = maglites.utils.ortho.makePlot(**options_basemap)
   
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

    #plt.draw()
    plt.pause(0.001)

def plotFields(fields=None,target_fields=None,completed_fields=None,options_basemap={},**kwargs):
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
        plotField(fields[i],target_fields,completed_fields,options_basemap,**kwargs)
        #plt.savefig('field_%08i.png'%i)
        if completed_fields is None: completed_fields = FieldArray(0)
        completed_fields = completed_fields + fields[[i]]

        #time.sleep(0.01)

def movieFields(outfile,fields=None,target_fields=None,completed_fields=None,**kwargs):
    if os.path.splitext(outfile)[-1] not in ['.gif']:
        msg = "Only animated gif currently supported."
        raise Exception(msg)

    tmpdir = tempfile.mkdtemp()

    if fields is None:
        fields = completed_fields[-1]

    if isinstance(fields,np.core.records.record):
        tmp = FieldArray(1)
        tmp[0] = fields
        fields = tmp

    for i,f in enumerate(fields):
        plotField(fields[i],target_fields,completed_fields,**kwargs)
        png = os.path.join(tmpdir,'field_%08i.png'%i)
        plt.savefig(png,bbox_inches='tight',dpi=100)
        if completed_fields is None: completed_fields = FieldArray(0)
        completed_fields = completed_fields + fields[[i]]

    cmd = 'convert -delay 10 -loop 0 %s/*.png %s'%(tmpdir,outfile)
    logging.info(cmd)
    subprocess.call(cmd,shell=True)
    shutil.rmtree(tmpdir)
    return outfile

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

    #plt.draw()
    plt.pause(0.001)
    #fig.canvas.draw()

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
    """
    date = ephem.Date(date).datetime()
    datestr = date.strftime('%Y/%m/%d %H:%M:%S')
    datestr += '{:.{precision}f}'.format(date.microsecond/1.e6,
                                         precision=precision)[1:]
    """

    if precision < 0:
        msg = "Precision must be positive."
        raise Exception(msg)

    # This is a bit annoying, but works better than converting to datetime
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],precision)
    minutes = datetuple[-2]
    hours = datetuple[-3]
    minutes += int(seconds//60)
    seconds = seconds%60.
    hours += int(minutes//60)
    minutes = minutes%60

    strtuple = datetuple[:-3]+(hours,minutes,seconds)
    width = precision+2 if precision == 0 else precision+3

    #datestr = '%d/%02d/%02d %02i:%02i:%07.4f'%strtuple
    datestr = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:0{width}.{precision}f}'
    datestr = datestr.format(*strtuple,precision=precision,width=width)

    return datestr

############################################################

def nite2utc(nite, observer=None):
    import dateutil.parser
    import dateutil.tz
    import datetime

    if observer is None:
        observer = ephem.Observer()
        observer.lon = constants.LON_CTIO
        observer.lat = constants.LAT_CTIO
        observer.elevation = constants.ELEVATION_CTIO
        
    if not isinstance(nite,datetime.datetime):
        nite = dateutil.parser.parse(str(nite))

    # Move to (local) noon
    # This depends on where the user is and isn't very safe
    nite = nite.replace(hour=12,tzinfo=dateutil.tz.tzlocal())
    utc = ephem.Date(nite - nite.utcoffset())

    # Maybe something like this instead...
    #offset = int( (observer.lon/(2*np.pi)) * 24. * 60) * 60
    #tzinfo = dateutil.tz.tzoffset('OBS',offset)
    #nite = nite.replace(hour=12,tzinfo=tzinfo)
    #utc = ephem.Date(nite - nite.utcoffset())

    observer.date = utc
    observer.horizon = '-14'
    #ret = observer.next_antitransit(ephem.Sun())
    ret = observer.next_setting(ephem.Sun(), use_center=True)

    observer.horizon = '0'
    return ret

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

def get_nite(utc=None, observer=None):
    """Convert from a UTC date and time to the 'nite'. 

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
    return utc2nite(utc, observer)


############################################################

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

