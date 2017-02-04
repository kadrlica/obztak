import os
from os.path import expandvars
import shutil
import time
import logging
import tempfile
import subprocess
import warnings
from collections import OrderedDict as odict

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import ephem

import obztak.utils.projector
from obztak.utils.projector import gal2cel, SphericalRotator

from obztak.utils import constants
from obztak.utils import fileio
from obztak.field import FieldArray
from obztak.utils.date import datestring,nite2utc,utc2nite,get_nite
from obztak.ctio import CTIO
from obztak.utils.constants import RA_LMC,DEC_LMC,RADIUS_LMC
from obztak.utils.constants import RA_SMC,DEC_SMC,RADIUS_SMC

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

COLORS = odict([
    ('none','black'),
    ('u','blue'),
    ('g','green'),
    ('r','red'),
    ('i','gold'),
    ('z','magenta'),
    ('Y','black'),
    ('VR','gray'),
])

############################################################

def setdefaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

class DECamBasemap(Basemap):

    def __init__(self, *args, **kwargs):
        super(DECamBasemap,self).__init__(self,*args,**kwargs)
        self.draw_parallels()
        self.draw_meridians()

    def draw_parallels(self,*args,**kwargs):
        defaults = dict()
        if not args: defaults.update(circles=np.arange(-90,120,30))
        setdefaults(kwargs,defaults)
        self.drawparallels(*args, **kwargs)

    def draw_meridians(self,*args,**kwargs):
        defaults = dict(labels=[1,0,0,1])
        if self.projection in ['ortho','geos','nsper','aeqd']:
            defaults.update(labels=[0,0,0,0])
        if not args: defaults.update(meridians=np.arange(0,420,60))
        setdefaults(kwargs,defaults)
        self.drawmeridians(*args,**kwargs)
        
    def proj(self,lon,lat):
        """ Remove points outside of projection """
        x, y = self(np.atleast_1d(lon),np.atleast_1d(lat))
        x[x > 1e29] = None
        y[y > 1e29] = None
        #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
        return x, y

    @staticmethod
    def roll(ra,dec):
        idx = np.abs(ra - 180).argmin()
        if   (ra[idx]<180) and (ra[idx+1]>180): idx += 1
        elif (ra[idx]>180) and (ra[idx+1]<180): idx += 1
        return np.roll(ra,-idx), np.roll(dec,-idx)

    @staticmethod
    def split(ra,angle=180):
        pass

    def draw_polygon(self,filename,**kwargs):
        """ Draw a polygon footprint on this Basemap instance.
        """
        defaults=dict(color='k', lw=2)
        setdefaults(kwargs,defaults)

        perim = np.loadtxt(filename,dtype=[('ra',float),('dec',float)])
        self.draw_polygon_radec(perim['ra'],perim['dec'],**kwargs)

    def draw_polygon_radec(self,ra,dec,**kwargs):
        xy = self.proj(ra,dec)
        self.plot(*xy,**kwargs)
        
    def draw_galaxy(self,width=10,**kwargs):
        defaults = dict(color='k',lw=1.5,ls='-')
        setdefaults(kwargs,defaults)

        glon = np.linspace(0,360,200)
        glat = np.zeros_like(glon)
        ra,dec = self.roll(*gal2cel(glon,glat))

        self.draw_polygon_radec(ra,dec,**kwargs)
        
        if width:
            kwargs.update(dict(ls='--',lw=1))
            for delta in [+width,-width]:
                ra,dec = self.roll(*gal2cel(glon,glat+delta))
                self.draw_polygon_radec(ra,dec,**kwargs)
            
        
    def draw_maglites(self,**kwargs):
        defaults=dict(color='blue', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(fileio.get_datadir(),'maglites-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_bliss(self,**kwargs):
        defaults=dict(color='magenta', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(fileio.get_datadir(),'bliss-poly1.txt')
        self.draw_polygon(filename,**kwargs)

        filename = os.path.join(fileio.get_datadir(),'bliss-poly2.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_des(self,**kwargs):
        """ Draw the DES footprint on this Basemap instance.
        """
        defaults=dict(color='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(fileio.get_datadir(),'round13-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_smash(self,**kwargs):
        defaults=dict(facecolor='none',color='k')
        setdefaults(kwargs,defaults)

        filename = os.path.join(fileio.get_datadir(),'smash_fields_final.txt')
        smash=np.genfromtxt(filename,dtype=[('ra',float),('dec',float)],usecols=[4,5])
        xy = self.proj(smash['ra'],smash['dec'])
        self.scatter(*xy,**kwargs)

    def draw_decals(self,**kwargs):
        defaults=dict(color='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = os.path.join(fileio.get_datadir(),'decals-perimeter.txt')
        decals = np.genfromtxt(filename,names=['poly','ra','dec'])
        poly1 = decals[decals['poly'] == 1]
        poly2 = decals[decals['poly'] == 2]
        #self.draw_polygon_radec(poly1['ra'],poly1['dec'],**kwargs)
        #self.draw_polygon_radec(poly2['ra'],poly2['dec'],**kwargs)
        self.scatter(*self.proj(poly1['ra'],poly1['dec']))
        self.scatter(*self.proj(poly2['ra'],poly2['dec']))
        

    def draw_airmass(self, observatory, airmass, npts=360, **kwargs):
        defaults = dict(color='green', lw=2)
        setdefaults(kwargs,defaults)

        altitude_radians = (0.5 * np.pi) - np.arccos(1. / airmass)
        ra_contour = np.zeros(npts)
        dec_contour = np.zeros(npts)
        for ii, azimuth in enumerate(np.linspace(0., 2. * np.pi, npts)):
            ra_radians, dec_radians = observatory.radec_of(azimuth, '%.2f'%(np.degrees(altitude_radians)))
            ra_contour[ii] = np.degrees(ra_radians)
            dec_contour[ii] = np.degrees(dec_radians)
        xy = self.proj(ra_contour, dec_contour)
        self.plot(*xy, **kwargs)

        self.draw_zenith(observatory,**kwargs)

    def draw_zenith(self, observatory,**kwargs):
        """
        Plot a to-scale representation of the focal plane size at the zenith.
        """
        defaults = dict(color='green',alpha=0.75,lw=1.5)
        setdefaults(kwargs,defaults)

        # RA and Dec of zenith
        ra_zenith, dec_zenith = np.degrees(observatory.radec_of(0, '90'))
        xy = self.proj(ra_zenith, dec_zenith)

        self.plot(*xy,marker='+',ms=10,mew=1.5, **kwargs)
        self.tissot(ra_zenith, dec_zenith, constants.DECAM, 100, fc='none',**kwargs)

    def draw_moon(self, date):
        moon = ephem.Moon()
        moon.compute(date)
        ra_moon = np.degrees(moon.ra)
        dec_moon = np.degrees(moon.dec)
     
        proj = self.proj(np.array([ra_moon]), np.array([dec_moon]))
     
        if np.isnan(proj[0]).all() or np.isnan(proj[1]).all(): return
     
        self.scatter(*proj, color='%.2f'%(0.01 * moon.phase), edgecolor='black', s=500)
        color = 'black' if moon.phase > 50. else 'white'
        plt.text(proj[0], proj[1], '%.2f'%(0.01 * moon.phase),
                 fontsize=10, ha='center', va='center', color=color)

    def draw_jethwa(self,filename=None,log=True,**kwargs):
        import healpy as hp
        if not filename: 
            datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
            filename = os.path.join(datadir,'jethwa_satellites_n256.fits.gz')
        hpxmap = hp.read_map(filename)
        if log:
            return self.draw_hpxmap(np.log10(hpxmap),**kwargs)
        else:
            return self.draw_hpxmap(hpxmap,**kwargs)

    def draw_planet9(self,**kwargs):
        from scipy.interpolate import interp1d
        defaults=dict(color='b',lw=3)
        setdefaults(kwargs,defaults)
        datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
        ra_lo,dec_lo=np.genfromtxt(datadir+'p9_lo.txt',usecols=(0,1)).T
        ra_lo,dec_lo = self.roll(ra_lo,dec_lo)
        ra_hi,dec_hi=np.genfromtxt(datadir+'p9_hi.txt',usecols=(0,1)).T
        ra_hi,dec_hi = self.roll(ra_hi,dec_hi)
        self.plot(ra_lo,dec_lo,latlon=True,**kwargs)
        self.plot(ra_hi,dec_hi,latlon=True,**kwargs)

        orb = fileio.csv2rec(datadir+'P9_orbit_Cassini.csv')[::7]
        kwargs = dict(marker='o',s=40,edgecolor='none',cmap='jet_r')
        self.scatter(*self.proj(orb['ra'],orb['dec']),c=orb['cassini'],**kwargs)

    def draw_ligo(self,filename=None, log=True,**kwargs):
        import healpy as hp
        from astropy.io import fits as pyfits
        if not filename:
            datadir = '/home/s1/kadrlica/projects/bliss/v0/data/'
            filename = datadir + 'obsbias_heatmap_semesterA.fits'
        hpxmap = pyfits.open(filename)[0].data
        if log: self.draw_hpxmap(np.log10(hpxmap))
        else:   self.draw_hpxmap(hpxmap)
        
    def draw_lmc(self):
        proj = self.proj(RA_LMC,DEC_LMC)
        self.tissot(RA_LMC,DEC_LMC,RADIUS_LMC,100,fc='0.7',ec='0.5')
        plt.text(proj[0],proj[1], 'LMC', weight='bold',
                 fontsize=10, ha='center', va='center', color='k')

    def draw_smc(self):
        proj = self.proj(RA_SMC,DEC_SMC)
        self.tissot(RA_SMC,DEC_SMC,RADIUS_SMC,100,fc='0.7',ec='0.5')
        plt.text(proj[0],proj[1], 'SMC', weight='bold',
                 fontsize=8, ha='center', va='center', color='k')

    def draw_fields(self,fields,**kwargs):
        defaults = dict(edgecolor='k',facecolor='none',s=15)
        if self.projection == 'ortho': defaults.update(s=50)
        else: defaults.update(s=15)
        setdefaults(kwargs,defaults)
        self.scatter(*self.proj(fields['RA'],fields['DEC']),**kwargs)

    def draw_hpxmap(self, hpxmap, xsize=800, **kwargs):
        """
        Use pcolormesh to draw healpix map
        """
        import healpy
        if not isinstance(hpxmap,np.ma.MaskedArray):
            mask = ~np.isfinite(hpxmap) | (hpxmap==healpy.UNSEEN)
            hpxmap = np.ma.MaskedArray(hpxmap,mask=mask)

        vmin,vmax = np.percentile(hpxmap.compressed(),[0.1,99.9])

        defaults = dict(latlon=True, rasterized=True, vmin=vmin, vmax=vmax)
        setdefaults(kwargs,defaults)

        ax = plt.gca()

        lon = np.linspace(0, 360., xsize)
        lat = np.linspace(-90., 90., xsize)
        lon, lat = np.meshgrid(lon, lat)

        nside = healpy.get_nside(hpxmap.data)
        try:
            pix = healpy.ang2pix(nside,lon,lat,lonlat=True)
        except TypeError:
            pix = healpy.ang2pix(nside,np.radians(90-lat),np.radians(lon))

        values = hpxmap[pix]
        #mask = ((values == healpy.UNSEEN) | (~np.isfinite(values)))
        #values = np.ma.array(values,mask=mask)
        if self.projection is 'ortho':
            im = self.pcolor(lon,lat,values,**kwargs)
        else:
            im = self.pcolormesh(lon,lat,values,**kwargs)

        return im

    def draw_focal_planes(self, ra, dec, **kwargs):
        defaults = dict(alpha=0.2,color='red',edgecolors='none',lw=0)
        setdefaults(kwargs,defaults)
        ra,dec = np.atleast_1d(ra,dec)
        if len(ra) != len(dec):
            msg = "Dimensions of 'ra' and 'dec' do not match"
            raise ValueError(msg)
        decam = DECamFocalPlane()
        # Should make sure axis exists....
        ax = plt.gca()
        for _ra,_dec in zip(ra,dec):
            corners = decam.project(self,_ra,_dec)
            collection = matplotlib.collections.PolyCollection(corners,**kwargs)
            ax.add_collection(collection)
        plt.draw()

class DECamMcBride(DECamBasemap):
    def __init__(self,*args,**kwargs):
        defaults = dict(projection='mbtfpq',lon_0=0,rsphere=1.0,celestial=True)
        setdefaults(kwargs,defaults)
        super(DECamMcBride,self).__init__(*args, **kwargs)

class DECamOrtho(DECamBasemap):
    def __init__(self,*args,**kwargs):
        defaults = dict(projection='ortho',celestial=True,lon_0=0,lat_0=-30.17,rsphere=1.0)
        setdefaults(kwargs,defaults)
        super(DECamOrtho,self).__init__(*args, **kwargs)

class DECamFocalPlane(object):
    """Class for storing and manipulating the corners of the DECam CCDs.
    """

    filename = os.path.join(fileio.get_datadir(),'ccd_corners_xy_fill.dat')

    def __init__(self):
        # This is not safe. Use yaml instead (extra dependency)
        self.ccd_dict = eval(''.join(open(self.filename).readlines()))

        # These are x,y coordinates
        self.corners = np.array(self.ccd_dict.values())

        # Since we don't know the original projection of the DECam
        # focal plane into x,y it is probably not worth trying to
        # deproject it right now...

        #x,y = self.ccd_array[:,:,0],self.ccd_array[:,:,1]
        #ra,dec = Projector(0,0).image2sphere(x.flat,y.flat)
        #self.corners[:,:,0] = ra.reshape(x.shape)
        #self.corners[:,:,1] = dec.reshape(y.shape)

    def rotate(self, ra, dec):
        """Rotate the corners of the DECam CCDs to a given sky location.

        Parameters:
        -----------
        ra      : The right ascension (deg) of the focal plane center
        dec     : The declination (deg) of the focal plane center

        Returns:
        --------
        corners : The rotated corner locations of the CCDs
        """
        corners = np.copy(self.corners)

        R = SphericalRotator(ra,dec)
        _ra,_dec = R.rotate(corners[:,:,0].flat,corners[:,:,1].flat,invert=True)

        corners[:,:,0] = _ra.reshape(corners.shape[:2])
        corners[:,:,1] = _dec.reshape(corners.shape[:2])
        return corners

    def project(self, basemap, ra, dec):
        """Apply the given basemap projection to the DECam focal plane at a
        location given by ra,dec.

        Parameters:
        -----------
        basemap : The DECamBasemap to project to.
        ra      : The right ascension (deg) of the focal plane center
        dec     : The declination (deg) of the focal plane center

        Returns:
        --------
        corners : Projected corner locations of the CCDs
        """
        corners = self.rotate(ra,dec)

        x,y = basemap.proj(corners[:,:,0],corners[:,:,1])

        # Remove CCDs that cross the map boundary
        x[(np.ptp(x,axis=1) > np.pi)] = np.nan

        corners[:,:,0] = x
        corners[:,:,1] = y
        return corners

############################################################
# Depricated module functions

def drawDES(basemap, color='red'):
    msg = "drawDES is depricated; use DECamBasemap.draw_des instead."
    warnings.warn(msg)
    basemap.draw_des(color=color)

def drawSMASH(basemap, color='none', edgecolor='black', marker='h', s=50):
    msg = "drawSMASH is depricated; use DECamBasemap.draw_smash instead."
    warnings.warn(msg)
    basemap.draw_smash(color=color, edgecolor=edgecolor, marker=marker, s=s)

def drawMAGLITES(basemap, color='blue'):
    msg = "drawMAGLITES is depricated; use DECamBasemap.draw_smash instead."
    warnings.warn(msg)
    basemap.draw_maglites(color=color)

def drawAirmassContour(basemap, observatory, airmass, n=360, s=50):
    msg = "drawAirmassContour is depricated; use DECamBasemap.draw_airmass instead."
    warnings.warn(msg)
    basemap.draw_airmass(observatory=observatory, airmass=airmass, n=n, s=s)

def drawZenith(basemap, observatory):
    """
    Plot a to-scale representation of the focal plane size at the zenith.
    """
    msg = "drawZenith is depricated; use DECamBasemap.draw_zenith instead."
    warnings.warn(msg)
    basemap.draw_zenith(observatory)

def drawMoon(basemap, date):
    msg = "drawMoon is depricated; use DECamBasemap.draw_moon instead."
    warnings.warn(msg)
    basemap.draw_moon(date)

############################################################

def makePlot(date=None, name=None, figsize=(10.5,8.5), dpi=80, s=50, center=None, airmass=True, moon=True, des=True, smash=True, maglites=True, bliss=True, galaxy=True):
    """
    Create map in orthographic projection
    """
    if date is None: date = ephem.now()
    if type(date) != ephem.Date:
        date = ephem.Date(date)

    observatory = CTIO()
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

    if des:      basemap.draw_des()
    if smash:    basemap.draw_smash(s=s)
    if maglites: basemap.draw_maglites()
    #if bliss:    basemap.draw_bliss()
    if airmass:
        airmass = 2.0 if isinstance(airmass,bool) else airmass
        basemap.draw_airmass(observatory, airmass)
    if moon:     basemap.draw_moon(date)
    if galaxy:   basemap.draw_galaxy()

    plt.title('%s UTC'%(datestring(date)))

    return fig, basemap

def plotField(field, target_fields=None, completed_fields=None, options_basemap={}, **kwargs):
    """
    Plot a specific target field.
    """
    defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    #defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='gray_r')
    setdefaults(kwargs,defaults)

    if isinstance(field,np.core.records.record):
        tmp = FieldArray(1)
        tmp[0] = field
        field = tmp

    msg="%s: id=%10s, "%(datestring(field['DATE'][0],0),field['ID'][0])
    msg +="ra=%(RA)-6.2f, dec=%(DEC)-6.2f, secz=%(AIRMASS)-4.2f"%field[0]
    #msg = "%s: date=%s, "%(field['ID'][0],datestring(field['DATE'][0],0))
    #msg +="ra=%(RA)-6.2f, dec=%(DEC)-6.2f, secz=%(AIRMASS)-4.2f"%field[0]
    logging.info(msg)

    #if plt.get_fignums(): plt.cla()

    defaults = dict(date=field['DATE'][0], name='ortho')
    options_basemap = dict(options_basemap)
    setdefaults(options_basemap,defaults)
    fig, basemap = makePlot(**options_basemap)

    # Plot target fields
    if target_fields is not None:
        proj = basemap.proj(target_fields['RA'], target_fields['DEC'])
        basemap.scatter(*proj, c='0.8', **kwargs)

    # Plot completed fields
    if completed_fields is not None:
        proj = basemap.proj(completed_fields['RA'],completed_fields['DEC'])
        basemap.scatter(*proj, c=completed_fields['TILING'], **kwargs)

    # Try to draw the colorbar
    try:
        if len(fig.axes) == 2:
            # Draw colorbar in existing axis
            colorbar = plt.colorbar(cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar()
        colorbar.set_label('Tiling')
    except TypeError:
        pass
    plt.sca(fig.axes[0])

    # Show the selected field
    proj = basemap.proj(field['RA'], field['DEC'])
    basemap.scatter(*proj, c=COLORS[field[0]['FILTER']],edgecolor='none',s=50)

    return basemap

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
        basemap = plotField(fields[i],target_fields,completed_fields,options_basemap,**kwargs)
        #plt.savefig('field_%08i.png'%i)
        if completed_fields is None: completed_fields = FieldArray()
        completed_fields = completed_fields + fields[[i]]
        plt.pause(0.001)
        #time.sleep(0.01)
    return basemap

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

    plt.ioff()
    for i,f in enumerate(fields):
        plotField(fields[i],target_fields,completed_fields,**kwargs)
        png = os.path.join(tmpdir,'field_%08i.png'%i)
        plt.savefig(png,bbox_inches='tight',dpi=80)
        if completed_fields is None: completed_fields = FieldArray()
        completed_fields = completed_fields + fields[[i]]
    plt.ion()

    cmd = 'convert -delay 10 -loop 0 %s/*.png %s'%(tmpdir,outfile)
    logging.info(cmd)
    subprocess.call(cmd,shell=True)
    shutil.rmtree(tmpdir)
    return outfile

def plotWeights(date, target_fields, weights,options_basemap={},**kwargs):
    defaults = dict(c=weights, edgecolor='none', s=50, vmin=np.min(weights), vmax=np.min(weights) + 300., cmap='Spectral')
    setdefaults(kwargs,defaults)

    defaults = dict(date=date, name='ortho')
    options_basemap = dict(options_basemap)
    setdefaults(options_basemap,defaults)
    fig, basemap = makePlot(**options_basemap)

    proj = basemap.proj(target_fields['RA'], target_fields['DEC'])
    basemap.scatter(*proj, **kwargs)

    # Try to draw the colorbar
    try:
        if len(fig.axes) == 2:
            # Draw colorbar in existing axis
            colorbar = plt.colorbar(cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar()
        colorbar.set_label('Tiling')
    except TypeError:
        pass
    plt.sca(fig.axes[0])

def plotWeight(field, target_fields, weight, **kwargs):
    if isinstance(field,FieldArray):
        field = field[-1]

    date = ephem.Date(field['DATE'])

    if plt.get_fignums(): plt.cla()
    fig, basemap = obztak.utils.ortho.makePlot(date,name='weight')

    index_sort = np.argsort(weight)[::-1]
    proj = basemap.proj(target_fields['RA'][index_sort], target_fields['DEC'][index_sort])
    weight_min = np.min(weight)
    basemap.scatter(*proj, c=weight[index_sort], edgecolor='none', s=50, vmin=weight_min, vmax=weight_min + 300., cmap='Spectral')

    #cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
    #proj = obztak.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
    #basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)

    """
    cut_accomplished = np.in1d(self.target_fields['ID'],self.accomplished_fields['ID'])
    proj = obztak.utils.ortho.safeProj(basemap,
                                         self.target_fields['RA'][~cut_accomplished],
                                         self.target_fields['DEC'][~cut_accomplished])
    basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')

    proj = obztak.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
    basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    """

    # Draw colorbar in existing axis
    if len(fig.axes) == 2:
        colorbar = plt.colorbar(cax=fig.axes[-1])
    else:
        colorbar = plt.colorbar()
    colorbar.set_label('Weight')

    # Show the selected field
    proj = basemap.proj([field['RA']], [field['DEC']])
    basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

    #plt.draw()
    plt.pause(0.001)
    #fig.canvas.draw()

############################################################

def plot_progress(outfile=None,**kwargs):
    defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    for k,v in defaults.items():
        kwargs.setdefault(k,v)

    fields = FieldArray.load_database()

    nites = [get_nite(date) for date in fields['DATE']]
    nite = ephem.Date(np.max(nites))
    date = '%d/%02d/%d 00:00:00'%(nite.tuple()[:3])

    fig,basemap = makePlot(date=date,moon=False,airmass=False,center=(0,-90),smash=False)
    proj = basemap.proj(fields['RA'],fields['DEC'])
    basemap.scatter(*proj, c=fields['TILING'],  **kwargs)
    colorbar = plt.colorbar()
    colorbar.set_label('Tiling')
    plt.title('Maglites Coverage (%d/%02d/%d)'%nite.tuple()[:3])

    if outfile is not None:
        plt.savefig(outfile,bbox_inches='tight')

    return fig,basemap

def plot_bliss_coverage(fields,outfile=None,**kwargs):
    BANDS = ['g','r','i','z']
    filename = os.path.join(fileio.get_datadir(),'bliss-target-fields.csv')
    target = FieldArray.read(filename)
    target = target[~np.in1d(target.unique_id,fields.unique_id)]

    fig,ax = plt.subplots(2,2,figsize=(16,9))
    plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,bottom=0.01,top=0.99)
    kwargs = dict(edgecolor='none', s=12, alpha=0.2)

    for i,b in enumerate(BANDS):
        plt.sca(ax.flat[i])

        f = fields[fields['FILTER'] == b]
        t = target[target['FILTER'] == b]

        bmap = DECamMcBride()
        bmap.draw_des()
        bmap.draw_galaxy(10)

        proj = bmap.proj(t['RA'],t['DEC'])
        bmap.scatter(*proj, c='0.7', **kwargs)

        proj = bmap.proj(f['RA'],f['DEC'])
        bmap.scatter(*proj, c=COLORS[b], **kwargs)
        plt.gca().set_title('BLISS %s-band'%b)

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

