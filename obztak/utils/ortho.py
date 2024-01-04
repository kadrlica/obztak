import os
from os.path import expandvars
import shutil
import time
import logging
import tempfile
import subprocess
import warnings
from collections import OrderedDict as odict

from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.cm
import pylab as plt
import numpy as np
import healpy as hp
import ephem

from obztak import get_survey
import obztak.utils.projector
from obztak.utils.projector import cel2gal, gal2cel, SphericalRotator
from obztak.utils.projector import pix2ang

from obztak.utils import constants
from obztak.utils import fileio
from obztak.field import FieldArray
from obztak.utils.date import datestring,nite2utc,utc2nite,get_nite,setdefaults,datestr
from obztak.ctio import CTIO
from obztak.utils.constants import RA_LMC,DEC_LMC,RADIUS_LMC
from obztak.utils.constants import RA_SMC,DEC_SMC,RADIUS_SMC
from obztak.utils.constants import COLORS, CMAPS, TCOLORS
from obztak.utils.constants import FIGSIZE, SCALE, DPI

# ADW: This is bad...
#plt.ion()

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

DPI = 80

############################################################

class DECamBasemap(Basemap):

    def __init__(self, *args, **kwargs):
        super(DECamBasemap,self).__init__(*args,**kwargs)
        self.draw_parallels()
        self.draw_meridians()

    def draw_parallels(self,*args,**kwargs):
        defaults = dict()
        if not args: defaults.update(circles=np.arange(-90,120,30))
        setdefaults(kwargs,defaults)
        self.pardict = self.drawparallels(*args, **kwargs)
        return self.pardict

    def draw_meridians(self,*args,**kwargs):
        defaults = dict(labels=[1,0,0,1])
        if self.projection in ['ortho','geos','nsper','aeqd']:
            defaults.update(labels=[0,0,0,0])
        if not args: defaults.update(meridians=np.arange(0,420,60))
        setdefaults(kwargs,defaults)
        self.merdict = self.drawmeridians(*args,**kwargs)
        return self.merdict
        
    def proj(self,lon,lat):
        """ Remove points outside of projection """
        x, y = self(np.atleast_1d(lon),np.atleast_1d(lat))
        x[x > 1e29] = None
        y[y > 1e29] = None
        #return np.ma.array(x,mask=x>1e2),np.ma.array(y,mask=y>1e2)
        return x, y

    @staticmethod
    def wrap_index(lon, lat, wrap=180.):
        """ Find the index where the array wraps.
        """
        # No wrap: ignore
        if wrap is None:  return None

        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)

        # No array: ignore
        if len(lon)==1 or len(lat)==1: return None

        # Map [0,360)
        lon = np.mod(lon,360)
        wrap = np.mod(wrap,360)

        # Find the index of the entry closest to the wrap angle
        idx = np.abs(lon - wrap).argmin()
        # First or last index: ignore
        if idx == 0 or idx+1 == len(lon): return None
        # Value exactly equals wrap, choose next value
        elif (lon[idx] == wrap): idx += 1
        # Wrap angle sandwiched
        elif (lon[idx]<wrap) and (lon[idx+1]>wrap): idx += 1
        elif (lon[idx]<wrap) and (lon[idx-1]>wrap): idx += 0
        elif (lon[idx]>wrap) and (lon[idx+1]<wrap): idx += 1
        elif (lon[idx]>wrap) and (lon[idx-1]<wrap): idx += 0
        # There is no wrap: ignore
        else: return None

        return idx

    @classmethod
    def roll(cls,lon,lat,wrap=180.):
        """ Roll an lon,lat combination to split 180 boundary
        Parameters:
        -----------
        lon : right ascension (deg)
        lat: declination (deg)
        wrap_angle : angle to wrap at (deg)
        """
        lon = np.atleast_1d(lon)
        lat = np.atleast_1d(lat)

        # Do nothing
        if wrap is None: return lon,lat
        if len(lon)==1 or len(lat)==1: return lon,lat

        idx = cls.wrap_index(lon,lat,wrap)
        if idx is None: return lon, lat

        return np.roll(lon,-idx), np.roll(lat,-idx)


    @staticmethod
    def split(ra,angle=180):
        pass

    def path_select(self,filename,nside=512):
        npix = hp.nside2npix(nside)
        sel = np.zeros(npix,dtype=bool)
        pix = np.arange(npix)
        radec = np.array(pix2ang(nside,pix)).T

        data = np.genfromtxt(filename,names=['ra','dec','poly'])
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            path = Path(zip(poly['ra'],poly['dec']))
            sel |= path.contains_points(radec)

        return sel

    def path_area(self,filename,nside=512):
        sel = self.path_select(filename,nside)
        return sel.sum() * hp.nside2pixarea(nside,degrees=True)

    def draw_polygons(self,filename,**kwargs):
        """ Draw a polygon footprint on this Basemap instance.
        """
        defaults=dict(color='k', lw=2)
        setdefaults(kwargs,defaults)

        data = np.genfromtxt(filename,names=['ra','dec','poly'])
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            self.draw_polygon_radec(poly['ra'],poly['dec'],**kwargs)

    def draw_polygon(self,filename,**kwargs):
        """ Draw a polygon footprint on this Basemap instance.
        """
        defaults=dict(color='k', lw=2)
        setdefaults(kwargs,defaults)

        perim = np.loadtxt(filename,dtype=[('ra',float),('dec',float)])
        self.draw_polygon_radec(perim['ra'],perim['dec'],**kwargs)

    def draw_polygon_radec(self,ra,dec,**kwargs):
        ra,dec = self.roll(ra,dec)
        xy = self.proj(ra,dec)
        self.plot(*xy,**kwargs)
        
    def draw_galaxy(self,width=10,**kwargs):
        defaults = dict(color='k',lw=1.5,ls='-')
        setdefaults(kwargs,defaults)

        glon = np.linspace(0,360,200)
        glat = np.zeros_like(glon)
        ra,dec = self.roll(*gal2cel(glon,glat))

        #self.plot(*gal2cel(0,0),marker='o',ms=25,color=kwargs['color'],latlon=True)
        self.draw_polygon_radec(ra,dec,**kwargs)
        
        if width:
            kwargs.update(dict(ls='--',lw=1))
            for delta in [+width,-width]:
                ra,dec = self.roll(*gal2cel(glon,glat+delta))
                self.draw_polygon_radec(ra,dec,**kwargs)
            
    def draw_magellanic_stream(self,**kwargs):
        import fitsio
        defaults = dict(xsize=800, vmin=17., vmax=25.0, rasterized=True,
                        cmap=plt.cm.binary)
        setdefaults(kwargs,defaults)

        filename = get_datafile('allms_coldens_gal_nside_1024.fits')
        galhpx = fitsio.read(filename)['coldens']
        celhpx = obztak.utils.projector.hpx_gal2cel(galhpx)
        return self.draw_hpxmap(celhpx,**kwargs)

    def draw_sfd(self,**kwargs):
        import healpy as hp
        defaults = dict(rasterized=True,cmap=plt.cm.binary)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('lambda_sfd_ebv.fits')

        galhpx = hp.read_map(filename)
        celhpx = obztak.utils.projector.hpx_gal2cel(galhpx)
        return self.draw_hpxmap(np.log10(celhpx),**kwargs)

    def draw_maglites(self,**kwargs):
        defaults=dict(color='blue', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('maglites-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_bliss(self,**kwargs):
        defaults=dict(color='magenta', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('bliss-poly.txt')
        data = np.genfromtxt(filename,names=['ra','dec','poly'])
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            self.draw_polygon_radec(poly['ra'],poly['dec'],**kwargs)

    def draw_blissII(self,**kwargs):
        defaults=dict(color='darkorange', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('blissII-poly.txt')
        self.draw_polygons(filename,**kwargs)

    def draw_des(self,**kwargs):
        """ Draw the DES footprint on this Basemap instance.
        """
        defaults=dict(color='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('des-round19-poly.txt')
        self.draw_polygon(filename,**kwargs)

    def draw_smash(self,**kwargs):
        defaults=dict(facecolor='none',color='k')
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('smash_fields_final.txt')
        smash=np.genfromtxt(filename,dtype=[('ra',float),('dec',float)],usecols=[4,5])
        xy = self.proj(smash['ra'],smash['dec'])
        self.scatter(*xy,**kwargs)

    def draw_decals(self,**kwargs):
        defaults=dict(color='red', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('decals-perimeter.txt')
        decals = np.genfromtxt(filename,names=['poly','ra','dec'])
        poly1 = decals[decals['poly'] == 1]
        poly2 = decals[decals['poly'] == 2]
        #self.draw_polygon_radec(poly1['ra'],poly1['dec'],**kwargs)
        #self.draw_polygon_radec(poly2['ra'],poly2['dec'],**kwargs)
        self.scatter(*self.proj(poly1['ra'],poly1['dec']))
        self.scatter(*self.proj(poly2['ra'],poly2['dec']))

    def draw_delve(self,**kwargs):
        defaults=dict(color='blue', lw=2)
        setdefaults(kwargs,defaults)

        deep = odict([
            ('SextansB', (150.00,   5.33, 2.5)),
            ('IC5152',   (330.67, -51.30, 3.5)),
            ('NGC300',   ( 13.72, -37.68, 3.5)),
            ('NGC55',    (  3.79, -39.22, 3.5)),
            #('Peg4',     ( 328.5,   26.5, 2.0)),
            #('LMi',     ( 164.266, 28.877, 2.0)),
        ])
        for ra,dec,radius in deep.values():
            # This doesn't deal with boundaries well
            #self.tissot(ra, dec, radius, 100, fc='none',**kwargs)
            x,y = self.proj(np.array([ra]), np.array([dec]))
            self.scatter(x,y,facecolor='none',edgecolor=kwargs['color'],s=400)

        filename = fileio.get_datafile('delve-mc.txt')
        self.draw_polygon(filename,**kwargs)

        # DELVE-z
        #self.draw_polygon_radec([240,240,270,260,240],
        #                        [-20,-5 ,-5 ,-20,-20],
        #                        c='b',lw=2)
        #self.draw_polygon_radec([280,290,310,340,340,310,280],
        #                        [-30,-15,-15,-15,-30,-30,-30],
        #                        c='b',lw=2)

        # EDFS
        #x,y = self.proj(np.array([61.24]), np.array([-48.42]))
        #self.scatter(x,y,facecolor='none',edgecolor='b',s=600)

        #self.tissot(RA_LMC,DEC_LMC,25,100,fc='none',**kwargs)
        #self.tissot(RA_SMC,DEC_SMC,10,100,fc='none',**kwargs)

    def draw_magic(self,**kwargs):
        """ Draw the MAGIC footprint on this Basemap instance.
        """
        defaults=dict(color='magenta', lw=2)
        setdefaults(kwargs,defaults)

        filename = fileio.get_datafile('magic-poly.txt')
        self.draw_polygon(filename,**kwargs)


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
     
        x,y = self.proj(np.array([ra_moon]), np.array([dec_moon]))
        if np.isnan(x).all() or np.isnan(y).all(): return
     
        self.scatter(x,y,color='%.2f'%(0.01*moon.phase),edgecolor='black',s=600,zorder=99)
        color = 'black' if moon.phase > 50. else 'white'
        #text = '%.2f'%(0.01 * moon.phase)
        text = '%2.0f%%'%(moon.phase)
        plt.text(x, y, text, fontsize=10, ha='center', va='center', color=color,zorder=99)

    def draw_jethwa(self,filename=None,log=True,**kwargs):
        import healpy as hp
        if not filename: 
            filename = fileio.get_datafile('jethwa_satellites_n256.fits.gz')
        hpxmap = hp.read_map(filename)
        if log:
            return self.draw_hpxmap(np.log10(hpxmap),**kwargs)
        else:
            return self.draw_hpxmap(hpxmap,**kwargs)

    def draw_planet9(self,**kwargs):
        from scipy.interpolate import interp1d
        from scipy.interpolate import UnivariateSpline
        defaults=dict(color='b',lw=2)
        setdefaults(kwargs,defaults)

        ra_lo,dec_lo=np.genfromtxt(fileio.get_datafile('p9_lo.txt'),usecols=(0,1)).T
        ra_lo,dec_lo = self.roll(ra_lo,dec_lo)
        ra_lo += -360*(ra_lo > 180)
        ra_lo,dec_lo = ra_lo[::-1],dec_lo[::-1]
        ra_hi,dec_hi=np.genfromtxt(fileio.get_datafile('p9_hi.txt'),usecols=(0,1)).T
        ra_hi,dec_hi = self.roll(ra_hi,dec_hi)
        ra_hi += -360*(ra_hi > 180)
        ra_hi,dec_hi = ra_hi[::-1],dec_hi[::-1]

        spl_lo = UnivariateSpline(ra_lo,dec_lo)
        ra_lo_smooth = np.linspace(ra_lo[0],ra_lo[-1],360)
        dec_lo_smooth = spl_lo(ra_lo_smooth)

        spl_hi = UnivariateSpline(ra_hi,dec_hi)
        ra_hi_smooth = np.linspace(ra_hi[0],ra_hi[-1],360)
        dec_hi_smooth = spl_hi(ra_hi_smooth)

        #self.plot(ra_lo_smooth,dec_lo_smooth,latlon=True,**kwargs)
        #self.plot(ra_hi_smooth,dec_hi_smooth,latlon=True,**kwargs)

        orb = fileio.csv2rec(fileio.get_datafile('P9_orbit_Cassini.csv'))[::7]
        #kwargs = dict(marker='o',s=40,edgecolor='none',cmap='jet_r')
        #self.scatter(*self.proj(orb['ra'],orb['dec']),c=orb['cassini'],**kwargs)

        ra,dec = self.roll(orb['ra'],orb['dec'])
        self.plot(ra,dec,latlon=True,**kwargs)

    def draw_ligo(self,filename=None, log=True,**kwargs):
        import healpy as hp
        from astropy.io import fits as pyfits
        if not filename:
            filename = fileio.get_datafile('obsbias_heatmap_semesterA.fits')
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
        defaults = dict(edgecolor='none',s=15)
        if self.projection == 'ortho': defaults.update(s=50)
        colors = [COLORS[b] for b in fields['FILTER']]
        defaults.update(c=colors)
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
        if self.projection == 'ortho':
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
        self.observatory = CTIO()
        defaults = dict(projection='ortho',celestial=True,rsphere=1.0,
                        lon_0=0,lat_0=self.observatory.get_lat())
        setdefaults(kwargs,defaults)

        if 'date' in kwargs:
            kwargs.update(lon_0=self.parse_date(kwargs.pop('date')))

        super(DECamOrtho,self).__init__(*args, **kwargs)

    def draw_meridians(self,*args,**kwargs):
        cardinal = kwargs.pop('cardinal',False)
        meridict = super(DECamOrtho,self).draw_meridians(*args,**kwargs)
        ax = plt.gca()
        for mer in meridict.keys():
            ax.annotate(r'$%i^{\circ}$'%mer,self.proj(mer,5),ha='center')
        if cardinal:
            ax.annotate('West',xy=(1.0,0.5),ha='left',xycoords='axes fraction')
            ax.annotate('East',xy=(0.0,0.5),ha='right',xycoords='axes fraction')
        return meridict


    def parse_date(self,date):
        date = ephem.Date(date) if date else ephem.now()
        self.observatory.date = date

        # RA and Dec of zenith
        lon_zen, lat_zen = np.degrees(self.observatory.radec_of(0,'90'))
        return -lon_zen


class DECamFocalPlane(object):
    """Class for storing and manipulating the corners of the DECam CCDs.
    """

    filename = fileio.get_datafile('ccd_corners_xy_fill.dat')

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
        with np.errstate(invalid='ignore'):
            x[(np.ptp(x,axis=1) > np.pi)] = np.nan

        corners[:,:,0] = x
        corners[:,:,1] = y
        return corners

############################################################

def makePlot(date=None, name=None, figsize=(10.5,8.5), dpi=80, s=50, center=None, airmass=True, moon=True,
             des=True, smash=False, maglites=None, bliss=None, delve=None, galaxy=True):
    """
    Create map in orthographic projection
    """
    if date is None: date = ephem.now()
    if type(date) != ephem.Date:
        date = ephem.Date(date)


    fig = plt.figure(name, figsize=figsize, dpi=dpi)
    plt.cla()

    proj_kwargs = dict()
    if center: proj_kwargs.update(lon_0=center[0], lat_0=center[1])
    smap = DECamOrtho(date=date, **proj_kwargs)
    observatory = smap.observatory

    survey = get_survey()
    if survey=='des' or des:
        smap.draw_des()
        if airmass is True: airmass = 1.4
    if survey=='smash' or smash:
        smap.draw_smash(s=s)
    if 'maglites' in survey or maglites:
        smap.draw_maglites()
        if airmass is True: airmass = 2.0
    if survey=='bliss' or bliss:
        smap.draw_bliss()
        if airmass is True: airmass = 1.4
    if survey=='delve' or delve:
        smap.draw_delve()
    if survey=='magic':
        smap.draw_magic()

    if airmass:
        airmass = 2.0 if isinstance(airmass,bool) else airmass
        smap.draw_airmass(observatory, airmass)
    if moon:     smap.draw_moon(date)
    if galaxy:   smap.draw_galaxy()

    plt.title('%s UTC'%(datestring(date)))

    return fig, smap

def plotField(field, target_fields=None, completed_fields=None, options_basemap={}, **kwargs):
    """
    Plot a specific target field.

    Parameters:
    -----------
    field            : The specific field of interest.
    target_fields    : The fields that will be observed
    completed_fields : The fields that have been observed
    options_basemap  : Keyword arguments to the basemap constructor
    kwargs           : Keyword arguments to the matplotlib.scatter function

    Returns:
    --------
    basemap : The basemap object
    """
    if isinstance(field,np.core.records.record):
        tmp = FieldArray(1)
        tmp[0] = field
        field = tmp
    band = field[0]['FILTER']
    cmap = matplotlib.cm.get_cmap(CMAPS[band])
    defaults = dict(marker='H',s=100,edgecolor='',vmin=-1,vmax=4,cmap=cmap)
    #defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    #defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='gray_r')
    defaults = dict(marker='H',s=100,edgecolor='none',vmin=-1,vmax=4,cmap=cmap)
    setdefaults(kwargs,defaults)

    msg="%s: id=%10s, "%(datestring(field['DATE'][0],0),field['ID'][0])
    msg +="ra=%(RA)-6.2f, dec=%(DEC)-6.2f, secz=%(AIRMASS)-4.2f"%field[0]
    logging.info(msg)

    defaults = dict(date=field['DATE'][0], name='ortho')
    options_basemap = dict(options_basemap)
    setdefaults(options_basemap,defaults)
    fig, basemap = makePlot(**options_basemap)
    plt.subplots_adjust(left=0.03,right=0.97,bottom=0.03,top=0.97)

    # Plot target fields
    if target_fields is not None and len(target_fields):
        sel = target_fields['FILTER']==band
        x,y = basemap.proj(target_fields['RA'], target_fields['DEC'])
        kw = dict(kwargs,c='w',edgecolor='0.6',s=0.8*kwargs['s'])
        basemap.scatter(x[sel], y[sel], **kw)
        kw = dict(kwargs,c='w',edgecolor='0.8',s=0.8*kwargs['s'])
        basemap.scatter(x[~sel], y[~sel], **kw)

    # Plot completed fields
    if completed_fields is not None and len(completed_fields):
        sel = completed_fields['FILTER']==band
        x,y = basemap.proj(completed_fields['RA'],completed_fields['DEC'])
        kw = dict(kwargs)
        basemap.scatter(x[~sel], y[~sel], facecolor='0.6', **kw)
        basemap.scatter(x[sel], y[sel], c=completed_fields['TILING'][sel], **kw)

    # Try to draw the colorbar
    try:
        if len(fig.axes) == 2:
            # Draw colorbar in existing axis
            colorbar = plt.colorbar(cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar()
        colorbar.set_label('Tiling (%s-band)'%band)
    except TypeError:
        pass
    plt.sca(fig.axes[0])

    # Show the selected field
    x,y = basemap.proj(field['RA'], field['DEC'])
    kw = dict(kwargs,edgecolor='k')
    basemap.scatter(x,y,c=COLORS[band],**kw)

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
        if completed_fields is None: completed_fields = FieldArray()
        completed_fields = completed_fields + fields[[i]]
        plt.pause(0.001)

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
        plt.savefig(png,dpi=DPI)
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

def plot_coverage(fields,nitestr,outfile=None,**kwargs):
    """ Plot the BLISS survey coverage

    Parameters:
    -----------
    fields  : the bliss fields to plot
    outfile : the output file to write to
    kwargs  : plotting keyword arguments

    Returns:
    --------
    None
    """
    from obztak import factory
    bands = ['g','r','i','z']
    defaults = dict(edgecolor='none', alpha=0.2, vmin=-1, vmax=2)
    setdefaults(kwargs,defaults)

    filename = factory.scheduler_factory()._defaults['targets']

    target = factory.field_factory().read(filename)
    target = target[~np.in1d(target.unique_id,fields.unique_id)]

    for pro in ['mbt','ort']:
        if pro == 'mbt':
            fig,ax = plt.subplots(2,2,figsize=(16,9))
            plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,
                                bottom=0.01,top=0.99)
            kwargs['s'] = 12 # roughly scaled to image
            def create_skymap(): return DECamMcBride()
        else:
            fig,ax = plt.subplots(2,2,figsize=(12,12))
            plt.subplots_adjust(wspace=0.01,hspace=0.05,left=0.01,right=0.99,
                                bottom=0.01,top=0.97)
            kwargs['s'] = 45 # scaled to image
            def create_skymap(): return DECamOrtho(date='2016/2/11 03:00',lon_0=0,lat_0=-90)

        for i,b in enumerate(bands):
            plt.sca(ax.flat[i])

            f = fields[fields['FILTER'] == b]
            t = target[target['FILTER'] == b]

            bmap = create_skymap()
            bmap.draw_des()
            bmap.draw_galaxy(10)

            proj = bmap.proj(t['RA'],t['DEC'])
            bmap.scatter(*proj, c='0.7', **kwargs)

            proj = bmap.proj(f['RA'],f['DEC'])
            bmap.scatter(*proj, c=f['TILING'], cmap=CMAPS[b], **kwargs)
            plt.gca().set_title('DECam %s-band'%b)

        plt.suptitle('Coverage (%s)'%nitestr,fontsize=16)
        plt.savefig('nightsum_summary_%s_%s.png'%(nitestr,pro))


def plot_completion(fields,tonight,nitestr,outfile=None,**kwargs):
    """ Plot the BLISS survey coverage

    Parameters:
    -----------
    fields  : all completed fields
    tonight : fields that were completed tonight
    outfile : the output file to write to
    kwargs  : plotting keyword arguments

    Returns:
    --------
    None
    """
    from obztak import factory
    bands = ['g','r','i','z']
    tiles = [1,2,3,4,9]
    DONE  = -1
    defaults = dict(alpha=1.0,rasterized=True,lw=0.75)
    setdefaults(kwargs,defaults)

    filename = factory.scheduler_factory()._defaults['targets']

    targets = factory.field_factory().read(filename)
    targets['PRIORITY'][np.in1d(targets.unique_id,fields.unique_id)] = DONE

    for pro in ['mbt','ort']:
        if pro == 'mbt':
            fig,ax = plt.subplots(2,2,figsize=(16,9))
            plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,
                                bottom=0.01,top=0.99)
            kwargs['s'] = 12 # roughly scaled to image
            def create_skymap(): return DECamMcBride()
        else:
            fig,ax = plt.subplots(2,2,figsize=(12,12))
            plt.subplots_adjust(wspace=0.01,hspace=0.05,left=0.01,right=0.99,
                                bottom=0.01,top=0.97)
            kwargs['s'] = 45 # scaled to image
            def create_skymap(): return DECamOrtho(date='2016/2/11 03:00',lon_0=0,lat_0=-90)

        for i,b in enumerate(bands):
            plt.sca(ax.flat[i])

            bmap = create_skymap()
            bmap.draw_des()
            bmap.draw_galaxy(10)
            bmap.draw_delve(color='gray')

            t = targets[targets['FILTER'] == b]
            n = tonight[tonight['FILTER'] == b]

            # Exposures done
            for tile in tiles:
                if tile == 1:
                    todo = (t['TILING'] == tile) & (t['PRIORITY'] >= 0)
                    proj = bmap.proj(t[todo]['RA'],t[todo]['DEC'])
                    bmap.scatter(*proj, c=TCOLORS[0], edgecolor='none', **kwargs)

                done = (t['TILING'] == tile) & (t['PRIORITY'] == DONE)
                proj = bmap.proj(t[done]['RA'],t[done]['DEC'])
                bmap.scatter(*proj, c=TCOLORS[tile], edgecolor='none', **kwargs)

            # Exposures done tonight
            for tile in tiles:
                sel = (n['TILING']==tile)
                proj = bmap.proj(n[sel]['RA'],n[sel]['DEC'])
                bmap.scatter(*proj, c=TCOLORS[tile], edgecolor='orangered', **kwargs)

            plt.gca().set_title('DECam %s-band'%b)

        # Plot legend
        plt.sca(ax.flat[0])
        for tile,color in TCOLORS.items():
            plt.plot(np.nan,np.nan,'o',mfc=color,mec='none',label="Tile %s"%tile)
        plt.plot(np.nan,np.nan,'o',mfc='none',mec='orangered',label='Tonight')
        plt.legend(loc='center',bbox_to_anchor=(1.0,0.0))

        plt.suptitle('Completion (%s)'%nitestr,fontsize=16)
        plt.savefig('nightsum_complete_%s_%s.png'%(nitestr,pro),dpi=200)

def plot_focal_planes(fields,nightstr):
    fig,axes = plt.subplots(1,2,figsize=(12,5))
    for i,d in enumerate(['2017/02/08 07:00:00','2017/02/08 19:00:00']):
        plt.sca(axes[i])
        bmap = DECamOrtho(date=d)
        for b in np.unique(fields['FILTER']):
            f = fields[fields['FILTER']==b]
            bmap.draw_focal_planes(f['RA'],f['DEC'],color=COLORS[b],alpha=0.3)
        bmap.draw_bliss()
        bmap.draw_galaxy()
        bmap.draw_des()

    plt.suptitle('Coverage (%s)'%nitestr,fontsize=16)
    plt.savefig('nightsum_summary_%s.png'%nitestr)

def plot_psf(new,old,nbins=35):
    kwargs = dict(normed=True)
    step_kwargs = dict(kwargs,histtype='step',lw=3.5)
    fill_kwargs = dict(kwargs,histtype='stepfilled',lw=1.0,alpha=0.7)
    step_kwargs['bins'] = np.linspace(0.5,2.5,nbins)
    fill_kwargs['bins'] = np.linspace(0.5,2.5,nbins)

    plt.hist(new['psf'],color='green',zorder=10, label='tonight', **fill_kwargs)
    plt.hist(new['psf'],color='green',zorder=10, **step_kwargs)
    plt.hist(old['psf'],color='0.5', label='previously', **fill_kwargs)
    plt.hist(old['psf'],color='0.5', **step_kwargs)
    plt.axvline(1.20,ls='--',lw=2,color='gray')
    plt.legend()
    plt.xlabel('FWHM (arcsec)')
    plt.ylabel('Normalized Number of Exposures')
    
def plot_teff(new,old,nbins=35):
    kwargs = dict(normed=True)
    step_kwargs = dict(kwargs,histtype='step',lw=3.5)
    fill_kwargs = dict(kwargs,histtype='stepfilled',lw=1.0,alpha=0.7)
    step_kwargs['bins'] = np.linspace(0,1.5,nbins)
    fill_kwargs['bins'] = np.linspace(0,1.5,nbins)

    plt.hist(new['teff'],color='green',zorder=10,label='tonight', **fill_kwargs)
    plt.hist(new['teff'],color='green',zorder=10, **step_kwargs)
    plt.hist(old['teff'],color='0.5',label='previously', **fill_kwargs)
    plt.hist(old['teff'],color='0.5', **step_kwargs)
    plt.axvline(0.25,ls='--',lw=2,color='gray')
    plt.legend()
    plt.xlabel('Teff')
    plt.ylabel('Normalized Number of Exposures')


def plot_nightsum(fields,nitestr,date):
    """ Plot the bliss night summary. 

    Parameters:
    -----------
    fields:  the fields observed tonight
    nitestr: the nite in strig format

    Returns:
    --------
    None
    """
    from obztak.utils.database import Database
    plt.ioff()


    # Select the fields from the database
    db = Database()
    db.connect()
    query = """
    select id, qc_fwhm as psf, qc_teff as teff, filter from exposure
    where delivered = True and propid = '%s'
    and flavor = 'object'
    and qc_teff is not NULL
    and qc_fwhm is not NULL
    and date %s
    """

    d = datestr(date)
    q = query%(fields.PROPID,"between (timestamp '%s') AND (timestamp '%s' + interval '12 hours')"%(d,d))
    logging.debug(q)
    new = db.query2recarray(q)
    try:
        q = query%(fields.PROPID,"< (timestamp '%s')"%d)
        logging.debug(q)
        old = db.query2recarray(q)
    except ValueError as e:
        print(e)
        old = np.recarray(0,dtype=new.dtype)

    for b in ['u','g','r','i','z','Y']:
        f = new[new['filter'] == b]
        print(' %s-band:'%b, len(f))

    if not len(new):
        logging.warn("No new exposures...")
        return

    new_sel = (np.array(map(utc2nite,fields['DATE'])) == nitestr)
    new_fields = fields[new_sel]
    old_fields = fields[~new_sel]

    ##########################
    #print("Plotting coverage...")
    #plot_coverage(fields,nitestr)

    ##########################
    #print("Plotting focal planes...")
    #plot_focal_planes(fields,nitestr)

    ##########################
    print("Plotting completion...")
    plot_completion(fields,new_fields,nitestr)

    ##########################

    fig,axes = plt.subplots(1,2,figsize=(12,5))
    plt.sca(axes[0])
    plt.plot(np.nan,np.nan,'-w',label='all')
    plot_psf(new,old)
    plt.title('Seeing (%s)'%nitestr)
    plt.sca(axes[1])
    plt.plot(np.nan,np.nan,'-w',label='all')
    plot_teff(new,old)
    plt.title('Effective Depth (%s)'%nitestr)
    plt.savefig('nightsum_psf_teff_%s.png'%nitestr,bbox_inches='tight')

    fig,axes = plt.subplots(2,2,figsize=(14,10))
    axes = axes.flatten()
    for i,b in enumerate(['g','r','i','z']):
        plt.sca(axes[i])
        plt.plot(np.nan,np.nan,'-w',label='%s-band'%b)
        plot_psf(new[new['filter'] == b],old[old['filter'] == b])
    plt.savefig('nightsum_psf_%s.png'%nitestr,bbox_inches='tight')

    fig,axes = plt.subplots(2,2,figsize=(14,10))
    axes = axes.flatten()
    for i,b in enumerate(['g','r','i','z']):
        plt.sca(axes[i])
        plt.plot(np.nan,np.nan,'-w',label='%s-band'%b)
        plot_teff(new[new['filter'] == b],old[old['filter'] == b])
    plt.savefig('nightsum_teff_%s.png'%nitestr,bbox_inches='tight')

############################################################

if __name__ == '__main__':
    makePlot('2016/2/10 03:00')

############################################################

