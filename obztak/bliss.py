#!/usr/bin/env python
"""
Code related to the Magellanic Satellites Survey (MagLiteS).
"""
import os,sys
import logging
import copy
from collections import OrderedDict as odict

import numpy as np
import fitsio

from obztak.field import FieldArray, SISPI_DICT, SEP
from obztak.survey import Survey
from obztak.scheduler import Scheduler

from obztak.utils.projector import cel2gal, angsep
from obztak.utils import constants
from obztak.utils import fileio
from obztak.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS,COLORS
from obztak.utils.constants import COLORS, CMAPS
from obztak.utils.date import datestring, setdefaults, nite2utc,utc2nite,datestr

NAME    = 'BLISS'
PROGRAM = NAME.lower()
PROPID  = '2017A-0260'
PROPOSER = 'Soares-Santos'
BANDS = ['g','r','i','z']
TILINGS = 4

class BlissSurvey(Survey):
    """ Survey sublcass for BLISS. """

    # 2017A SCHEDULED (bliss-windows.csv is actually used)
    nights = [
        ['2017/02/07', 'second'], # phase=90%, set=07:40 (1.5h dark)
        ['2017/02/08', 'second'], # phase=95%, set=08:40 (0.5h dark)
        ['2017/02/13', 'second'], # phase=90%, set=nope
        ['2017/02/14', 'second'], # phase=82%, set=nope
        ['2017/02/19', 'second'], # phase=37%, set=nope
        ['2017/02/20', 'second'], # phase=28%, rise=05:51 (1h dark)
        ['2017/03/07', 'first'],  # phase=77%, set=nope
        ['2017/03/17', 'second'], # phase=73%, set=nope
        ['2017/04/13', 'second'], # phase=92%, set=nope
        ['2017/04/14', 'second'], # phase=86%, set=nope
        ['2017/05/02', 'first'],  # phase=49%, set=04:20 (0.2h dark)
        ['2017/05/30', 'first'],  # phase=34%, set=03:12 (1.3h dark)
        ['2017/05/31', 'first'],  # phase=45%, set=04:11 (0.3h dark)
        ['2017/06/01', 'first'],  # phase=55%, set=nope
        ['2017/06/02', 'full'],   # phase=65%, set=06:04 (5h dark)
        ['2017/06/03', 'full'],   # phase=75%, set=06:58 (4h dark)
        ['2017/06/04', 'full'],   # phase=82%, set=07:52 (3h dark)
        ['2017/07/01', 'second'], # phase=62%, set=05:47 (5h dark)
        ['2017/07/02', 'second'], # phase=71%, set=06:40 (4h dark)
        ['2017/07/15', 'second'], # phase=57%, set=nope
        ]

    extra_night = [
        ['2017/08/05', 'second'], # phase=99%, set=06:40 (4h dark)
        ['2017/08/06', 'full'],   # phase=100%, set=nope
        ['2017/08/07', 'full'],   # phase=99%, set=nope
    ]

    alfredo_nights = [
        ['2017/03/06', 'full'],  # level=9
        ['2017/03/07','second'], # level=10
        ['2017/06/21','first'],
        ]
    
    def prepare_fields(self, infile=None, outfile=None, mode='bliss_rotate', plot=True, smcnod=False):
        """Create the list of fields to be targeted by the BLISS survey.

        Selection starts with 3 regions:
        - P9 - Planet 9 region above DES footprint (priority 1)
        - LIGO - Region targeted based on LIGO sensitivity maps (and good for MW)
        - Alfredo - Overlap with Alfredo's eRosita companion survey.

        Fields that have been previously covered by DECam are removed
        from the LIGO and Alfredo fooptrint regions.

        Parameters:
        -----------
        infile : File containing all possible field locations.
        outfile: Output file of selected fields
        mode   : Mode for dithering. default: 'bliss_rotate'
        plot   : Create an output plot of selected fields.

        Returns:
        --------
        fields : A FieldArray of the selected fields.
        """
        # Import the dither function here...
        #def dither(ra,dec,dx,dy):
        #    return ra,dec

        if mode is None or mode.lower() == 'none':
            def dither(ra,dec,dx,dy):
                return ra,dec
            OFFSETS = TILINGS*[(0,0)]
        elif mode.lower() == 'smash_dither':
            OFFSETS = [(0,0), (1.0,0.0),
                       (-1.0,0.0), (0.0,-0.75)][:TILINGS]
            dither = self.smash_dither
        elif mode.lower() == 'smash_rotate':
            OFFSETS = [(0,0), (0.75,0.75),
                       (-0.75,0.75), (0.0,-0.75)][:TILINGS]
            dither = self.smash_rotate
        elif mode.lower() == 'decam_dither':
            OFFSETS = [(0., 0.),(8/3.*CCD_X, -11/3.*CCD_Y),
                       (8/3.*CCD_X, 8/3.*CCD_Y),(-8/3.*CCD_X, 0.)][:TILINGS]
            dither = self.decam_dither
        elif mode.lower() == 'coord_rotate':
            OFFSETS = [(0,0), (153.0,-17.0),
                       (-1.0,1.0), (0.0,-1.0)][:TILINGS]
            dither = self.coord_rotate
        elif mode.lower() == 'decals_rotate':
            OFFSETS = [(0,0),(-0.2917, 0.0833),
                       (-0.5861, 0.1333)][:TILINGS]
            dither = self.decals_rotate
        elif mode.lower() == 'bliss_rotate':
            OFFSETS = [(0, 0), (8/3.*CCD_X, -11/3.*CCD_Y),
                       (8/3.*CCD_X, 8/3.*CCD_Y), (-8/3.*CCD_X, 0.)][:TILINGS]
            dither = self.decals_rotate
        else:
            msg = "Unrecognized dither mode: %s"%mode
            raise ValueError(msg)
        logging.info("Dither mode: %s"%mode.lower())

        if infile is None:
            #infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
            #infile = os.path.join(fileio.get_datadir(),'ctr-healpy-32-13131.txt')
            infile = os.path.join(fileio.get_datadir(),'decam-tiles_obstatus.fits')
        #data = np.recfromtxt(infile, names=True)
        raw_data = fitsio.read(infile)
        data = raw_data[(raw_data['PASS'] == 1)]

        # Apply footprint selection after tiling/dither
        #sel = obztak.utils.projector.footprint(data['RA'],data['DEC'])

        # This is currently a non-op
        decals_id = data['TILEID']
        ra        = data['RA']
        dec       = data['DEC']

        nhexes = len(data)
        #ntilings = len(DECAM_DITHERS)
        ntilings = TILINGS
        nbands = len(BANDS)
        nfields = nhexes*nbands*ntilings

        logging.info("Number of hexes: %d"%nhexes)
        logging.info("Number of tilings: %d"%ntilings)
        logging.info("Number of filters: %d"%nbands)

        fields = FieldArray(nfields)
        fields['HEX'] = np.tile(np.repeat(decals_id,nbands),ntilings)
        fields['PRIORITY'].fill(1)
        fields['TILING'] = np.repeat(np.arange(1,ntilings+1),nhexes*nbands)
        fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)

        #for i in range(ntilings):
        for i,offset in enumerate(OFFSETS):
            idx0 = i*nhexes*nbands
            idx1 = idx0+nhexes*nbands
            ra_dither,dec_dither = dither(ra,dec,offset[0],offset[1])
            #ra_dither = raw_data[raw_data['PASS'] == i+1]['RA']
            #dec_dither = raw_data[raw_data['PASS'] == i+1]['DEC']
            fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
            fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

        # Apply footprint selection after tiling/dither
        #sel = self.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION

        # Apply footprint selection after tiling/dither
        p9 = self.planet9(fields['RA'],fields['DEC'])
        ligo = self.ligo_mw(fields['RA'],fields['DEC'])
        alfredo = self.alfredo(fields['RA'],fields['DEC'])
        p9v2 = self.planet9v2(fields['RA'],fields['DEC'])

        fields['PRIORITY'][p9] = 1
        fields['PRIORITY'][ligo] = 2
        fields['PRIORITY'][alfredo] = 3

        #sel = (p9 | ligo | alfredo)
        sel = (p9 | ligo | alfredo | p9v2)

        # Apply telescope constraints
        sel &= (fields['DEC'] > constants.SOUTHERN_REACH)

        # Apply covered fields
        #sel &= self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]
        # Apply covered fields (but take everything in P9 region)
        uncovered = self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]
        sel &= ( (p9v2) | (p9 & (fields['RA'] > 180)) | uncovered )

        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))
        logging.debug("Unique priorities: ",np.unique(fields['PRIORITY']))

        if plot:
            import pylab as plt
            from obztak.utils.ortho import makePlot
            from obztak.utils.ortho import DECamOrtho, DECamMcBride

            kwargs = dict(edgecolor='none',cmap='viridis_r',vmin=0,vmax=ntilings)

            fig,ax = plt.subplots(2,2,figsize=(16,9))
            plt.subplots_adjust(wspace=0.01,hspace=0.02,left=0.01,right=0.99,
                                bottom=0.01,top=0.99)
            for i,b in enumerate(BANDS):
                plt.sca(ax.flat[i])
                bmap = DECamMcBride()
                bmap.draw_galaxy()
                bmap.draw_des()
                f = fields[fields['FILTER'] == b]
                bmap.scatter(*bmap.proj(f['RA'],f['DEC']),c=COLORS[b],s=15,**kwargs)

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'_mbt.png'
                plt.savefig(outfig,bbox_inches='tight')

            fig,ax = plt.subplots(2,2,figsize=(10,10))
            plt.subplots_adjust(wspace=0.05,hspace=0.05,left=0.05,right=0.95,
                                bottom=0.05,top=0.95)
            for i,b in enumerate(BANDS):
                plt.sca(ax.flat[i])
                bmap = DECamOrtho(date='2017/06/02 03:00')
                bmap.draw_galaxy()
                bmap.draw_des()
                f = fields[fields['FILTER'] == b]
                bmap.scatter(*bmap.proj(f['RA'],f['DEC']),c=COLORS[b],s=50,**kwargs)

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'_ortho.png'
                plt.savefig(outfig,bbox_inches='tight')

            if not sys.flags.interactive:
                plt.show(block=True)

        # Just 3rd tiling
        #fields = fields[fields['TILING'] == 3]

        if outfile: fields.write(outfile)

        return fields

    @staticmethod
    def footprint(ra,dec):
        l, b = cel2gal(ra, dec)
        # Starting to the north of the Galactic plane
        sel = ((ra > 120) & (ra < 360)) | (ra < 10)
        # Band in declination
        sel &= (dec < -30 ) & (dec > -40)
        # LIGO/DSPHS high probability region south of nominal stripe
        sel |= ((ra > 140) & (ra < 185) & (dec < -30) & (dec > -60))
        # Alfredo's eRosita survey
        sel |= BlissSurvey.alfredo(ra,dec)
        # 10 deg from the Galactic plane
        sel &= (np.fabs(b) > 10.)
        return sel

    @staticmethod
    def footprint(ra,dec):
        # The high-probability region for LIGO/MW dwarfs
        sel = BlissSurvey.ligo_mw(ra,dec)
        # Alfredo's eRosita survey
        sel |= BlissSurvey.alfredo(ra,dec)
        # Planet 9 Region
        sel |= BlissSurvey.planet9(ra,dec)

        return sel

    @staticmethod
    def ligo_mw(ra,dec):
        """The high-probability region for LIGO/MW dwarfs"""
        l, b = cel2gal(ra, dec)
        # Starting to the north of the Galactic plane
        sel = (ra > 120) & (ra < 270)
        # Upper limit in declination
        sel &= (dec < -30 )
        # 10 deg above the Galactic plane
        sel &= (b > 10.)
        return sel

    @staticmethod
    def planet9(ra,dec):
        """The high-probability region for Planet 9"""
        sel = ((ra > 305)|(ra < 15)) & (dec > -40) & (dec < -30) # v0, v4
        #sel = ((ra > 305)|(ra < 15)) & (dec > -40) & (dec < -25) # v1
        #sel = ((ra > 305)|(ra < 15)) & (dec > -35) & (dec < -25) # v2
        #sel = ((ra > 305)|(ra < 15)) & (dec > -40) & (dec < -28) # v3
        return sel

    @staticmethod
    def planet9v2(ra,dec):
        """The other high-probability region for Planet 9"""
        from matplotlib.path import Path
        data = np.genfromtxt(fileio.get_datafile('blissII-poly.txt'),
                             names=['RA','DEC','POLY'])
        poly = data[data['POLY'] == 1]
        path = Path(zip(poly['RA'],poly['DEC']))
        sel = path.contains_points(np.vstack([ra,dec]).T)

        poly = data[data['POLY'] == 4]
        path = Path(zip(poly['RA'],poly['DEC']))
        sel |= path.contains_points(np.vstack([ra,dec]).T)

        return sel

    @staticmethod
    def alfredo(ra,dec):
        """Alfredo's eRosita survey extension"""
        l, b = cel2gal(ra, dec)
        sel = ((b > 20) & (b < 30))
        sel &= ((l < 360) & (l > 240))
        #sel &= ((ra > 135) & (ra < 245))
        sel &= ((ra > 135) & (ra < 180))
        sel &= ((dec > -30) & (dec < -10))
        return sel

    @staticmethod
    def uncovered(ra,dec,band):
        """
        Determine which fields haven't been previously covered by DECam

        Parameters:
        -----------
        ra  : right ascension of field
        dec : declination of field
        band: observing band

        Returns:
        --------
        sel, frac : selection of fields and coverage fraction
        """
        import healpy as hp
        #dirname = '/home/s1/kadrlica/projects/bliss/v0/data'
        #basename = 'decam_coverage_90s_%s_n1024.fits.gz'
        dirname = '/Users/kadrlica/bliss/observing/data'
        basename1 = 'decam_max_expmap_%s_n1024.fits.gz'
        basename2 = 'decam_sum_expmap_%s_n1024.fits.gz'

        sel = np.ones_like(ra,dtype=bool)
        frac = np.zeros_like(ra,dtype=float)
        ra,dec,band=np.atleast_1d(ra),np.atleast_1d(dec),np.atleast_1d(band)

        for b in np.unique(band):
            idx = (band==b)
            filename1 = os.path.join(dirname,basename1%b)
            logging.info("Reading %s..."%os.path.basename(filename1))
            skymap1 = hp.read_map(filename1,verbose=False)

            filename2 = os.path.join(dirname,basename2%b)
            logging.info("Reading %s..."%os.path.basename(filename2))
            skymap2 = hp.read_map(filename2,verbose=False)

            # Apply minimum exposure time selection
            #skymap = (skymap > 30)
            #skymap = (skymap > 45)
            #skymap = (skymap > 60)
            #skymap = (skymap > 90)

            skymap = (skymap1 > 30) & (skymap2 > 45)

            nside = hp.get_nside(skymap)
            vec = hp.ang2vec(np.radians(90.-dec[idx]),np.radians(ra[idx]))
            f = []
            for i,v in enumerate(vec):
                print '\r%s/%s'%(i+1,len(vec)),
                sys.stdout.flush()
                pix = hp.query_disc(nside,v,np.radians(constants.DECAM))
                f.append(skymap[pix].sum()/float(len(pix)))
            print
            frac[idx] = np.array(f)
            
        sel = (frac < 2/3.)
        return sel,frac

class BlissFieldArray(FieldArray):
    """ Array of BLISS fields """
    PROGRAM  = PROGRAM
    PROPID   = PROPID
    PROPOSER = PROPOSER

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID
    SISPI_DICT["proposer"] = PROPOSER

    OBJECT_FMT = 'BLISS field' + SEP + ' %s'
    SEQID_FMT = 'BLISS scheduled' + SEP + ' %(DATE)s'
    BANDS = BANDS

    @classmethod
    def query(cls, **kwargs):
        """ Generate the database query.

        Parameters:
        -----------
        kwargs : Keyword arguments to fill the query.

        Returns:
        --------
        query  : The query string.
        """
        defaults = dict(propid=cls.SISPI_DICT['propid'], limit='',
                        object_fmt = cls.OBJECT_FMT%'')
        kwargs = setdefaults(kwargs,copy.deepcopy(defaults))

        # Should pull this out to be accessible (self.query())?
        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW
        FROM exposure where propid = '%(propid)s' and exptime > 89
        and discard = False and delivered = True and flavor = 'object'
        and object like '%(object_fmt)s%%'
        -- i-band AOS failures 'sqrt(pow(qc_fwhm,2)-pow(dimm2see,2)) > 0.9'
        and id not in (652771,652794,652795,652796,652797,652799,652800,652803,652804,652806,652807,652808,652809,652810,652812,652813,652814,652815,652816,652817,652818,652819,652820,652821,652822,652823,652824,652825,652826,652827,652828,652829,652830,652831,652832,652833,652834,652835,652836,652837,652838,652839)
        -- z-band AOS failers 'sqrt(pow(qc_fwhm,2)-pow(dimm2see,2)) > 0.7'
        and id not in (652692,652693,652694,652695,652702,652703,652704,652705,652706,652707,652709,652752,652753,652760,652762,652763,652764,652765,652773,652774,652775,652776,652777,652781,652782,652784,652785,652786,652787,652788,652789,652790,652791,652792,652801,652802,652811)
        -- t_eff values (careful about nulls)
        and not (qc_teff < 0.1 and date < '2017/08/07 12:00:00')
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class AlfredoFieldArray(FieldArray):
    """ Array of BLISS fields """
    PROGRAM  = 'eRosita'
    PROPID   = '2017A-0388'
    PROPOSER = 'Alfredo Zenteno'

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID
    SISPI_DICT["proposer"] = PROPOSER

    OBJECT_FMT = '%s'
    SEQID_FMT = '%(DATE)s'
    BANDS = BANDS

    @classmethod
    def query(cls, **kwargs):
        """ Generate the database query.

        Parameters:
        -----------
        kwargs : Keyword arguments to fill the query.

        Returns:
        --------
        query  : The query string.
        """
        from obztak.utils.date import setdefaults
        defaults = dict(propid=cls.SISPI_DICT['propid'], limit='')
        kwargs = setdefaults(kwargs,copy.deepcopy(defaults))

        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW
        FROM exposure where propid = '%(propid)s' and exptime > 59
        and discard = False and delivered = True and flavor = 'object'
        and qc_teff > 0.3
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class BlissScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"bliss-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"bliss-target-fields.csv")),
    ])
    FieldType = BlissFieldArray

# Below are several utility functions for plotting the survey
