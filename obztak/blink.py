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
from obztak.utils.date import setdefaults

NAME    = 'BLINK'
PROGRAM = NAME.lower()
PROPID  = '2018A-0914'
PROPOSER = 'Makler'
BANDS = ['i']
TILINGS = 2

class BlinkSurvey(Survey):
    """ Survey sublcass for BLINK. """

    # 2018A SCHEDULED (blink-windows.csv is actually used for scheduling)
    nights = [
        ['2018/02/04', 'second'], 
        ['2018/02/05', 'second'], 
        ['2018/04/11', 'second'], 
        ['2018/04/12', 'second'], 
        ['2018/04/13', 'second'],
        ['2018/04/17', 'second']
        ]
    
    def prepare_fields(self, infile=None, outfile=None, mode='bliss_rotate', plot=True, smcnod=False):
        """Create the list of fields to be targeted by the BLINK survey.

        Selection starts with 1 region:
        

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
        _id      = data['TILEID']
        ra       = data['RA']
        dec      = data['DEC']

        nhexes = len(data)
        #ntilings = len(DECAM_DITHERS)
        ntilings = TILINGS
        nbands = len(BANDS)
        nfields = nhexes*nbands*ntilings

        logging.info("Number of hexes: %d"%nhexes)
        logging.info("Number of tilings: %d"%ntilings)
        logging.info("Number of filters: %d"%nbands)

        fields = FieldArray(nfields)
        fields['HEX'] = np.tile(np.repeat(_id,nbands),ntilings)
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
        sel = self.blink(fields['RA'],fields['DEC'])
        fields['PRIORITY'][sel] = 1

        # Apply telescope constraints
        sel &= (fields['DEC'] > constants.SOUTHERN_REACH)

        # Apply covered fields
        #sel &= self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]
        # Apply covered fields (but take everything in P9 region)
        uncovered = self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]
        sel &= uncovered

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
                outfig = os.path.splitext(outfile)[0]+'_mcbride.png'
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
        # The high-probability region for LIGO/MW dwarfs
        return BlinkSurvey.blink(ra,dec)

    @staticmethod
    def blink(ra,dec):
        """ Initial guess at BLINK survey region """
        sel  = (( dec > 3) & (dec < 10)) #(( dec > -1) & (dec < 8)) | ((dec > 3) & (dec < 8)) | ((dec > 4) & (dec < 11))
        sel &= (( ra > 133.7) & (ra < 230)) #(( ra > 141.5) & (ra < 156)) | ((ra > 133.7) & (ra < 156)) | ((ra > 170) & (ra < 238))
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
        #dirname = '/Users/kadrlica/bliss/observing/data'
        dirname = fileio.get_datadir() #'/home/s1/kadrlica/software/obztak/obztak/data'
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

class BlinkFieldArray(FieldArray):
    """ Array of BLINK fields """
    PROGRAM  = PROGRAM
    PROPID   = PROPID
    PROPOSER = PROPOSER

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID
    SISPI_DICT["proposer"] = PROPOSER

    OBJECT_FMT = 'BLINK field' + SEP + ' %s'
    SEQID_FMT = 'BLINK scheduled' + SEP + ' %(DATE)s'
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
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class BlinkScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"blink-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"blink-target-fields.csv")),
    ])
    FieldType = BlinkFieldArray
