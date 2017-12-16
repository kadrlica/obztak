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

NAME    = 'PALS'
PROGRAM = NAME.lower()
PROPID  = '2018A-0273'
PROPOSER = 'Dawson'
BANDS = ['g','r']
TILINGS = 124

class PalsSurvey(Survey):
    """ Survey sublcass for BLISS. """

    # 2018A SCHEDULED (data/pals-windows.csv is actually used)
    nights = [
        ['2018/02/15', 'full'],
        ['2018/02/16', 'full'],
        ['2018/03/11', 'full'],
        ['2018/03/12', 'full'],
        ['2018/07/18', 'full'],
        ['2018/07/19', 'full'],
        ]

    def prepare_fields(self, infile=None, outfile=None, mode='none', plot=True, smcnod=False):
        """Create the list of fields to be targeted by the PALS survey.

        Selection starts with 3 regions:

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
        lmc    = self.lmc(fields['RA'],fields['DEC'])
        smc    = self.smc(fields['RA'],fields['DEC'])
        bulge  = self.bulge(fields['RA'],fields['DEC'])

        sel = (lmc | smc | bulge)

        # Apply telescope constraints
        sel &= (fields['DEC'] > constants.SOUTHERN_REACH)

        # Apply covered fields
        #sel &= self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]
        # Apply covered fields (but take everything in P9 region)

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
        sel = PalsSurvey.lmc(ra,dec)
        # Alfredo's eRosita survey
        sel |= PalsSurvey.smc(ra,dec)
        # Planet 9 Region
        sel |= PalsSurvey.gc(ra,dec)

        return sel

    @staticmethod
    def bulge(ra,dec):
        """Bulge fields"""
        RA,DEC = 273.91666667, -28.5
        dRA,dDEC = 7.5, 7.5
        sel = (np.abs(ra - RA) < dRA/2.) & (np.abs(dec - DEC) < dDEC/2.)
        return sel

    @staticmethod
    def lmc(ra,dec):
        """LMC fields"""
        RA,DEC = 93.61166667, -68.30319444
        dRA,dDEC = 7.0, 8.0
        sel = (np.abs(ra - RA) < dRA/2.) & (np.abs(dec - DEC) < dDEC/2.)
        return sel

    @staticmethod
    def smc(ra,dec):
        """SMC fields"""
        RA,DEC = 13.15833333, -72.80027778
        dRA,dDEC = 3.0, 4.5
        sel = (np.abs(ra - RA) < dRA/2.) & (np.abs(dec - DEC) < dDEC/2.)
        return sel

class PalsFieldArray(FieldArray):
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
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class PalsScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"pals-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"pals-target-fields.csv")),
    ])
    FieldType = PalsFieldArray
