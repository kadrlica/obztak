#!/usr/bin/env python
"""
Code related to the Magellanic Satellites Survey (MagLiteS).
"""
import os,sys
import logging
import copy
from collections import OrderedDict as odict

import numpy as np

from obztak.field import FieldArray, SISPI_DICT, SEP
from obztak.survey import Survey
from obztak.scheduler import Scheduler

from obztak.utils.projector import cel2gal, angsep
from obztak.utils import constants
from obztak.utils import fileio

PROGRAM = 'bliss'
PROPID  = '2017A-0260'
BANDS = ['g','r','i','z']
TILINGS = 2

class BlissSurvey(Survey):
    """ Survey sublcass for BLISS. """

    # 2017A ACTUAL
    nights = [
        ['2017/02/07', 'second'],
        ['2017/02/08', 'second'],
        ['2017/02/13', 'second'],
        ['2017/02/14', 'second'],
        ['2017/02/19', 'second'],
        ['2017/02/20', 'second'],
        ['2017/03/07', 'first'],
        ['2017/03/17', 'second'],
        ['2017/04/13', 'second'],
        ['2017/04/14', 'second'],
        ['2017/05/02', 'first'],
        ['2017/05/30', 'first'],
        ['2017/05/31', 'first'],              
        ['2017/06/01', 'first'],
        ['2017/06/02', 'full'],
        ['2017/06/03', 'full'],
        ['2017/06/04', 'full'],
        ['2017/07/01', 'second'],
        ['2017/07/02', 'second'],
        ['2017/07/15', 'second'],
        ]

    alfredo_nights = [
        ['2017/03/06', 'full'],
        ['2017/03/07','second'],
        ['2017/06/21','first'],
        ]
    
    def prepare_fields(self, infile=None, outfile=None, mode='smash_dither', plot=True, smcnod=False):
        """ Create the list of fields to be targeted by this survey.

        Parameters:
        -----------
        infile : File containing all possible field locations.
        outfile: Output file of selected fields
        mode   : Mode for dithering: 'smash_dither', 'smash_rotate', 'decam_dither', 'none'
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

        if infile is None:
            infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
        data = np.recfromtxt(infile, names=True)

        # Apply footprint selection after tiling/dither
        #sel = obztak.utils.projector.footprint(data['RA'],data['DEC'])

        # This is currently a non-op
        smash_id = data['ID']
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
        fields['HEX'] = np.tile(np.repeat(smash_id,nbands),ntilings)
        fields['PRIORITY'].fill(1)
        fields['TILING'] = np.repeat(np.arange(1,ntilings+1),nhexes*nbands)
        fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)

        #for i in range(ntilings):
        for i,offset in enumerate(OFFSETS):
            idx0 = i*nhexes*nbands
            idx1 = idx0+nhexes*nbands
            ra_dither,dec_dither = dither(ra,dec,offset[0],offset[1])
            fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
            fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

        # Apply footprint selection after tiling/dither
        sel = self.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION
        sel &= (fields['DEC'] > constants.SOUTHERN_REACH)

        # Apply covered fields
        sel &= self.uncovered(fields['RA'],fields['DEC'],fields['FILTER'])[0]

        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            from obztak.utils.ortho import makePlot, safeProj
            from obztak.utils.ortho import DECamOrtho, DECamMcBride

            kwargs = dict(edgecolor='none',cmap='Spectral',vmin=0,vmax=ntilings)

            fig, bmap = makePlot('2017/02/08 03:00',center=(150,-30),
                                    airmass=False,moon=False)

            bmap.scatter(*bmap.proj(fields['RA'],fields['DEC']),
                          c=fields['TILING'],s=50,**kwargs)
            bmap.draw_galaxy(10)
            colorbar = plt.colorbar(label='Tiling')

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'_ortho.png'
                fig.savefig(outfig,bbox_inches='tight')

            plt.figure()
            bmap = DECamMcBride()
            bmap.draw_galaxy(10)
            bmap.draw_des()
            bmap.scatter(*bmap.proj(fields['RA'],fields['DEC']),
                          c=fields['TILING'],s=15,**kwargs)
            plt.colorbar(label='Tiling',orientation='horizontal',fraction=0.1,
                         pad=0.1,aspect=30,shrink=0.8)
            if outfile:
                outfig = os.path.splitext(outfile)[0]+'_mcbride.png'
                fig.savefig(outfig,bbox_inches='tight')

            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile: fields.write(outfile)

        return fields

    @staticmethod
    def footprint(ra,dec):
        l, b = cel2gal(ra, dec)
        # Starting to the north of the Galactic plane
        sel = (ra > 120) & (ra < 360)
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
    def alfredo(ra,dec):
        l, b = cel2gal(ra, dec)
        sel = ((b > 20) & (b < 30))
        sel &= ((l < 360) & (l > 240))
        sel &= ((ra > 135) & (ra < 245))
        sel &= ((dec > -30) & (dec < -10))
        return sel

    @staticmethod
    def uncovered(ra,dec,band):
        import healpy as hp
        dirname = '/home/s1/kadrlica/projects/bliss/v0/data'
        basename = 'decam_coverage_90s_%s_n1024.fits.gz'
        
        sel = np.ones_like(ra,dtype=bool)
        frac = np.zeros_like(ra,dtype=float)
        ra,dec,band=np.atleast_1d(ra),np.atleast_1d(dec),np.atleast_1d(band)

        for b in np.unique(band):
            idx = (band==b)
            filename = os.path.join(dirname,basename%b)
            logging.info("Reading %s..."%os.path.basename(filename))
            skymap = hp.read_map(filename,verbose=True)
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
    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID

    OBJECT_FMT = 'BLISS field' + SEP + ' %s'
    SEQID_FMT = 'BLISS scheduled' + SEP + ' %(DATE)s'
    BANDS = BANDS


class BlissScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"bliss-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"bliss-target-fields.csv")),
    ])
    FieldType = BlissFieldArray
