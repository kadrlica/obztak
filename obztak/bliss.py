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

from obztak.utils.projector import cel2gal
from obztak.utils import constants
from obztak.utils import fileio

PROGRAM = 'bliss'
PROPID  = '2017A-0260'
BANDS = ['g','r','i','z']

class BlissSurvey(Survey):
    """ Survey sublcass for BLISS. """

    # 2017A ACTUAL
    nights = [['2017/02/07', 'second'],
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
            TILINGS = [(0,0),(0,0),(0,0)]
        elif mode.lower() == 'smash_dither':
            TILINGS = [(0,0), (1.0,0.0), (-1.0,0.0)]
            dither = self.smash_dither
        elif mode.lower() == 'smash_rotate':
            TILINGS = [(0,0), (0.75,0.75), (-0.75,0.75)]
            dither = self.smash_rotate
        elif mode.lower() == 'decam_dither':
            TILINGS = [(0., 0.),(8/3.*CCD_X, -11/3.*CCD_Y),
                       (8/3.*CCD_X, 8/3.*CCD_Y)]
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
        ntilings = len(TILINGS)
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
        for i,tiling in enumerate(TILINGS):
            idx0 = i*nhexes*nbands
            idx1 = idx0+nhexes*nbands
            ra_dither,dec_dither = dither(ra,dec,tiling[0],tiling[1])
            fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
            fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

        # Apply footprint selection after tiling/dither
        sel = self.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION
        sel = sel & (fields['DEC'] > constants.SOUTHERN_REACH)
        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            from obztak.utils.ortho import makePlot, safeProj

            fig, basemap = makePlot('2016/2/11 03:00',center=(180,-30),airmass=False,moon=False)

            proj = safeProj(basemap,fields['RA'],fields['DEC'])
            basemap.scatter(*proj, c=fields['TILING'], edgecolor='none', s=50, cmap='Spectral',vmin=0,vmax=len(TILINGS))
            colorbar = plt.colorbar(label='Tiling')

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'.png'
                fig.savefig(outfig,bbox_inches='tight')
            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile: fields.write(outfile)

        return fields

    @staticmethod
    def footprint(ra,dec):
        l, b = cel2gal(ra, dec)
        sel  = (np.fabs(b) > 10.)
        sel &= (ra > 120) & (ra < 180)
        sel &= (dec < 10 ) & (dec > -30)
        return sel

class BlissFieldArray(FieldArray):
    """ Array of BLISS fields """
    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID

    OBJECT_FMT = 'BLISS field' + SEP + ' %s'
    SEQID_FMT = 'BLISS scheduled' + SEP + ' %(DATE)s'
    BANDS = BANDS

class BlissScheduler(Scheduler):
    _defaults = odict([
        ('windows',os.path.join(fileio.get_datadir(),"bliss-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"bliss-target-fields.csv")),
    ])
    FieldType = BlissFieldArray
