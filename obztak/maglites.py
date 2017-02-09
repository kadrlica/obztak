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

from obztak.utils import constants
from obztak.utils.constants import SMASH_POLE,CCD_X,CCD_Y,STANDARDS
from obztak.utils.projector import cel2gal, angsep
from obztak.utils.date import datestring, setdefaults
from obztak.utils import fileio

NAME = 'MagLiteS'
PROGRAM = NAME.lower()
PROPID  = '2016A-0366'
PROPOSER = 'Bechtol'
BANDS = ['g','r']

class MaglitesSurvey(Survey):
    """ Survey sublcass for MagLiteS. """

    # 2016A ACTUAL
    nights_2016A = [
        ['2016/02/11', 'second'],
        ['2016/02/12', 'second'],
        ['2016/02/13', 'second'],
        ['2016/02/14', 'second'],
        ['2016/02/15', 'second'],
        ['2016/02/16', 'second'],
        ['2016/06/27', 'full'],
        ['2016/06/28', 'full'],
        ['2016/06/29', 'full'],
        ]

    # 2017A ACTUAL
    nights_2017A = [
        ['2017/02/21', 'full'],
        ['2017/02/22', 'full'],
        ['2017/02/23', 'full'],
        ['2017/06/18', 'full'],
        ['2017/06/19', 'full'],
        ['2017/06/20', 'full']
        ]

    nights = nights_2016A + nights_2017A

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
            TILINGS = [(0,0),(0,0),(0,0),(0,0)]
        elif mode.lower() == 'smash_dither':
            TILINGS = [(0,0), (1.0,0.0), (-1.0,0.0), (0.0,-0.75)]
            dither = self.smash_dither
        elif mode.lower() == 'smash_rotate':
            TILINGS = [(0,0), (0.75,0.75), (-0.75,0.75), (0.0,-0.75)]
            dither = self.smash_rotate
        elif mode.lower() == 'decam_dither':
            TILINGS = [(0., 0.),(8/3.*CCD_X, -11/3.*CCD_Y),
                       (8/3.*CCD_X, 8/3.*CCD_Y),(-8/3.*CCD_X, 0.)]
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
        if smcnod:
            # Include SMC northern overdensity fields
            sel_smcnod = self.footprintSMCNOD(fields) # SMCNOD OPERATION
            sel = sel | sel_smcnod
            #sel = sel_smcnod
            fields['PRIORITY'][sel_smcnod] = 99
        #if True:
        #    # Include 'bridge' region between Magellanic Clouds
        #    sel_bridge = self.footprintBridge(fields['RA'],fields['DEC'])
        #    sel = sel | sel_bridge
        sel = sel & (fields['DEC'] > constants.SOUTHERN_REACH)
        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            import obztak.utils.ortho

            plt.ion()

            fig, basemap = obztak.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

            proj = obztak.utils.ortho.safeProj(basemap,fields['RA'],fields['DEC'])
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

        angsep_lmc = angsep(constants.RA_LMC, constants.DEC_LMC, ra, dec)
        angsep_smc = angsep(constants.RA_SMC, constants.DEC_SMC, ra, dec)
        sel = (np.fabs(b) > 10.) \
              & ((angsep_lmc < 30.) | (angsep_smc < 30.)) \
              & (dec < -55.) & (ra > 100.) & (ra < 300.)
        #sel = sel | ((dec < -65.) & (angsep_lmc > 5.) & (angsep_smc > 5.))
        sel = sel | ((dec < -65.) & (ra > 300.) & (ra < 360.)) # SMC
        sel = sel | (dec < -80.)

        return sel

    @staticmethod
    def footprintSMCNOD(fields):
        """
        Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
        """
        sel = np.in1d(fields['HEX'], constants.HEX_SMCNOD) \
              & np.in1d(fields['TILING'], constants.TILING_SMCNOD)
        return sel

    @staticmethod
    def footprintBridge(ra, dec):
        """
        Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
        """
        sel = (ra > 30.) & (ra < 60.) & (dec < -65.)
        return sel

class MaglitesFieldArray(FieldArray):
    PROGRAM  = PROGRAM
    PROPID   = PROPID
    PROPOSER = PROPOSER

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID
    SISPI_DICT["proposer"] = PROPOSER

    OBJECT_FMT = NAME.upper() + ' field'+SEP+' %s'
    SEQID_FMT  = NAME.upper() + ' scheduled'+SEP+' %(DATE)s'
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
        -- Discard exposures from Y1 with teff < 0.1
        and not (qc_teff < 0.1 and date < '2017/01/01')
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class MaglitesScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"maglites-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"maglites-target-fields.csv")),
    ])

    FieldType = MaglitesFieldArray

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
