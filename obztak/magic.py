#!/usr/bin/env python
"""
Code related to the DECam Dwarf Galaxy Survey.
"""
from __future__ import print_function

import os,sys
import logging
import copy
from collections import OrderedDict as odict

import pandas as pd
import numpy as np
import fitsio
from matplotlib.path import Path

from obztak.field import FieldArray, SISPI_DICT, SEP
from obztak.survey import Survey
from obztak.scheduler import Scheduler
from obztak.tactician import Tactician

from obztak.utils.projector import cel2gal, angsep
from obztak.utils import constants
from obztak.utils import fileio
from obztak.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS,DECAM
from obztak.utils.constants import COLORS, CMAPS
from obztak.utils.date import datestring, setdefaults, nite2utc,utc2nite,datestr

NAME    = 'MAGIC'
PROGRAM = NAME.lower()
PROPID  = '2023B-646244'
PROPOSER = 'Chiti'
BANDS = ['N395']
TILINGS = [1]
DONE = -1

class MagicSurvey(Survey):
    """ Survey sublcass for BLISS. """

    nights_2023B = [
        ['2023/08/10', 'first'],
        ['2023/08/11', 'full'],
        ['2023/08/13', 'full'],
        ['2023/08/14', 'full'],
        ['2023/08/16', 'first'],
        ['2023/08/17', 'first'],
        ['2023/08/19', 'first'],
        ['2023/08/20', 'first'],
        ['2023/08/22', 'second'],
        ['2023/08/23', 'second'],
    ]

    extra_nights = []

    nights = nights_2023B + extra_nights

    def prepare_fields(self, infile=None, outfile=None, plot=True, **kwargs):
        """ Create the list of fields to be targeted by this survey.

        Parameters:
        -----------
        infile : File containing all possible field locations.
        outfile: Output file of selected fields
        plot   : Create an output plot of selected fields.

        Returns:
        --------
        fields : A FieldArray of the selected fields.
        """

        if infile is None:
            infile = fileio.get_datafile('decam-tiles-bliss-v1.fits.gz')
        logging.info("Reading tiles from: %s"%os.path.basename(infile))
        data = fitsio.read(infile)

        fields  = self.create_magic_fields(data)

        logging.info("Masking bright stars...")
        mask = self.bright_stars(fields['RA'],fields['DEC'])
        fields['PRIORITY'][mask] = DONE

        # Can't go lower
        mask = fields['DEC'] < -89.0
        fields['PRIORITY'][mask] = DONE

        # Exclusion
        exclude = []
        fields['PRIORITY'][np.in1d(fields.unique_id,exclude)] = DONE

        if plot:
            import pylab as plt
            import skymap.survey
            plt.ion()

            sel = [fields['PRIORITY'] > 0]

            plt.figure()
            smap = skymap.survey.MaglitesSkymap()
            smap.draw_fields(fields[sel],alpha=0.3,edgecolor='none')
            smap.draw_des(color='r')
            smap.draw_milky_way()
            #smap.draw_smash()

            plt.figure()
            smap = skymap.survey.SurveyMcBryde()
            smap.draw_fields(fields[sel],alpha=0.3,edgecolor='none')
            smap.draw_des(color='r')
            smap.draw_milky_way()
            #smap.draw_smash()

            if outfile:
                plt.savefig(os.path.splitext(outfile)[0]+'.png',bbox_inches='tight')
            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile:
            print("Writing %s..."%outfile)
            fields.write(outfile)

        return fields

    @classmethod
    def update_covered_fields(cls, fields):
        """ Update the priority of covered fields. """
        # Currently NOOP
        fields = copy.deepcopy(fields)
        frac, depth = cls.covered(fields)
        done = (fields['PRIORITY'] == DONE)
        print("Found %i exposures already done."%done.sum())

        return fields

    def create_magic_fields(self, data, plot=False):
        """ Create the wide field observations """
        logging.info("Creating MAGIC fields...")
        EXPTIME = [720]
        TILINGS = [1]

        nhexes = len(np.unique(data['TILEID']))
        nbands = len(BANDS)

        nfields = len(data)*nbands

        logging.info("  Number of hexes: %d"%nhexes)
        logging.info("  Filters: %s"%BANDS)
        logging.info("  Exposure time: %s"%EXPTIME)
        logging.info("  Tilings: %s"%TILINGS)

        fields = FieldArray(nfields)
        fields['PROGRAM'] = PROGRAM
        fields['HEX'] = np.repeat(data['TILEID'],nbands)
        fields['TILING'] = np.repeat(data['PASS'],nbands)
        fields['RA'] = np.repeat(data['RA'],nbands)
        fields['DEC'] = np.repeat(data['DEC'],nbands)

        fields['FILTER'] = np.tile(BANDS,len(data))
        fields['EXPTIME'] = np.tile(EXPTIME,len(data))
        fields['PRIORITY'] = fields['TILING']

        sel  = self.footprintMAGIC(fields['RA'],fields['DEC'])
        sel &= (fields['TILING'] == 1)

        fields = fields[sel]

        if plot: self.plot_depth(fields,depth,f'{PROGRAM}-%s-gt%i.png')

        logging.info("Number of target fields: %d"%len(fields))

        outfile = f'{PROGRAM}-fields.fits.fz'
        logging.info("Writing %s..."%outfile)
        fields.write(outfile,clobber=True)

        return fields

    @staticmethod
    def footprintMAGIC(ra,dec):
        """ Selecting wide-field exposures plane """
        ra,dec = np.copy(ra), np.copy(dec)
        ra -= 360 * (ra > 180)

        vertices = [
            [90,  -30],
            [ 0,  -30],
            [-50, -30],
            [-60, -32],
            [-80, -45],
            [-80, -60],
            [-80, -80],
            [  0, -80],
            [120, -80],
            [120, -60],
            [105, -45],
            [ 90, -30],
        ]
        path = Path(vertices)
        points = np.vstack([ra,dec]).T
        sel = path.contains_points(points)

        return sel

    def plot_depth(self, fields, depth, outbase, proj='mcbryde', **kwargs):
        import skymap, skymap.survey
        import pylab as plt
        bands = np.unique(fields['FILTER'])
        ra,dec = fields['RA'],fields['DEC']
        for b in bands:
            sel = fields['FILTER']==b
            for d in np.unique(fields[sel]['EXPTIME']*fields[sel]['TILING'])[:-1]:
                plt.figure()
                if proj == 'mcbryde': smap = skymap.McBrydeSkymap()
                elif proj == 'maglites': smap = skymap.survey.MaglitesSkymap()
                smap.scatter(*smap(ra[sel],dec[sel]),c=depth[sel],vmax=d,edgecolor='none',s=3)
                smap.draw_lmc(fc='none')
                smap.draw_smc(fc='none')
                plt.colorbar()
                plt.savefig(outbase%(b,d),bbox_inches='tight')
                plt.close()


class MagicFieldArray(FieldArray):
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
        #date_column = 'date' or 'to_timestamp(utc_beg)'
        defaults['date_column'] = 'date'
        kwargs = setdefaults(kwargs,copy.deepcopy(defaults))

        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(%(date_column)s, 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW, PROGRAM
        FROM exposure where propid in ('%(propid)s')
        and exptime > 89
        and discard = False and delivered = True and flavor = 'object'
        and object LIKE '%(object_fmt)s%%'
        and (
             (COALESCE(qc_teff,-1) NOT BETWEEN 0 and 0.2
             AND COALESCE(qc_fwhm,1) BETWEEN 0.5 and 1.5)
             OR %(date_column)s  > (now() - interval '14 hours')
        )
        ORDER BY %(date_column)s %(limit)s
        """%kwargs
        return query


class MagicScheduler(Scheduler):
    _defaults = odict(list(Scheduler._defaults.items()) + [
        ('tactician','coverage'),
        ('windows',fileio.get_datafile("magic-windows-20230714.csv.gz")),
        ('targets',fileio.get_datafile("magic-fields-20230714.csv.gz")),
    ])

    FieldType = MagicFieldArray


class MagicTactician(Tactician):
    CONDITIONS = odict([
        (None,       [1.0, 1.4]),
        ('magic',    [1.0, 1.4]),
    ])

    def __init__(self, *args, **kwargs):
        super(MagicTactician,self).__init__(*args,**kwargs)
        #Default to mode 'wide' if no mode in kwargs
        self.mode = kwargs.get('mode','magic')

    @property
    def viable_fields(self):
        viable = super(MagicTactician,self).viable_fields
        viable &= (self.fields['PRIORITY'] >= 0)
        return viable

    def skybright_select(self):
        """Select fields based on skybrightness and band.
        All MAGIC observations should be dark.

        Parameters:
        -----------
        None

        Returns:
        --------
        sel : boolean selection
        """
        sel = np.ones(len(self.fields),dtype=bool)
        return sel

    @property
    def weight(self):
        """ Calculate the weight from set of programs. """
        return self.weight_magic()

    def weight_magic(self):
        """ Calculate the field weight for the WIDE survey.

        Parameters
        ----------
        None

        Returns
        -------
        weight : array of weights per field
        """
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        weight = np.zeros(len(sel))

        # Prefer DES region
        sel &= (self.fields['DEC'] > -50) & (self.fields['DEC'] < -40)
        sel &= ((self.fields['RA'] > 290) | (self.fields['RA'] < 90))
        inDES = Survey.footprintDES(self.fields['RA'],self.fields['DEC'])
        weight += 1e3 * ~inDES

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS['magic']

        if False: # Don't use seeing
            sel &= ((airmass > airmass_min) & (airmass < airmass_max))
        elif self.fwhm <= 1.4:
            sel &= ((airmass > airmass_min) & (airmass < airmass_max))
        else:
            sel &= ((airmass > airmass_min) & (airmass < 1.3))

        # Higher weight for fields close to the moon (when up)
        # angle = 50 -> weight = 6.4
        # Moon angle constraints (viable fields sets moon_angle > 20.)
        if (self.moon.alt > -0.04) and (self.moon.phase >= 10):
            #moon_limit = np.min(20 + self.moon.phase/2., 40)
            moon_limit = 50. + (self.moon.phase/10.)
            sel &= (moon_angle > moon_limit)

            #weight += 100 * (35./moon_angle)**3
            #weight += 10 * (35./moon_angle)**3
            weight += 1 * (35./moon_angle)**3

        # Higher weight for rising fields (higher hour angle)
        # HA [min,max] = [-53,54] (for airmass 1.4)
        weight += 10.0 * self.hour_angle
        #weight += 1.0 * self.hour_angle
        #weight += 0.1 * self.hour_angle

        # Higher weight for larger slews
        # slew = 10 deg -> weight = 1e2
        weight += self.slew**2
        #weight += self.slew
        #weight += 1e3 * self.slew

        # Higher weight for higher airmass
        # airmass = 1.4 -> weight = 6.4
        weight += 100. * (airmass - 1.)**3
        #weight += 1e3 * (airmass - 1.)**2

        # Set infinite weight to all disallowed fields
        weight[~sel] = np.inf

        return weight

    def select_index(self):
        weight = self.weight
        index = np.array([np.argmin(weight)],dtype=int)
        if np.any(~np.isfinite(weight[index])):
            plot = (logging.getLogger().getEffectiveLevel()==logging.DEBUG)
            msg = "Infinite weight selected..."
            logging.warn(msg)
            logging.info(">>> To plot fields enter 'plot=True'")
            logging.info(">>> Enter 'c' to continue")
            import pdb; pdb.set_trace()
            if plot:
                import obztak.utils.ortho, pylab as plt
                airmass = self.CONDITIONS[self.mode][1]
                if len(self.completed_fields):
                    last_field =self.completed_fields[-1]
                else:
                    last_field = self.fields[-1]
                smap = obztak.utils.ortho.plotFields(last_field,self.fields,self.completed_fields,options_basemap=dict(airmass=airmass))
                plt.show()
                logging.info(">>> Enter 'c' to continue")
                pdb.set_trace()
            raise ValueError(msg)

        return index
