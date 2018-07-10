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
import ephem

from obztak.field import FieldArray, SISPI_DICT, SEP
from obztak.survey import Survey
from obztak.scheduler import Scheduler
from obztak.tactician import Tactician

from obztak.utils import constants
from obztak.utils.constants import SMASH_POLE,CCD_X,CCD_Y,STANDARDS
from obztak.utils.projector import cel2gal, angsep
from obztak.utils.date import datestring, setdefaults
from obztak.utils import fileio

NAME = 'MagLiteS-II'
PROGRAM = NAME.lower()
PROPID  = '2018A-0242'
PROPOSER = 'Bechtol'
BANDS = ['r','i']
TILINGS = [1,2,3,4]

class Maglites2Survey(Survey):
    """ Survey sublcass for MagLiteS. """

    # 2018A ACTUAL
    nights_2018A = [
        ['2018/07/20', 'full'],
        ['2018/07/21', 'full'],
        ['2018/07/22', 'full'],
        ['2018/07/23', 'full'],
        ]

    # 2019A PREDICTED
    nights_2019A = [
        ['2019/02/24', 'full'],
        ['2019/02/25', 'full'],
        ['2019/02/26', 'full'],
        ['2019/02/27', 'full'],
        ['2019/06/24', 'full'],
        ['2019/06/25', 'full']
        ]

    nights = nights_2018A + nights_2019A

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
            #infile = os.path.join(fileio.get_datadir(),'decam-tiles_obstatus.fits')
            infile = os.path.join(fileio.get_datadir(),'decam-tiles-bliss.fits.gz')

        # For decals input file
        #data['TILEID'] = np.tile(data['TILEID'][data['PASS'] == 1],len(TILINGS))

        data = fitsio.read(infile)

        # Only take tilings of interest
        data = data[np.in1d(data['PASS'],TILINGS)]

        # Number of unique tiles
        nhexes = len(np.unique(data['TILEID']))
        ntilings = len(TILINGS)
        nbands = len(BANDS)
        nfields = nhexes*nbands*ntilings

        assert nfields == len(data)*nbands

        logging.info("Number of hexes: %d"%nhexes)
        logging.info("Number of tilings: %d"%ntilings)
        logging.info("Number of filters: %d"%nbands)

        fields = FieldArray(nfields)
        fields['HEX'] = np.repeat(data['TILEID'],nbands)
        fields['TILING'] = np.repeat(data['PASS'],nbands)
        fields['PRIORITY'] = fields['TILING']
        fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)
        fields['RA'] = np.repeat(data['RA'],nbands)
        fields['DEC'] = np.repeat(data['DEC'],nbands)

        # Apply footprint selection after tiling/dither
        sel = self.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION
        sel = sel & (fields['DEC'] > constants.SOUTHERN_REACH)
        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            import obztak.utils.ortho

            plt.ion()

            fig, basemap = obztak.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

            proj = basemap.proj(fields['RA'],fields['DEC'])
            basemap.scatter(*proj, c=fields['TILING'], edgecolor='none', s=50, cmap='Spectral',vmin=0,vmax=ntilings)
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
        import matplotlib.path
        filename = os.path.join(fileio.get_datadir(),'maglitesII-poly.txt')
        data = np.genfromtxt(filename,names=['ra','dec','poly'])
        paths = []
        ra -= 360 * (ra > 180)
        for p in np.unique(data['poly']):
            poly = data[data['poly'] == p]
            vertices = np.vstack(np.vstack([poly['ra'],poly['dec']])).T
            paths.append(matplotlib.path.Path(vertices))
        sel = np.sum([p.contains_points(np.vstack([ra,dec]).T) for p in paths],axis=0) > 0
        return sel

    @staticmethod
    def footprintBridge(ra, dec):
        """
        Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
        """
        sel = (ra > 30.) & (ra < 60.) & (dec < -65.)
        return sel

class Maglites2FieldArray(FieldArray):
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

        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW
        FROM exposure where propid = '%(propid)s' and exptime > 89
        and discard = False and delivered = True and flavor = 'object'
        and object like '%(object_fmt)s%%'
        -- Discard exposures with teff < 0.1
        --and not (qc_teff < 0.1 and date < '2017/01/01')
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class Maglites2Scheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"maglites2-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"maglites2-target-fields.csv")),
    ])

    FieldType = Maglites2FieldArray


class Maglites2Tactician(Tactician):
    CONDITIONS = odict([
        (None,        [1.0, 2.0]),
        ('maglites2', [1.0, 2.0]),
    ])

    def __init__(self, *args, **kwargs):
        super(Maglites2Tactician,self).__init__(*args,**kwargs)
        self.mode = kwargs.get('mode',None)

    @property
    def weight(self):
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        weight = np.zeros(len(sel))

        # Moon angle constraints
        moon_limit = 30.
        sel &= (moon_angle > moon_limit)

        # Moon band constraints
        if (self.moon.phase >= 80) and (self.moon.alt > 0.0):
            # Moon is very bright; only do i
            sel &= (np.char.count('i',self.fields['FILTER']) > 0)
            # Allow i,z but prefer z
            #sel &= (np.char.count('iz',self.fields['FILTER']) > 0)
            #weight += 1e2 * (np.char.count('i',self.fields['FILTER']) > 0)
        elif (self.moon.phase >= 45) and (self.moon.alt > 0.0):
            # Moon is more than half full; do i,z
            sel &= (np.char.count('ri',self.fields['FILTER']) > 0)
        else:
            # Moon is faint or down; do g,r (unless none available)
            sel &= (np.char.count('r',self.fields['FILTER']) > 0)
            #weight += 1e8 * (np.char.count('iz',self.fields['FILTER']) > 0)

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS[self.mode]
        sel &= ((airmass > airmass_min) & (airmass < airmass_max))

        # Set the weights for each field. Lower weight means more favorable.

        # Higher weight for rising fields (higher hour angle)
        # HA [min,max] = [-53,54] (for airmass 1.4)
        #weight += 5.0 * self.hour_angle
        weight += 1.0 * self.hour_angle
        #weight += 0.1 * self.hour_angle

        # Higher weight for larger slews
        # slew = 10 deg -> weight = 1e2
        weight += self.slew**2
        #weight += self.slew
        #weight += 1e3 * self.slew

        # Higher weight for higher airmass
        # airmass = 1.4 -> weight = 6.4
        weight += 100. * (airmass - 1.)**3

        # Higher weight for fields close to the moon (when up)
        # angle = 50 -> weight = 6.4
        if (self.moon.alt > -0.04):
            #weight += 100 * (35./moon_angle)**3
            #weight += 10 * (35./moon_angle)**3
            weight += 1 * (35./moon_angle)**3

        # Try hard to do the first tiling
        weight += 1e6 * (self.fields['TILING'] - 1)

        # Set infinite weight to all disallowed fields
        weight[~sel] = np.inf

        return weight

    def select_index(self):
        weight = self.weight
        index = np.array([np.argmin(weight)],dtype=int)
        if np.any(~np.isfinite(weight[index])):
        #if True:
            msg = "Infinite weight selected"
            print(msg)
            import obztak.utils.ortho, pylab as plt
            airmass_min, airmass_max = self.CONDITIONS[self.mode]
            bmap = obztak.utils.ortho.plotFields(self.completed_fields[-1],self.fields,self.completed_fields,options_basemap=dict(airmass=airmass_max))
            import pdb; pdb.set_trace()
            raise ValueError(msg)

        return index

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
