#!/usr/bin/env python
"""
Code related to the DECam Dwarf Galaxy Survey.
"""
import os,sys
import logging
import copy
from collections import OrderedDict as odict

import pandas as pd
import numpy as np
import fitsio

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

NAME    = 'DELVE'
PROGRAM = NAME.lower()
PROPID  = '2019A-0305'
PROPOSER = 'Drlica-Wagner'
BANDS = ['g','r','i']
TILINGS = [1,2,3,4]
DONE = -1

class DelveSurvey(Survey):
    """ Survey sublcass for BLISS. """

    # 2019A SCHEDULED (bliss-windows.csv is actually used)
    nights = [
        ['2019/02/07', 'second'], # phase=%, set=
        ['2019/02/08', 'second'], # phase=%, set=
        ['2019/02/09', 'second'], # phase=%, set=
        ['2019/02/12', 'full  '], # phase=%, set=
        ['2019/02/13', 'full  '], # phase=%, set=
        ['2019/02/14', 'second'], # phase=%, set=
        ['2019/02/15', 'full  '], # phase=%, set=
        ['2019/02/24', 'second'], # phase=%, set=
        ['2019/02/25', 'second'], # phase=%, set=
        ['2019/02/26', 'second'], # phase=%, set=
        ['2019/02/27', 'second'], # phase=%, set=
        ['2019/02/28', 'second'], # phase=%, set=
        ['2019/03/01', 'second'], # phase=%, set=
        ['2019/05/12', 'full  '], # phase=%, set=
        ['2019/05/13', 'full  '], # phase=%, set=
        ['2019/05/28', 'second'], # phase=%, set=
        ['2019/05/29', 'second'], # phase=%, set=
        ['2019/05/30', 'second'], # phase=%, set=
        ['2019/05/31', 'second'], # phase=%, set=
        ['2019/06/01', 'second'], # phase=%, set=
        ['2019/06/02', 'second'], # phase=%, set=
        ['2019/06/03', 'second'], # phase=%, set=
        ['2019/06/04', 'second'], # phase=%, set=
        ['2019/06/05', 'full  '], # phase=%, set=
        ['2019/06/06', 'full  '], # phase=%, set=
        ['2019/06/07', 'full  '], # phase=%, set=
        ['2019/06/08', 'full  '], # phase=%, set=
        ['2019/06/09', 'full  '], # phase=%, set=
        ['2019/06/23', 'second'], # phase=%, set=
        ['2019/06/24', 'second'], # phase=%, set=
        ['2019/06/25', 'second'], # phase=%, set=
        ['2019/06/26', 'second'], # phase=%, set=
        ['2019/06/27', 'second'], # phase=%, set=
        ['2019/06/28', 'second'], # phase=%, set=
        ]

    extra_nights = []

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
            #infile = fileio.get_datafile('decam-tiles-smash-v1.fits.gz')
            #infile = fileio.get_datafile('decam-tiles-decals-v1.fits.gz')
        logging.info("Reading tiles from: %s"%os.path.basename(infile))
        data = fitsio.read(infile)

        deep_fields = self.create_deep_fields(data)
        mc_fields   = self.create_mc_fields(data)
        wide_fields = self.create_wide_fields(data)

        fields = wide_fields + mc_fields + deep_fields

        if plot:
            import pylab as plt
            import skymap.survey
            plt.ion()

            sel = [fields['PRIORITY'] > 0]

            plt.figure()
            smap = skymap.survey.MaglitesSkymap()
            smap.draw_fields(fields[sel],alpha=0.3,edgecolor='none')
            smap.draw_des(c='r')
            smap.draw_milky_way()
            smap.draw_smash()

            plt.figure()
            smap = skymap.survey.SurveyMcBryde()
            smap.draw_fields(fields[sel],alpha=0.3,edgecolor='none')
            smap.draw_des(c='r')
            smap.draw_milky_way()
            smap.draw_smash()

            if outfile:
                plt.savefig(os.path.splitext(outfile)[0]+'.png',bbox_inches='tight')
            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile: fields.write(outfile)

        return fields

    def create_wide_fields(self, data, plot=False):
        """ Create the wide field observations """
        logging.info("Creating DEEP fields...")
        BANDS = ['g','i']
        EXPTIME = [90,90]
        TILINGS = [4,4]
        TEFF_MIN = pd.DataFrame(dict(FILTER=['g','i'],TEFF=[0.4,0.5]))

        nhexes = len(np.unique(data['TILEID']))
        nbands = len(BANDS)

        nfields = len(data)*nbands

        logging.info("  Number of hexes: %d"%nhexes)
        logging.info("  Filters: %s"%BANDS)
        logging.info("  Exposure time: %s"%EXPTIME)
        logging.info("  Tilings: %s"%TILINGS)

        fields = FieldArray(nfields)
        fields['PROGRAM'] = PROGRAM+'-wide'
        fields['HEX'] = np.repeat(data['TILEID'],nbands)
        fields['TILING'] = np.repeat(data['PASS'],nbands)
        fields['RA'] = np.repeat(data['RA'],nbands)
        fields['DEC'] = np.repeat(data['DEC'],nbands)

        fields['FILTER'] = np.tile(BANDS,len(data))
        fields['EXPTIME'] = np.tile(EXPTIME,len(data))
        fields['PRIORITY'] = fields['TILING']

        sel = self.footprintWIDE(fields['RA'],fields['DEC'])
        sel &= (~self.footprintMilkyWay(fields['RA'],fields['DEC']))
        sel &= (~self.footprintDES(fields['RA'],fields['DEC']))
        sel &= (~self.footprintSMASH(fields['RA'],fields['DEC'],angsep=0.75*DECAM))
        sel &= (~self.footprintMC(fields['RA'],fields['DEC']))
        # Avoid DEEP fields? No.
        #sel &= (~self.footprintDEEP(fields['RA'],fields['DEC']))

        fields = fields[sel]

        frac, depth = self.covered(fields)

        teffmin = pd.DataFrame(fields).merge(TEFF_MIN,on='FILTER').to_records()['TEFF']

        fields['PRIORITY'][depth > teffmin*fields['TILING']*fields['EXPTIME']] = DONE
        # Avoid MagLiteS-II for now
        fields['PRIORITY'][self.footprintMaglites2(fields['RA'],fields['DEC'])] = DONE

        if plot: self.plot_depth(fields,depth,'delve-wide-%s-gt%i.png')

        logging.info("Number of target fields: %d"%len(fields))

        outfile = 'delve-wide-fields.fits.fz'
        logging.info("Writing %s..."%outfile)
        fields.write(outfile,clobber=True)

        return fields


    def create_mc_fields(self, data, plot=False):
        """ Select fields around the LMC """
        logging.info("Creating MC fields...")
        BANDS = ['g','r','i']
        EXPTIME = [267,267,333]
        TILINGS = [4, 4, 4]
        TEFF_MIN = pd.DataFrame(dict(FILTER=['g','r','i'],TEFF=[0.3,0.3,0.45]))

        nhexes = len(np.unique(data['TILEID']))
        nbands = len(BANDS)

        nfields = len(data)*nbands

        logging.info("  Number of hexes: %d"%nhexes)
        logging.info("  Filters: %s"%BANDS)
        logging.info("  Exposure time: %s"%EXPTIME)
        logging.info("  Tilings: %s"%TILINGS)

        fields = FieldArray(nfields)
        fields['PROGRAM'] = PROGRAM+'-mc'
        fields['HEX'] = np.repeat(data['TILEID'],nbands)
        fields['TILING'] = np.repeat(data['PASS'],nbands)
        fields['RA'] = np.repeat(data['RA'],nbands)
        fields['DEC'] = np.repeat(data['DEC'],nbands)

        fields['FILTER'] = np.tile(BANDS,len(data))
        fields['EXPTIME'] =np.tile(EXPTIME,len(data))
        fields['PRIORITY'] = fields['TILING']

        sel = self.footprintMC(fields['RA'],fields['DEC'])
        sel &= (~self.footprintDES(fields['RA'],fields['DEC']))
        sel &= (~self.footprintSMASH(fields['RA'],fields['DEC'],angsep=0.75*DECAM))
        sel &= (~self.footprintMilkyWay(fields['RA'],fields['DEC']))

        fields = fields[sel]

        frac, depth = self.covered(fields)

        teffmin = pd.DataFrame(fields).merge(TEFF_MIN,on='FILTER').to_records()['TEFF']

        fields['PRIORITY'][depth > teffmin*fields['TILING']*fields['EXPTIME']] = DONE

        # Avoid MagLiteS-II for now
        fields['PRIORITY'][self.footprintMaglites2(fields['RA'],fields['DEC'])] = DONE

        if plot: self.plot_depth(fields,depth,'delve-mc-%s-gt%i.png',proj='maglites')

        logging.info("Number of target fields: %d"%len(fields))

        outfile = 'delve-mc-fields.fits.fz'
        logging.info("Writing %s..."%outfile)
        fields.write(outfile,clobber=True)

        return fields

    def create_deep_fields(self, data, plot=False):
        logging.info("Creating DEEP fields...")
        BANDS = ['g','i']
        EXPTIME = [300,300]
        TILINGS = [15,10]

        dirname = '/Users/kadrlica/delve/observing/data'
        hexbase = 100000 # hex offset
        # target number and filename
        basenames = odict([
            (000, 'sextansB_fields.txt'),
            (100, 'ic5152_fields.txt'),
            (200, 'ngc300_fields.txt'),
            (300, 'ngc55_fields.txt')
        ])

        fields = FieldArray()
        for num,basename in basenames.items():
            filename = os.path.join(dirname,basename)
            target = np.genfromtxt(filename,names=True,dtype=None)
            f = FieldArray(len(target))

            filters = np.tile(BANDS,len(f))
            exptimes = np.tile(EXPTIME,len(f))

            f['RA']     = target['ra']
            f['DEC']    = target['dec']
            f['HEX']    = hexbase + num + target['field']
            f['TILING'] = target['tiling']
            f['PRIORITY'] = num + target['priority']
            f['PROGRAM'] = PROGRAM+'-deep'

            # Group the fields by hex/tiling
            f = np.repeat(f,len(BANDS))
            f['FILTER'] = filters
            f['EXPTIME'] = exptimes

            for (b,t) in zip(BANDS, TILINGS):
                f = f[~((f['FILTER'] == b) & (f['TILING'] > t))]

            fields = fields + f

        nhexes = len(np.unique(fields['HEX']))
        logging.info("  Number of hexes: %d"%nhexes)
        logging.info("  Filters: %s"%BANDS)
        logging.info("  Exposure time: %s"%EXPTIME)

        if plot: self.plot_depth(fields,depth,'delve-deep-%s-gt%i.png')

        logging.info("Number of target fields: %d"%len(fields))

        outfile = 'delve-deep-fields.fits.fz'
        logging.info("Writing %s..."%outfile)
        fields.write(outfile,clobber=True)

        return fields

    def create_deep_fields2(self, data, plot=False):
        """ DEPRECATED: Create the deep field observations """

        logging.info("Creating DEEP fields...")
        BANDS = ['g','i']
        EXPTIME = [300,300]
        TILINGS = [15,10]

        d = data[data['PASS'] == 1]
        sel = self.footprintDEEP(d['RA'],d['DEC'])
        data = np.copy(d[sel])

        nhexes = len(np.unique(data['TILEID']))
        ntilings = np.sum(TILINGS)
        nbands = len(BANDS)

        nfields = np.sum( np.sum(sel) * np.array(TILINGS))

        logging.info("  Number of hexes: %d"%nhexes)
        logging.info("  Filters: %s"%BANDS)
        logging.info("  Exposure time: %s"%EXPTIME)
        logging.info("  Tilings: %s"%TILINGS)

        tilings = np.array(range(1,TILINGS[0]+1)+range(1,TILINGS[1]+1))
        filters = np.repeat(BANDS,TILINGS)
        exptimes = np.repeat(EXPTIME,TILINGS)

        fields = FieldArray(nfields)
        fields['PROGRAM'] = PROGRAM+'-deep'
        fields['HEX'] = np.repeat(data['TILEID'],ntilings)
        fields['RA'] = np.repeat(data['RA'],ntilings)
        fields['DEC'] = np.repeat(data['DEC'],ntilings)

        fields['EXPTIME'] = np.tile(exptimes,nhexes)
        fields['TILING'] = np.tile(tilings,nhexes)
        fields['FILTER'] = np.tile(filters,nhexes)
        fields['PRIORITY'] = fields['TILING']

        frac, depth = self.covered(fields)
        fields['PRIORITY'][depth > fields['TILING']*fields['EXPTIME']] = DONE

        if plot: self.plot_depth(fields,depth,'delve-deep-%s-gt%i.png')

        logging.info("Number of target fields: %d"%len(fields))

        outfile = 'delve-deep-fields.fits.fz'
        logging.info("Writing %s..."%outfile)
        fields.write(outfile,clobber=True)

        return fields


    @staticmethod
    def footprintDEEP(ra,dec):
        """ Selecting exposures around the deep drilling fields """
        ra,dec = np.copy(ra), np.copy(dec)
        sel = np.zeros(len(ra),dtype=bool)
        filename = fileio.get_datafile('LV_MC_analogs_DECam.txt')
        targets = np.genfromtxt(filename,names=True,dtype=None)
        for t in targets:
            sel |= (angsep(t['RA'],t['Dec'],ra,dec) < t['r_vir'])
        return sel

    @staticmethod
    def footprintMC(ra,dec):
        """ Selecting exposures around the Magellanic Clouds """
        ra,dec = np.copy(ra), np.copy(dec)
        sel = angsep(constants.RA_LMC,constants.DEC_LMC,ra,dec) < 25.0
        sel |= angsep(constants.RA_SMC,constants.DEC_SMC,ra,dec) < 15.0

        return sel

    @staticmethod
    def footprintWIDE(ra,dec):
        """ Selecting wide-field exposures plane """
        ra,dec = np.copy(ra), np.copy(dec)
        sel = (dec < 0)
        return sel

    @staticmethod
    def footprintMaglites2(ra,dec):
        from obztak.maglites2 import Maglites2Survey
        return Maglites2Survey.footprint(ra,dec)

    @staticmethod
    def covered(fields, percent=67.):
        """
        Determine which fields haven't been previously covered by DECam

        Parameters:
        -----------
        fields : field information
        percent: fraction of the field that is covered

        Returns:
        --------
        sel, frac : selection of fields and coverage fraction
        """
        import healpy as hp
        # These maps are SUM(teff * exptime)
        dirname = '/Users/kadrlica/delve/observing/data'
        basename = 'decam_sum_expmap_%s_n1024.fits.gz'

        sel = np.ones(len(fields),dtype=bool)
        frac  = np.zeros(len(fields),dtype=float)
        depth = np.zeros(len(fields),dtype=float)
        ra,dec,band=fields['RA'],fields['DEC'],fields['FILTER']

        for b in np.unique(band):
            idx = (band==b)
            filename = os.path.join(dirname,basename%b)
            logging.info("Reading %s..."%os.path.basename(filename))
            skymap = hp.read_map(filename,verbose=False)

            nside = hp.get_nside(skymap)
            vec = hp.ang2vec(np.radians(90.-dec[idx]),np.radians(ra[idx]))

            f,d = [],[]
            for i,v in enumerate(vec):
                print '\r%s/%s'%(i+1,len(vec)),
                sys.stdout.flush()
                pix = hp.query_disc(nside,v,np.radians(constants.DECAM))
                # Find the 33rd percentile sum of effective exposure time
                # i.e., 67% of the pixels have a larger SUM(teff * exptime)
                d.append( np.percentile(skymap[pix],100-percent))
                # Find the detection fraction
                f.append((skymap[pix] > d[-1]).sum()/float(len(pix)))

            print
            frac[idx] = np.array(f)
            depth[idx] = np.array(d)

        return frac,depth

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


class DelveFieldArray(FieldArray):
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
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW, PROGRAM
        FROM exposure where propid = '%(propid)s' and exptime > 89
        and discard = False and delivered = True and flavor = 'object'
        and object like '%(object_fmt)s%%'
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class DelveScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',fileio.get_datafile("delve-windows.csv.gz")),
        ('targets',fileio.get_datafile("delve-target-fields.csv.gz")),
    ])

    FieldType = DelveFieldArray


class DelveTactician(Tactician):
    CONDITIONS = odict([
        (None,       [1.0, 2.0]),
        ('wide',     [1.0, 1.4]),
        ('deep',     [1.0, 1.4]),
        ('mc',       [1.0, 2.0]),
    ])

    def __init__(self, *args, **kwargs):
        super(DelveTactician,self).__init__(*args,**kwargs)
        self.mode = kwargs.get('mode','wide')

    @property
    def viable_fields(self):
        viable = super(DelveTactician,self).viable_fields
        viable &= (self.fields['PRIORITY'] >= 0)
        return viable

    def skybright_select(self):
        """Select fields based on skybrightness and band.

        Parameters:
        -----------
        None

        Returns:
        --------
        sel : boolean selection
        """
        sel = np.ones(len(self.fields),dtype=bool)

        if (self.sun.alt > -0.28):
            # i-band if Sun altitude > -16 deg
            sel &= (np.char.count('i',self.fields['FILTER']) > 0)
        # Moon band constraints (alt = 0.175 rad = 10 deg)
        elif (self.moon.phase >= 80) and (self.moon.alt > 0.175):
            # Moon is very bright; only do i
            sel &= (np.char.count('i',self.fields['FILTER']) > 0)
            # Allow i,z but prefer z
            #sel &= (np.char.count('iz',self.fields['FILTER']) > 0)
            #weight += 1e2 * (np.char.count('i',self.fields['FILTER']) > 0)
        #elif (self.moon.phase >= 45) and (self.moon.alt > 0.175):
        elif (self.moon.phase >= 45) and (self.moon.alt > 0.0):
            # Moon is more than half full; do r,i
            sel &= (np.char.count('ri',self.fields['FILTER']) > 0)
        else:
            # Moon is faint or down; do g,r (unless none available)
            sel &= (np.char.count('gr',self.fields['FILTER']) > 0)
            #weight += 1e8 * (np.char.count('iz',self.fields['FILTER']) > 0)
        return sel

    @property
    def weight(self):

        if self.mode is None:
            # First priority is deep
            weights = self.weight_deep()
            if np.isfinite(weights).sum(): return weights
            # Then mc
            weights = self.weight_mc()
            if np.isfinite(weights).sum(): return weights
            # Then wide
            weights = self.weight_wide()
            if np.isfinite(weights).sum(): return weights
        elif self.mode == 'deep':
            return self.weight_deep()
        elif self.mode == 'mc':
            return self.weight_mc()
        elif self.mode == 'wide':
            return self.weight_wide()
        else:
            raise ValueError("Unrecognized mode: %s"%self.mode)

        raise ValueError("No viable fields")

    def weight_wide(self):
        """ Calculate the field weight for the WIDE survey. """
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        sel &= (self.fields['PROGRAM'] == 'delve-wide')
        weight = np.zeros(len(sel))

        # Moon angle constraints
        moon_limit = 50.
        sel &= (moon_angle > moon_limit)

        # Sky brightness selection
        sel &= self.skybright_select()

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS['wide']
        sel &= ((airmass > airmass_min) & (airmass < airmass_max))

        # Higher weight for rising fields (higher hour angle)
        # HA [min,max] = [-53,54] (for airmass 1.4)
        #weight += 5.0 * self.hour_angle
        #weight += 1.0 * self.hour_angle
        weight += 0.1 * self.hour_angle

        # Higher weight for larger slews
        # slew = 10 deg -> weight = 1e2
        weight += self.slew**2
        #weight += self.slew
        #weight += 1e3 * self.slew

        # Higher weight for higher airmass
        # airmass = 1.4 -> weight = 6.4
        #weight += 100. * (airmass - 1.)**3
        weight += 1e3 * (airmass - 1.)**2

        # Higher weight for fields close to the moon (when up)
        # angle = 50 -> weight = 6.4
        if (self.moon.alt > -0.04) and (self.moon.phase >= 10):
            #weight += 100 * (35./moon_angle)**3
            #weight += 10 * (35./moon_angle)**3
            weight += 1 * (40./moon_angle)**3

        ## Try hard to do high priority fields
        weight += 1e3 * (self.fields['PRIORITY'] - 1)

        # Set infinite weight to all disallowed fields
        weight[~sel] = np.inf

        return weight

    def weight_deep(self):
        """ Calculate the field weight for the WIDE survey. """
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        sel &= (self.fields['PROGRAM'] == 'delve-deep')
        weight = np.zeros(len(sel))

        # Moon angle constraints
        moon_limit = 30.
        sel &= (moon_angle > moon_limit)

        # Sky brightness selection
        sel &= self.skybright_select()

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS['deep']
        sel &= ((airmass > airmass_min) & (airmass < airmass_max))

        ## Try hard to do high priority fields
        weight += 1e2 * self.fields['PRIORITY']

        # Set infinite weight to all disallowed fields
        weight[~sel] = np.inf

        return weight

    def weight_mc(self):
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        sel &= (self.fields['PROGRAM'] == 'delve-mc')
        weight = np.zeros(len(sel))

        # Moon angle constraints
        moon_limit = 30.
        sel &= (moon_angle > moon_limit)

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS['wide']
        sel &= ((airmass > airmass_min) & (airmass < airmass_max))

        # Sky brightness selection
        sel &= self.skybright_select()

        # Only a single tiling
        #sel &= (self.fields['PRIORITY'] == 3)

        weight = 2.0 * self.hour_angle

        # Prioritize fields
        weight += 3. * 360. * self.fields['PRIORITY']

        # Slew weighting
        weight += self.slew**3
        # Try hard to do the same field
        weight += 1e5 * (self.slew != 0)

        # Higher weight for higher airmass
        # airmass = 1.4 -> weight = 6.4
        weight += 100. * (airmass - 1.)**3

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

