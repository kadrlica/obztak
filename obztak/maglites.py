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
from obztak.utils.date import datestring, setdefaults, nite2utc,utc2nite
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
            infile = fileio.get_datafile('smash_fields_alltiles.txt')
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
        -- Final observing push
        --and (
        --  -- All exposures from tonight
        --  (date > '2017/06/20 18:00:00')
        --  -- Or pass teff cut or NULL from past nights
        --  or ((qc_teff is NULL or qc_teff > 0.1) and date < '2017/06/20 18:00:00')
        --  -- Or teff = 0 from February
        --  or (qc_teff = 0 and date < '2017/03/01' and date > '2017/02/01')
        --)
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


class MaglitesScheduler(Scheduler):
    _defaults = odict(Scheduler._defaults.items() + [
        ('tactician','coverage'),
        ('windows',fileio.get_datafile("maglites-windows.csv")),
        ('targets',fileio.get_datafile("maglites-target-fields.csv")),
    ])

    FieldType = MaglitesFieldArray

# Some plotting functions

def plot_nightsum(fields,nitestr):
    """ Plot the night summary for MagLiteS.

    Parameters:
    -----------
    fields:  the fields observed tonight
    nitestr: the nite in strig format

    Returns:
    --------
    None
    """
    import pylab as plt
    from obztak.utils.database import Database
    from obztak.utils.ortho import makePlot

    #fields = FieldArray.load_database()
    #new = np.char.startswith(fields['DATE'],date)

    date = nite2utc(nitestr)
    new = (np.array(map(utc2nite,fields['DATE'])) == nitestr)
    new_fields = fields[new]
    old_fields = fields[~new]

    kwargs = dict(edgecolor='none', s=50, vmin=0, vmax=4)
    fig,basemap = makePlot(date=nitestr,name='nightsum',moon=False,airmass=False,center=(0,-90),bliss=False)
    plt.title('Coverage (%s)'%nitestr)
    kwargs['cmap'] = 'gray_r'
    proj = basemap.proj(old_fields['RA'], old_fields['DEC'])
    basemap.scatter(*proj, c=old_fields['TILING'],**kwargs)

    kwargs['cmap'] = 'summer_r'
    proj = basemap.proj(new_fields['RA'], new_fields['DEC'])
    basemap.scatter(*proj, c=new_fields['TILING'],  **kwargs)
    colorbar = plt.colorbar()
    colorbar.set_label('Tiling')

    plt.plot(np.nan, np.nan,'o',color='green',mec='green',label='Observed tonight')
    plt.plot(np.nan, np.nan,'o',color='0.7',mec='0.7',label='Observed previously')
    plt.legend(fontsize=10,loc='lower left',scatterpoints=1)
    plt.savefig('nightsum_coverage_%s.png'%nitestr,bbox_inches='tight')

    db = Database()
    db.connect()

    query = """
    select id, qc_fwhm as psf, qc_teff as teff from exposure
    where exptime = 90 and delivered = True 
    and propid = '%s'
    and qc_teff is not NULL and qc_fwhm is not NULL
    and to_timestamp(utc_beg) %s '%s'
    """

    new = db.query2recarray(query%(PROPID,'>',date))
    old = db.query2recarray(query%(PROPID,'<',date))

    nbins = 35
    kwargs = dict(normed=True)
    step_kwargs = dict(kwargs,histtype='step',lw=3.5)
    fill_kwargs = dict(kwargs,histtype='stepfilled',lw=1.0,alpha=0.7)

    plt.figure()
    step_kwargs['bins'] = np.linspace(0.5,2.5,nbins)
    fill_kwargs['bins'] = np.linspace(0.5,2.5,nbins)
    plt.hist(new['psf'],color='green',zorder=10, label='Observed tonight', **fill_kwargs)
    plt.hist(new['psf'],color='green',zorder=10, **step_kwargs)
    plt.hist(old['psf'],color='0.5', label='Observed previously', **fill_kwargs)
    plt.hist(old['psf'],color='0.5', **step_kwargs)
    plt.axvline(1.20,ls='--',lw=2,color='gray')
    plt.legend()
    plt.title('Seeing (%s)'%nitestr)
    plt.xlabel('FWHM (arcsec)')
    plt.ylabel('Normalized Number of Exposures')
    plt.savefig('nightsum_psf_%s.png'%nitestr,bbox_inches='tight')

    plt.figure()
    step_kwargs['bins'] = np.linspace(0,1.5,nbins)
    fill_kwargs['bins'] = np.linspace(0,1.5,nbins)
    plt.hist(new['teff'],color='green',zorder=10,label='Observed tonight', **fill_kwargs)
    plt.hist(new['teff'],color='green',zorder=10, **step_kwargs)
    plt.hist(old['teff'],color='0.5',label='Observed previously', **fill_kwargs)
    plt.hist(old['teff'],color='0.5', **step_kwargs)
    plt.axvline(0.25,ls='--',lw=2,color='gray')
    plt.legend()
    plt.title('Effective Depth (%s)'%nitestr)
    plt.xlabel('Teff')
    plt.ylabel('Normalized Number of Exposures')
    plt.savefig('nightsum_teff_%s.png'%nitestr,bbox_inches='tight')

def plot_progress(outfile=None,**kwargs):
    """ ADW 2018-07-22: DEPRECATED? """
    defaults = dict(edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    for k,v in defaults.items():
        kwargs.setdefault(k,v)

    fields = FieldArray.load_database()

    nites = [get_nite(date) for date in fields['DATE']]
    nite = ephem.Date(np.max(nites))
    date = '%d/%02d/%d 00:00:00'%(nite.tuple()[:3])

    fig,basemap = makePlot(date=date,moon=False,airmass=False,center=(0,-90),smash=False)
    proj = basemap.proj(fields['RA'],fields['DEC'])
    basemap.scatter(*proj, c=fields['TILING'],  **kwargs)
    colorbar = plt.colorbar()
    colorbar.set_label('Tiling')
    plt.title('Coverage (%d/%02d/%d)'%nite.tuple()[:3])

    if outfile is not None:
        plt.savefig(outfile,bbox_inches='tight')

    return fig,basemap

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
