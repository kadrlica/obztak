#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
from collections import OrderedDict as odict
import logging
import copy

import numpy as np
import pandas as pd
import dateutil.parser
import ephem

from obztak.utils import fileio
from obztak.utils.date import datestring
from obztak.utils.database import Database


# These are nominal transformation values from Eric Neilsen
# WAVE[x] = (lambda[x]/lambda[i])**0.2
WAVE = odict([
    ( 'u'  , 0.86603  ), # u (380nm) -> i (780nm)
    ( 'g'  , 0.9067   ), # g (480nm) -> i (780nm)
    ( 'r'  , 0.9609   ), # r (640nm) -> i (780nm)
    ( 'i'  ,    1.0   ), # i (780nm) -> i (780nm)
    ( 'z'  ,  1.036   ), # z (920nm) -> i (780nm)
    ( 'Y'  , 1.0523   ), # Y (990nm) -> i (780nm)
    ('dimm', 1/1.0916 ), # dimm (500 nm)->i (780nm)
    ('VR'  , 0.9551   ), # VR (620 nm)->i (780nm)
])
WAVE_DF = pd.DataFrame({'filter':WAVE.keys(),'trans':WAVE.values()})
DECAMINST = 0.5 # DECam instrumental contribution to the PSF [arcsec]
DIMMINST = 0.0  # DIMM instrumental contribution to the PSF [arcsec]

def convert(fwhm_1,
            band_1='dimm', airmass_1=1.0, inst_1=DIMMINST,
            band_2='i',    airmass_2=1.0, inst_2=DECAMINST):

    """
    Convert observed seeing value to another band and airmass.

    Parameters:
    -----------
    fwhm_1   : observed fwhm [arcsec]
    band_1   : observed band ['g','r','i','z','Y','dimm']
    airmass_1: observed airmass
    inst_1   : instrumental contribution to the observed psf [arcsec]

    band_2    : output band ['g','r','i','z','Y','dimm']
    airmass_2 : output airmass
    inst_2    : instrumental contribution to the output psf [arcsec]

    Returns:
    --------
    fwhm_2    : output fwhm [arcsec]
    """
    fwhm = np.sqrt(fwhm_1**2 - inst_1**2)
    if np.isscalar(band_1):
        wave_1 = WAVE[band_1]
    else:
        wave_1 = WAVE_DF.merge(pd.DataFrame({'filter':band_1}), on='filter').to_records()['trans']

    if np.isscalar(band_2):
        wave_2 = WAVE[band_2]
    else:
        wave_2 = WAVE_DF.merge(pd.DataFrame({'filter':band_2}), on='filter').to_records()['trans']

    fwhm_2 = fwhm * (wave_1/wave_2) * (airmass_2/airmass_1)**(0.6)

    return np.hypot(fwhm_2, inst_2)

class Seeing():
    """Class to manage seeing data. Seeign data is stored in two member variables:

    self.raw : the raw data before transformation
    self.data: seeing data transformed atmospheric i-band zenith

    The two values differ in that self.raw can have any source and
    includes the instrumental contribution. In contrast, self.data is
    the "atmospheric" i-band FWHM (arcsec). To get a prediction of the
    observed PSF, use `get_fwhm`.
    """
    DTYPE = [('date','<M8[ns]'),('fwhm',float),('airmass',float),('filter','S4')]

    def __init__(self, date=None, db='fnal', filename=None):
        self.set_date(date)
        self.df = self.read_file(filename)
        self.db = 'db-'+db

    def set_date(self, date):
        if date is None:
            #NOOP (consistent with Tactician)
            return
        elif date == 'now':
            self.date = dateutil.parser.parse(datestring(ephem.now()))
        else:
            self.date = dateutil.parser.parse(date)


    def get_fwhm(self, timedelta='10m', band='i', airmass=1.0, inst=DECAMINST):
        """Calculate the predict PSF FWHM (arcsec).

        Parameters:
        -----------
        date      : date to estimate the psf (defualt: now)
        timedelta : time range to use to estimate the psf
        band      : output band
        airmass   : output airmass
        inst      : output instrument contribution

        Returns:
        --------
        fwhm      : predicted fwhm (arcsec)
        """
        timedelta = pd.Timedelta(timedelta)

        self.load_data(timedelta=max(3*timedelta,pd.Timedelta('1h')))

        dt = pd.DatetimeIndex(self.data['date'])
        previous = slice(-1,None) # most recent exposure
        recent = (dt < self.date) & (dt > (self.date - timedelta))
        ancient = (dt < (self.date - timedelta)) & (dt > (self.date - 2*timedelta))

        # Nominal atmospheric psf i-band zenith fwhm = 0.9"
        xmu = np.log10(0.74833) # sqrt(0.9**2 - 0.5**2)

        if not len(self.data):
            # No data, use the mean and print a warning
            logging.warn("No fwhm data available; using median fwhm")
            xpred = xmu
        elif not len(recent):
            # Log of the i-band zenith fwhm from the previous exposure
            xpred = np.log10(self.data[previous]['fwhm'])
        elif not len(ancient):
            # Median of the log of the observed atmospheric psf i-band zenith
            xpred = np.log10(np.median(self.data[recent]['fwhm']))
        else:
            # Weighted median of recent and ancient exposures
            # Log of the observed atmospheric psf i-band zenith
            x = np.log10([np.median(self.data[recent]['fwhm']),
                          np.median(self.data[ancient]['fwhm'])])

            # Predicted log of the atmospheric psf
            # NB: These constants were derived for timedelta=5min
            # they may not hold for arbitrary time windows.
            xpred = xmu + 0.8 * (x[0] - xmu) + 0.14 * (x[1] - xmu)

        fwhm_pred = convert(10**xpred,
                            band_1='i' , airmass_1=1.0    , inst_1=0.0,
                            band_2=band, airmass_2=airmass, inst_2=inst)

        return fwhm_pred


class DimmSeeing(Seeing):
    """Estimate seeing from the DIMM."""

    @classmethod
    def read_file(cls, filename):
        if filename is None: return None
        df = pd.read_csv(filename,names=['date','fwhm'],
                         parse_dates=['date'],index_col=['date'])
        return df

    def get_data(self, date=None, timedelta='30m'):
        self.set_date(date)
        tmax = self.date
        tmin = self.date - pd.Timedelta(timedelta)
        if self.df is None:
            # Don't want to create the DB each time?
            db = Database(self.db)
            db.connect()
            query ="""
            select date, dimm2see as fwhm from exposure
            where date > '%s' and date < '%s'
            and dimm2see is not NULL
            """%(tmin, tmax)
            logging.debug(query)
            raw = db.query2rec(query)
        else:
            sel = (self.df.index > tmin) & (self.df.index < tmax)
            raw = self.df[sel].to_records()
        return raw

    def load_data(self, date=None, timedelta='30m'):
        raw = self.get_data(date, timedelta)

        # Save the raw dimm values
        self.raw = np.recarray(len(raw),dtype=self.DTYPE)
        self.raw['date'] = raw['date']
        self.raw['fwhm'] = raw['fwhm']
        self.raw['airmass'] = 1.0
        self.raw['filter'] = 'dimm'

        # Convert to i-band zenith
        self.data = copy.deepcopy(self.raw)
        self.data['filter'] = 'i'
        self.data['airmass'] = 1.0

        kwargs = dict(band_1='dimm', inst_1=DIMMINST, airmass_1=self.raw['airmass'])
        kwargs.update(band_2='i',    inst_2=0.0     , airmass_2=self.data['airmass'])
        self.data['fwhm'] = convert(self.raw['fwhm'],**kwargs)

        return self.data

class QcSeeing(Seeing):
    """Estimate seeing from the DECam QC values."""

    @classmethod
    def read_file(cls, filename):
        if filename is None: return None
        df = pd.read_csv(filename,names=['date','fwhm','airmass','filter'],
                         parse_dates=['date'],index_col=['date'])
        return df

    def get_data(self, date=None, timedelta='30m'):
        self.set_date(date)
        tmax = self.date
        tmin = self.date - pd.Timedelta(timedelta)
        if self.df is None:
            # Don't want to create the DB each time?
            db = Database()
            db.connect()
            query ="""
            select date, qc_fwhm as fwhm, airmass, filter from exposure
            where date > '%s' and date < '%s'
            --and filter != 'VR' and qc_fwhm is not NULL
            and qc_fwhm is not NULL
            """%(tmin, tmax)
            logging.debug(query)
            raw = db.query2rec(query)
        else:
            sel = (self.df.index > tmin) & (self.df.index < tmax)
            raw = self.df[sel].to_records()

        return raw

    def load_data(self, date=None, timedelta='30m'):
        raw = self.get_data(date,timedelta)

        # Save the raw dimm values
        self.raw = np.recarray(len(raw),dtype=self.DTYPE)
        self.raw['date']    = raw['date']
        self.raw['fwhm']    = raw['fwhm']
        self.raw['airmass'] = raw['airmass']
        self.raw['filter']  = raw['filter']

        # Convert to i-band zenith
        self.data = copy.deepcopy(self.raw)
        self.data['filter'] = 'i'
        self.data['airmass'] = 1.0

        kwargs = dict(band_1=self.raw['filter'], inst_1=DECAMINST, airmass_1=self.raw['airmass'])
        kwargs.update(band_2='i',                inst_2=0.0      , airmass_2=self.data['airmass'])
        self.data['fwhm'] = convert(self.raw['fwhm'],**kwargs)

        return self.data

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
