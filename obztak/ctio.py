#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
from distutils.version import LooseVersion

import numpy as np

import ephem
import pytz

from obztak.utils import constants, fileio
from obztak.utils import projector as proj

# WARNING: copy.deepcopy doesn't work for ephem.Observer
# need to use ephem.Observer.copy

class CTIO(ephem.Observer):
    """ Utility class for defining CTIO """
    def __init__(self):
        super(ephem.Observer,self).__init__()
        self.lon = constants.LON_CTIO
        self.lat = constants.LAT_CTIO
        self.elevation = constants.ELEVATION_CTIO
        self.tz = pytz.timezone('Chile/Continental')
        self.twilight = '-14'
        self._load_constraints()

    def get_lon(self):
        """ Return longitude in degrees. """
        if LooseVersion(ephem.__version__) > LooseVersion('4'):
            return np.degrees(self.lon)
        else:
            return self.lon

    def get_lat(self):
        """ Return latitude in degrees. """
        if LooseVersion(ephem.__version__) > LooseVersion('4'):
            return np.degrees(self.lat)
        else:
            return self.lat


    def utc2local(self, utc):
        """ Convert ephem in UTC to local observatory time """
        utc_dt = pytz.utc.localize(ephem.Date(utc).datetime())
        return utc_dt.astimezone(self.tz)

    def local2utc(self, local):
        """ Convert string of local time to UTC """
        local_dt = self.tz.localize(ephem.Date(local).datetime())
        return utc_dt.astimezone(self.tz)

    def next_sunrise(self, date=None):
        obs = self.copy()
        if date: obs.date = ephem.Date(date)
        obs.horizon = self.twilight
        return obs.next_rising(ephem.Sun(), use_center=True)

    def next_sunset(self, date=None):
        obs = self.copy()
        if date: obs.date = ephem.Date(date)
        obs.horizon = self.twilight
        return obs.next_setting(ephem.Sun(), use_center=True)

    def previous_sunrise(self, date=None):
        obs = self.copy()
        if date: obs.date = ephem.Date(date)
        obs.horizon = self.twilight
        return obs.next_rising(ephem.Sun(), use_center=True)

    def previous_sunset(self, date=None):
        obs = self.copy()
        if date: obs.date = ephem.Date(date)
        obs.horizon = self.twilight
        return obs.next_setting(ephem.Sun(), use_center=True)

    def _load_constraints(self, filename=None):
        """ Load Blanco constraint data """
        if filename is None:
            filename = os.path.join(fileio.get_datadir(),'blanco_hour_angle_limits.dat')
        self.constraints = np.genfromtxt(filename, names=True, dtype=None,
                                         encoding=None)

        # ADW: This is not very pythonic....
        ha_degrees = np.array([proj.hms2dec(ha) for ha in self.constraints['HA']])

        # Buffer to protect us from the chicken
        ha_degrees -= 1.25 
        
        self.ha_degrees = ha_degrees
        return self.ha_degrees

    def hour_angle_limit(self, dec):
        # Avoid relying on scipy.interpolate
        return np.interp(dec,self.constraints['Dec'],self.ha_degrees,left=-1,right=-1)

    def airmass_limit(self, dec):
        # Avoid relying on scipy.interpolate
        return np.interp(dec,self.constraints['Dec'], self.constraints['AirmassLimit'],left=-1,right=-1)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
