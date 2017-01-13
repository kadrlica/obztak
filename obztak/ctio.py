#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import ephem
import pytz

import numpy as np

from obztak.utils import constants
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
            from obztak.utils import fileio
            filename = os.path.join(fileio.get_datadir(),'blanco_hour_angle_limits.dat')
        self.constraints = np.recfromtxt(filename, names=True)

        # ADW: This is not very pythonic....
        ha_degrees = np.tile(0., len(self.constraints['HA']))
        for ii in range(0, len(self.constraints['HA'])):
            ha_degrees[ii] = proj.hms2dec(self.constraints['HA'][ii])

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
