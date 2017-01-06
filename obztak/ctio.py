#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import ephem
import pytz

from obztak.utils import constants

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

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
