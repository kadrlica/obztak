#!/usr/bin/env python
"""
Date and time utilities.
"""

import numpy as np
import ephem
import time
import logging

import obztak.utils.constants as constants
from obztak.ctio import CTIO

def setdefaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

def datestring(date,precision=4):
    """
    Convert an ephem.Date object to a string with increased precision

    Parameters:
    -----------
    date      : ephem.Date object
    precision : Output precision

    Returns:
    --------
    datestr   : String representation of the date
    """
    """
    date = ephem.Date(date).datetime()
    datestr = date.strftime('%Y/%m/%d %H:%M:%S')
    datestr += '{:.{precision}f}'.format(date.microsecond/1.e6,
                                         precision=precision)[1:]
    """
    if precision < 0:
        msg = "Precision must be positive."
        raise Exception(msg)

    # This is a bit annoying, but works better (more accurate) than
    # using the built-in datetime conversion
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],precision)
    minutes = datetuple[-2]
    hours = datetuple[-3]
    minutes += int(seconds//60)
    seconds = seconds%60.
    hours += int(minutes//60)
    minutes = minutes%60

    strtuple = datetuple[:-3]+(hours,minutes,seconds)
    width = precision+2 if precision == 0 else precision+3

    #datestr = '%d/%02d/%02d %02i:%02i:%07.4f'%strtuple
    datestr = '{:4d}/{:02d}/{:02d} {:02d}:{:02d}:{:0{width}.{precision}f}'
    datestr = datestr.format(*strtuple,precision=precision,width=width)

    return datestr

def datestr(date,precision=0): return datestring(date,precision)

def nitestr(nite,sep=''):
    """
    Convert an ephem.Date object to a nite string.

    Parameters:
    -----------
    nite     : ephem.Date object
    sep      : Output separator

    Returns:
    --------
    nitestr  : String representation of the nite
    """
    import dateutil.parser
    if isinstance(nite,basestring):
        nite = dateutil.parser.parse(nite)
    nite = ephem.Date(nite)
    strtuple = nite.tuple()[:3]
    nitestr = '{:4d}{sep}{:02d}{sep}{:02d}'
    nitestr = nitestr.format(*strtuple,sep=sep)
    return nitestr

nitestring = nitestr

def nite2utc(nite, observer=None):
    """ Convert a nite string to an ephem.Date of the UTC at sunset.

    Parameters:
    -----------
    nite     : The nite (string or datetime)
    observer : An ephem.Observer object (defaults to CTIO)

    Returns:
    --------
    utc      : The UTC at sunset on this nite (ephem.Date)
    """
    import datetime, pytz, dateutil.parser

    if observer is None: observer = CTIO()
    else:                observer = observer.copy()
    observer.horizon = observer.twilight

    # Parse the nite string (or ephem.Date)
    if not isinstance(nite,datetime.datetime):
        nite = dateutil.parser.parse(str(nite))

    # Find local noon at the observer
    noon = observer.tz.localize(nite).replace(hour=12)
    utc_noon = ephem.Date(noon.astimezone(pytz.utc))

    observer.date = utc_noon
    sunset = observer.next_setting(ephem.Sun(), use_center=True)
    return sunset

def utc2nite(utc, observer=None):
    """Convert the input time (UTC) to a corresponding nite string.
    This derives the nite from the local date of the previous sunrise.

    Parameters:
    -----------
    utc     : The input UTC time (ephem.Date)
    observer: The observer location (defaults to CTIO)

    Returns:
    --------
    nite    : The nite string

    """
    import pytz

    if observer is None: observer = CTIO()
    else:                observer = observer.copy()

    observer.date = ephem.Date(utc)
    utc_sunrise = observer.previous_rising(ephem.Sun(), use_center=True)
    utc_sunrise = pytz.utc.localize(utc_sunrise.datetime())
    sunrise = utc_sunrise.astimezone(observer.tz)
    return nitestring(sunrise)


def get_nite(utc=None, observer=None):
    """ ADW 20170205: IS THIS DEPRICATED???

    Convert from a UTC date and time to the 'nite'.

    A 'nite' is defined by the day (UTC) at noon local time in Chile
    before observing started. This follows the usual convention of
    naming the nite after the day (local time) on which it starts.

    Parameters:
    -----------
    utc : The date/time (UTC) to calculate the nite from.

    Returns:
    --------
    nite : An ephem.Date object containing the nite (at sunset)

    """
    if not utc: utc = ephem.now()
    return utc2nite(utc, observer)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
