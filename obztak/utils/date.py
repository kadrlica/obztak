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

    # This is a bit annoying, but works better than converting to datetime
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

def nitestring(nite,sep=''):
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

def nite2utc(nite, observer=None):
    """ Convert a nite string to an ephem.Date of the UTC at sunset. 

    Parameters:
    -----------
    nite     : The nite (string or datetime)
    observer : An ephem.Observer object (defaults to CTIO)

    Returns:
    --------
    utc      : The UTC of sunset (ephem.Date)
    """
    
    import dateutil.parser
    import dateutil.tz
    import datetime

    if observer is None: observer = CTIO()
    else:                observer = observer.copy()
    observer.horizon = '-14'

    if not isinstance(nite,datetime.datetime):
        nite = dateutil.parser.parse(str(nite))

    # Move to (local) noon
    # This depends on where the user is and isn't very safe
    nite = nite.replace(hour=12,tzinfo=dateutil.tz.tzlocal())
    utc = ephem.Date(nite - nite.utcoffset())

    # Maybe something like this instead...
    #offset = int( (observer.lon/(2*np.pi)) * 24. * 60) * 60
    #tzinfo = dateutil.tz.tzoffset('OBS',offset)
    #nite = nite.replace(hour=12,tzinfo=tzinfo)
    #utc = ephem.Date(nite - nite.utcoffset())

    observer.date = utc
    #ret = observer.next_antitransit(ephem.Sun())
    ret = observer.next_setting(ephem.Sun(), use_center=True)

    return ret

def utc2nite(utc, observer=None):
    """ Convert the input time (UTC) to the UTC of the start of the
    corresponding nite. This will return the UTC of the previous
    sunset `utc` is at night, or the UTC of the next sunset of `utc`
    is during the day.

    Parameters:
    -----------
    utc     : The time to determine sunset for (ephem.Date)
    observer: The observer location (defaults to CTIO)
    
    Returns:
    --------
    nite    : The UTC of sunset on the corresponding nite (ephem.Date)
    """
    sun = ephem.Sun()

    if observer is None: observer = CTIO()
    else:                observer = observer.copy()
    observer.date = utc
    observer.horizon = '-14'

    if observer.previous_setting(sun) > observer.previous_rising(sun):
        # It's night time, use the date of the previous setting
        nite = ephem.localtime(observer.previous_setting(sun, use_center=True))
    else:
        # It's daytime, use the next setting
        nite = ephem.localtime(observer.next_setting(sun, use_center=True))

    return ephem.Date(nite)

def get_nite(utc=None, observer=None):
    """Convert from a UTC date and time to the 'nite'.

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
