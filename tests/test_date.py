#!/usr/bin/env python
"""
Testing for the date utils.
"""

import numpy as np
import ephem
import time
import datetime
import dateutil.parser

from obztak.utils.date import datestring,nitestring
from obztak.utils.date import nite2utc,utc2nite,get_nite

def test_datestring():
    string = '2016/02/14 01:30:00.0000'
    dt = dateutil.parser.parse(string)
    eph = ephem.Date(string)

    for d in [dt,eph,string]:
        np.testing.assert_equal(datestring(d),string)
        for i,j in [(4,None),(3,-1),(2,-2),(1,-3),(0,-5)]:
            np.testing.assert_equal(datestring(d,i),string[slice(j)])

def test_nitestring():
    string = '20160214'
    dt = dateutil.parser.parse(string)
    eph = ephem.Date(dt)
    np.testing.assert_equal(nitestring(string),string)
    np.testing.assert_equal(nitestring(dt),string)
    np.testing.assert_equal(nitestring(eph),string)

    np.testing.assert_equal(nitestring(eph,'/'),'2016/02/14')
    np.testing.assert_equal(nitestring(eph,'-'),'2016-02-14')

def test_nite2utc():
    values = [
        ('20160214','2016/2/15 00:36:52'),
        ('20160628','2016/6/28 22:59:18'),
        ('20170222','2017/2/23 00:27:05'),
        ]
    for string,test in values:
        dt = dateutil.parser.parse(val1)
        eph = ephem.Date(dt)

        np.testing.assert_almost_equal(nite2utc(string),ephem.Date(test)n,5)
        np.testing.assert_almost_equal(nite2utc(eph),ephem.Date(test),5)
        np.testing.assert_almost_equal(nite2utc(dt),ephem.Date(test),5)

def test_utc2nite():
    values = [
        ('2016/02/14 18:13:00','20160214'),
        ('2016/06/29 00:10:00','20160628'),
        ('2017/02/23 06:00:01','20170222'),
        ]
    for string,test in values:
        dt = dateutil.parser.parse(val1)
        eph = ephem.Date(val1)

        np.testing.assert_equal(utc2nite(string),test)
        np.testing.assert_equal(utc2nite(eph),test)
        np.testing.assert_equal(utc2nite(dt),test)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
