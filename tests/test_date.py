#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

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
    for val1,val2 in values:
        utc = nite2utc(val1)
        test = ephem.Date(val2)
        np.testing.assert_almost_equal(utc,test,5)

def test_utc2nite():
    values = [
        ('2016/02/14 18:13:00','2016/2/14 18:36:51'),
        ('2016/06/29 00:10:00','2016/6/28 17:59:17'),
        ('2017/02/23 06:00:01','2017/2/22 18:27:05'),
        ]
    for val1,val2 in values:
        nite = utc2nite(val1)
        test = ephem.Date(val2)
        np.testing.assert_almost_equal(nite,test,5)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
