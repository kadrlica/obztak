#!/usr/bin/env python
"""
Generic python script.
"""
import os
import glob
import json

import numpy as np
import ephem

from obztak.utils.testing import call, check_dict, make_options
import obztak.tactician
from obztak.utils.date import datestring
from obztak.field import FieldArray
from obztak.utils import fileio

def create_fields():
    nfields = 10
    fields = FieldArray(nfields)
    fields['RA'] = np.linspace(3,359,nfields)
    fields['DEC'] = np.linspace(-89,89,nfields)
    fields['TILING'] = np.array([1,1,1,2,2,2,3,3,4,4])
    return fields

def create_tactician():
    fields = create_fields()
    tac = obztak.tactician.Tactician(fields)
    tac.set_date('2017/02/08 04:00:00')
    return tac

def test_create_tactician():
    tac = create_tactician()

def test_calculate_zenith():
    tac = create_tactician()
    ra,dec = tac.zenith_angle
    np.testing.assert_almost_equal(ra,127.46761160249451)
    np.testing.assert_almost_equal(dec,-30.108338409778298)

def test_calculate_airmass():
    tac = create_tactician()

    #test_airmass = np.array([2.0283324,  2.0152044,  1.2880471, 1.0038942,
    #                         1.2583179,  6.3051279,  999.     , 999.     ,
    #                         999.      , 999.       ])
    test_airmass = np.array([2.0283324,  2.0152044,  1.2880471, 1.0038942,
                             1.2583179,  6.3051279,  np.inf   , np.inf   ,
                             np.inf   ,  np.inf      ])
    np.testing.assert_almost_equal(tac.airmass,test_airmass,
                                   err_msg='airmass')

    #test_airmass_next = np.array([1.9934748, 4.77775207, 2.97209942, 999.,
    #                              999., 999., 2.86569785, 1.63423598, 999.,
    #                              999.])
    #np.testing.assert_almost_equal(tac.airmass_next,test_airmass_next,
    #                               err_msg='airmass_next')

def test_calculate_moon():
    tac = create_tactician()

    #test_moon_angle = np.array([108.77523275, 81.48432138, 22.26324414, 43.01144318, 82.20330354, 76.6695198, 63.40331764, 98.01411913, 136.98408244, 103.7186223 ])

    # Without topocentric correction (http://rhodesmill.org/pyephem/radec)
    test_moon_angle = np.array([108.9487001, 97.8196126, 70.7322332, 51.6236834, 64.0174593, 93.9227775, 116.3650182, 111.7185502, 88.2697742, 71.4721274])
    # With topocentric correction
    test_moon_angle = np.array([109.6670211, 98.396651, 71.3473021, 52.4287292, 64.6183409, 94.124153, 116.0532799, 111.0075416, 87.4853729, 70.7408394])
    # Newer version of ephem (2023)
    test_moon_angle = np.array([109.6670216,  98.3966407,  71.3472944,  52.4287432,  64.6183752, 94.1241885, 116.0533029, 111.0075417,  87.4853641,  70.7408376])

    np.testing.assert_almost_equal(tac.moon_angle,test_moon_angle,decimal=5,
                                   err_msg='moon_angle')

    test_moon_phase = 89.82009887695312
    # Newer version of ephem (2023)
    test_moon_phase = 89.82008361816406
    np.testing.assert_almost_equal(tac.moon_phase,test_moon_phase, decimal=4,
                                   err_msg='moon_phase')

def test_calculate_slew():
    tac = create_tactician()

    # Use a previously completed field
    completed_fields = FieldArray(1)
    completed_fields[0]['RA'] = 137
    completed_fields[0]['DEC'] = -43
    completed_fields[0]['DATE'] = datestring(tac.date)

    #test_slew = np.array([47., 66.6686451, 72.6195061, 103.43996017, 149.98635165, 122.82243022, 61.24007137, 38.8279613, 86.08700171, 122.10797116])
    test_slew = np.array([47.6988041, 51.8633918, 37.6586078,18.1169301, 39.2690875, 78.3913867, 118.9542319, 153.814957, 153.7471813, 133.7394395])
    zero_slew = np.zeros(len(tac.fields))

    tac.set_completed_fields(completed_fields)
    np.testing.assert_almost_equal(tac.slew, test_slew,
                            err_msg='slew from previous field')

    # Set previous field to be < 30 min before date
    completed_fields[0]['DATE'] = datestring(tac.date - 29*ephem.minute)
    tac.set_completed_fields(completed_fields)
    np.testing.assert_almost_equal(tac.slew, test_slew,
                            err_msg='slew not calculated for <30 min gap')

    # Set previous field to be > 30 min before date
    completed_fields[0]['DATE'] = datestring(tac.date - 31*ephem.minute)
    tac.set_completed_fields(completed_fields)
    np.testing.assert_equal(tac.slew, zero_slew,
                            err_msg='non-zero slew for >30 min gap')

    # Set previous field to None
    tac.set_completed_fields(None)
    np.testing.assert_equal(tac.slew, zero_slew,
                            err_msg='non-zero slew without previous field')

def test_calculate_hour_angle():
    tac = create_tactician()

    test_hour_angle = np.array([-124.4676116, -84.91205605, -45.35650049, -5.80094494, 33.75461062, 73.31016618, 112.86572173, 152.42127729, -168.02316716, -128.4676116])
    np.testing.assert_almost_equal(tac.hour_angle,test_hour_angle)

    fields = create_fields()
    np.testing.assert_equal(tac.fields,fields,
                            err_msg ="fields have been altered")

def test_viable_fields():
    tac = create_tactician()
    sel = tac.viable_fields

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
