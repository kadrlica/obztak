#!/usr/bin/env python
"""
Test script for field parsing.
"""
__author__ = "Alex Drlica-Wagner"

import subprocess
import json
import tempfile

import numpy as np

import obztak.field
from obztak.utils.testing import call, make_options

FILENAME='field_test.json'
open(FILENAME,'w').write("""
[
    {
        "object": "MAGLITES field: 6260-01-g",
        "seqnum": 1,
        "seqtot": 2,
        "seqid": "MAGLITES scheduled: 2017/02/22 06:00:00.0000",
        "expTime": 90.0000,
        "RA": 113.4230,
        "dec": -55.4093,
        "filter": "g",
        "count": 1,
        "expType": "object",
        "program": "maglites",
        "wait": "False",
        "propid": "2016A-0366",
        "comment": "MAGLITES v3.0.dev1+0.gf198e23.dirty: PRIORITY=1, AIRMASS=1.4844, SLEW=0.0000, MOONANGLE=105.2010, HOURANGLE=-57.8900"
    },
    {
        "object": "MAGLITES field: 6260-01-r",
        "seqnum": 2,
        "seqtot": 2,
        "seqid": "MAGLITES scheduled: 2017/02/22 06:02:00.0000",
        "expTime": 90.0000,
        "RA": 113.4230,
        "dec": -55.4093,
        "filter": "r",
        "count": 1,
        "expType": "object",
        "program": "maglites",
        "wait": "False",
        "propid": "2016A-0366",
        "comment": "MAGLITES v3.0.dev1+0.gf198e23.dirty: PRIORITY=1, AIRMASS=1.4844, SLEW=0.0000, MOONANGLE=105.2010, HOURANGLE=-57.8900"
    }
]
""")

def test_empty_fields():
    fields = obztak.field.FieldArray()

    assert len(fields) == 0

    np.testing.assert_equal(fields.dtype.names,['HEX','RA','DEC','FILTER','EXPTIME','TILING','PRIORITY','DATE','AIRMASS','SLEW','MOONANGLE','HOURANGLE','PROGRAM'])

def test_sispi_fields():
    sispi = json.loads(open(FILENAME).read())
    fields = obztak.field.FieldArray.read(FILENAME)

    np.testing.assert_equal(len(sispi),len(fields))
    for k in ['RA','dec','expTime']:
        np.testing.assert_allclose([s[k] for s in sispi],fields[k.upper()])
    for k in ['filter']:
        np.testing.assert_equal([s[k] for s in sispi],fields[k.upper()])

    assert np.all(np.char.count([s['object'] for s in sispi],fields.unique_id))
    np.testing.assert_equal([s['object'] for s in sispi],fields.object)

    assert np.all(np.char.count([s['seqid'] for s in sispi],fields['DATE']))
    np.testing.assert_equal([s['seqid'] for s in sispi],fields.seqid)

    sispi_comment = [str(s['comment']).split(':',1)[-1] for s in sispi]
    field_comment = np.char.partition(fields.comment,':')[:,-1]
    np.testing.assert_equal(sispi_comment,field_comment)

def test_database_fields():
    fields = obztak.field.FieldArray.load_database('db-fnal')

    # No database connection
    if len(fields) == 0: return

    idx = [0,17,357,1703]
    test = np.rec.array([
        (5143,101.331075,-55.618165,'g',90.0,1,1,'2016/02/11 04:56:51.972',1.29,87.4149,97.1638, 42.6655  , ''),
        (5345,106.011075,-58.41872, 'r',90.0,1,1,'2016/02/11 05:30:32.472',1.35,0.0015, 98.7751, 46.357625, ''),
        (6062,157.487654,-76.852859,'g',90.0,1,1,'2016/02/13 09:00:12.732',1.62,1.4089, 106.6357,49.368042, ''),
        (2206,345.810525,-70.197803,'g',90.0,2,1,'2016/06/30 08:37:52.759',1.31,1.6166, 91.3154, -8.633625, '')],
        dtype=[('HEX','<i8'),('RA','<f8'),('DEC','<f8'),('FILTER','S1'),
               ('EXPTIME','<f8'),('TILING','<i8'),('PRIORITY','<i8'),
               ('DATE','S30'), ('AIRMASS','<f8'),('SLEW','<f8'),
               ('MOONANGLE','<f8'),('HOURANGLE','<f8'),('PROGRAM','S30')])

    np.testing.assert_equal(fields[idx],test)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
