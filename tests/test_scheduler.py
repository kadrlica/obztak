#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os
import glob
import json
import subprocess
from subprocess import Popen, PIPE, check_call

def call(cmd,**kwargs):
    p = Popen(cmd,stdout=PIPE,stderr=PIPE,**kwargs)
    stdout,stderr = p.communicate()
    if p.returncode > 0:
        msg = str(cmd + '\n' + str(stdout) + '\n' + str(stderr))
        print(msg)
        raise subprocess.CalledProcessError(msg)
    print(stdout)
    return stdout

def check_dict(value,test):
    for k,v in value.items():
        val = test[k]
        if val != v:
            msg = '%s: %s (in), %s (out)'%(k,v,val)
            raise ValueError(msg)

    return test

def make_options(kwargs):
    return ' '.join(['--%s %s'%(k,v) for k,v in kwargs.items()])

def test_schedule_field():
    kwargs = dict(utc='2017-02-22T06:00:00',hex=6679,tiling=1,band='g',outfile='field_test.json')
    opts = make_options(kwargs)
    cmd = 'schedule_field %s'%(opts)
    call(cmd,shell=True)

    out = json.loads(open(kwargs['outfile']).read())
    assert len(out)==1

    index = 0
    value = {u'count': 1, u'program': u'maglites', u'seqtot': 2, u'seqnum': 1, u'expType': u'object', u'object': u'MAGLITES field: 6679-01-g', u'filter': u'g', u'seqid': u'MAGLITES scheduled: 2017/02/22 06:00:00.0000', u'RA': 121.083, u'propid': u'2016A-0366', u'dec': -57.0986, u'expTime': 90.0, u'wait': u'False'}
    check_dict(value,out[index])

def test_schedule_chunk():
    kwargs = dict(utc='2017-02-22T06:00:00',chunk=60,outfile='chunk_test.json',complete='None')
    opts = make_options(kwargs)
    cmd = 'schedule_chunk %s'%opts
    call(cmd,shell=True)

    out = json.loads(open(kwargs['outfile']).read())
    assert len(out)==30

    index = 17
    value = {u'count': 1, u'program': u'maglites', u'seqtot': 2, u'seqnum': 2, u'expType': u'object', u'object': u'MAGLITES field: 6679-01-r', u'filter': u'r', u'seqid': u'MAGLITES scheduled: 2017/02/22 06:34:00.0000', u'RA': 121.083, u'propid': u'2016A-0366', u'dec': -57.0986, u'expTime': 90.0, u'wait': u'False'}
    check_dict(value,out[index])

def test_schedule_night():
    kwargs = dict(nite='20170221',outfile='night_test/night_test.json',complete='None')
    opts = make_options(kwargs)
    cmd = 'schedule_night %s'%opts
    call(cmd,shell=True)

    files = sorted(glob.glob(os.path.dirname(kwargs['outfile'])+'/*.json'))
    assert len(files) == 5

    out = json.loads(open(files[2],'r').read())
    assert len(out)==30

    index = 17
    value = {u'count': 1, u'program': u'maglites', u'seqtot': 2, u'seqnum': 2, u'expType': u'object', u'object': u'MAGLITES field: 6891-01-r', u'filter': u'r', u'seqid': u'MAGLITES scheduled: 2017/02/22 07:40:00.0000', u'RA': 136.323, u'propid': u'2016A-0366', u'dec': -64.2085, u'expTime': 90.0, u'wait': u'False'}
    check_dict(value,out[index])

def NOTEST_schedule_survey():
    kwargs = dict(outfile='survey_test')
    opts = make_options(kwargs)
    cmd = 'schedule_survey %s'%opts
    call(cmd,shell=True)

    assert os.path.isdir(kwargs['outfile'])

    files = sorted(glob.glob(kwargs['outfile']+'/*/*.json'))
    assert len(files) == 131

    index = 17
    out = json.loads(open(files[index]).read())
    assert len(out) == 30

    index = 17
    value = {u'count': 1, u'program': u'maglites', u'seqtot': 2, u'seqnum': 2, u'expType': u'object', u'object': u'MAGLITES field: 6687-01-r', u'filter': u'r', u'seqid': u'MAGLITES scheduled: 2016/02/14 08:46:30.0000', u'RA': 176.263, u'propid': u'2016A-0366', u'dec': -74.1112, u'expTime': 90.0, u'wait': u'False'}
    check_dict(value,out[index])

if __name__ == "__main__":
    test_schedule_night()
