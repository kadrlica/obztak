#!/usr/bin/env python
"""
Deal with file input/output
"""
from collections import OrderedDict as odict
from matplotlib import mlab
import numpy as np
import json

from maglites.utils.constants import FLOAT_FMT

class FormatFloatForce(mlab.FormatFormatStr): 
    """
    mlab not doing well...
    """
    def __init__(self,fmt=FLOAT_FMT): 
        mlab.FormatFormatStr.__init__(self,fmt) 
    def toval(self, x): 
        return x 
    def fromstr(self, s): 
        return float(s) 

def csv2rec(infile, **kwargs):
    #mlab.csv2rec(infile)
    data = np.recfromcsv(infile,**kwargs)
    data.dtype.names = map(str.upper,data.dtype.names)
    return data

def rec2csv(outfile,data,**kwargs):
    formatd = dict()
    for name,(dtype,size) in data.dtype.fields.items():
        if dtype.kind == 'f': formatd[name] = FormatFloatForce()

    formatd.update(kwargs.pop('formatd',dict()))

    mlab.rec2csv(data,outfile,formatd=formatd,**kwargs)

def write_json(outfile,data,**kwargs):
    kwargs.setdefault('indent',4)
    json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')
    with open(outfile,'w') as out:
        out.write(json.dumps(data,**kwargs))

def read_json(filename,**kwargs):
    with open(filename,'r') as f:
        return json.loads(f.read(),**kwargs)
            
    
    
if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
