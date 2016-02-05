#!/usr/bin/env python
"""
Deal with fileio...
"""
from collections import OrderedDict as odict
from matplotlib import mlab
import numpy as np

class FormatFloatForce(mlab.FormatFormatStr): 
    """
    mlab not doing well...
    """
    def __init__(self,fmt="%.4f"): 
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

def write_json(outfile,data,indent=4):
    # SISPI-compatible json output
    
    template = odict([
        
    ])



if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
