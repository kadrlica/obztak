#!/usr/bin/env python
"""
Deal with file input/output
"""
import os,pwd
from os.path import splitext, exists, join
from collections import OrderedDict as odict
from matplotlib import mlab
import numpy as np
import json
import logging

from maglites import __version__
from maglites.utils.constants import FLOAT_FMT

#from maglites.field import FieldArray

def get_username():
    import os,pwd
    return pwd.getpwuid( os.getuid() )[ 0 ]

def get_hostname():
    import platform
    return platform.node()

def get_datadir():
    from os.path import abspath,dirname,join
    return join(dirname(dirname(abspath(__file__))),'data')

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

def csv2rec(filename, **kwargs):
    #mlab.csv2rec(infile)
    #data = np.recfromcsv(filename,**kwargs)
    #data.dtype.names = map(str.upper,data.dtype.names)
    
    import pandas as pd
    from distutils.version import LooseVersion
    kwargs.setdefault('parse_dates',False)
    kwargs.setdefault('comment','#')

    #if int(pd.__version__.replace('.','')) > 90:
    if LooseVersion(pd.__version__) > LooseVersion('0.9.0'):
        kwargs.setdefault('skip_blank_lines',True)
        #kwargs.setdefault('as_recarray',True)
        return pd.read_csv(filename,**kwargs).to_records(index=False)
    else:
        lines = open(filename,'r').readlines()
        comments = np.char.startswith(lines,'#')
        skiprows = np.argmin(comments)
        kwargs.setdefault('skiprows',skiprows)
        data = pd.read_csv(filename,**kwargs).to_records(index=False)
        return data
        
def rec2csv(filename,data,**kwargs):
    """
    Wrapper around numpy.savetxt

    Also see mlab.rec2csv (which is terrible...)
    """
    #formatd = dict()
    #for name,(dtype,size) in data.dtype.fields.items():
    #    if dtype.kind == 'f': formatd[name] = FormatFloatForce()
    #formatd.update(kwargs.pop('formatd',dict()))
    #
    #mlab.rec2csv(data,out,formatd=formatd,**kwargs)        

    import pandas as pd
    df = pd.DataFrame(data)

    kwargs.setdefault('float_format','%.4f')
    kwargs.setdefault('index',False)
    kwargs.setdefault('mode','w')
    kwargs.setdefault('na_rep','nan')
    
    with open(filename,'wb') as out:
        out.write(header())
        df.to_csv(out,**kwargs)

    #mlab.rec2csv(data,outfile,formatd=formatd,**kwargs)
    

def write_json(outfile,data,**kwargs):
    kwargs.setdefault('indent',4)
    json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    with open(outfile,'wb') as out:
        # It'd be nice to have a header
        #out.write(header())
        out.write(json.dumps(data,**kwargs))

def read_json(filename,**kwargs):
    with open(filename,'r') as f:
        return json.loads(f.read(),**kwargs)
            
def fields2sispi(infile,outfile=None,force=False):
    if not outfile: outfile = splitext(infile)[0]+'.json'
    fields = FieldArray.read(filename)
    if exists(outfile) and not force:
        msg = "Output file already exists: %s"%(outfile)
        raise IOError(msg)
    logging.debug("Writing %s..."%outfile)
    fields.write(outfile)
    return outfile

def header():    
    import ephem
    now = ephem.now()
    header  = "# author: %s@%s\n"%(get_username(),get_hostname())
    header += "# date: %s UTC\n"%(ephem.now())
    header += "# version: maglites v%s\n"%(__version__)
    return header
    
if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
