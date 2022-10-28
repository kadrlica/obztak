#!/usr/bin/env python
"""
Deal with file input/output.
"""
import os,pwd
from os.path import splitext, exists, join
from collections import OrderedDict as odict
import subprocess
import six

import numpy as np
import pandas as pd
import json
import logging

from obztak import __version__
from obztak.utils.constants import FLOAT_FMT
from obztak.utils.date import isstring

def get_username():
    import os,pwd
    return pwd.getpwuid( os.getuid() )[ 0 ]

def get_hostname():
    import platform
    return platform.node()

def get_datadir():
    from os.path import abspath,dirname,join
    return join(dirname(dirname(abspath(__file__))),'data')

def get_datafile(filename):
    dirname = get_datadir()
    filepath = os.path.join(dirname,filename)

    if not os.path.exists(filepath):
        msg = "File does not exists: %s"%filepath
        logging.warn(msg)
    else:
        return filepath

def read_csv(filename, **kwargs):
    """Wrapper around pandas.read_csv"""
    from distutils.version import LooseVersion
    kwargs.setdefault('parse_dates',False)
    kwargs.setdefault('comment','#')
    kwargs.setdefault('encoding',None)

    #if int(pd.__version__.replace('.','')) > 90:
    if LooseVersion(pd.__version__) > LooseVersion('0.9.0'):
        kwargs.setdefault('skip_blank_lines',True)
        #kwargs.setdefault('as_recarray',True)
        return pd.read_csv(filename,**kwargs)
    else:
        lines = open(filename,'r').readlines()
        comments = np.char.startswith(lines,'#')
        skiprows = np.argmin(comments)
        kwargs.setdefault('skiprows',skiprows)
        return pd.read_csv(filename,**kwargs)

def to_csv(filename,data,**kwargs):
    """ Call to pandas.DataFrame.to_csv
    """
    df = pd.DataFrame(data)

    kwargs.setdefault('float_format','%.4f')
    kwargs.setdefault('index',False)
    kwargs.setdefault('mode','w')
    kwargs.setdefault('na_rep','nan')

    if os.path.exists(filename):
        os.remove(filename)

    basename,ext = os.path.splitext(filename)
    if ext == '.gz': filename = basename

    with open(filename,'w') as out:
        out.write(header())
        df.to_csv(out,**kwargs)

    if ext == '.gz':
        cmd = 'gzip %s'%filename
        subprocess.check_call(cmd,shell=True)

def csv2rec(filename, **kwargs):
    """Read DataFrame from csv and return recarray
    """
    return read_csv(filename, **kwargs).to_records(index=False)

def rec2csv(filename,data,**kwargs):
    """Convert to DataFrame and write to csv.
    """
    return to_csv(filename, data, **kwargs)

def write_json(filename,data,**kwargs):
    kwargs.setdefault('indent',4)
    json.encoder.FLOAT_REPR = lambda o: format(o, '.4f')

    with open(filename,'w') as out:
        # It'd be nice to have a header
        #out.write(header())
        out.write(json.dumps(data,**kwargs))

def read_json(filename,**kwargs):
    with open(filename,'r') as f:
        return json.loads(f.read(),**kwargs)
            
def fields2sispi(filename,outfile=None,force=False):
    """ Convert a file of fields to a sispi json file.

    Parameters:
    -----------
    filename : input filename
    outfile  : output filename
    force    : overwrite output

    Returns:
    --------
    outfile  : output filename
    """
    if not outfile: outfile = splitext(filename)[0]+'.json'
    fields = FieldArray.read(filename)
    if exists(outfile) and not force:
        msg = "Output file already exists: %s"%(outfile)
        raise IOError(msg)
    logging.debug("Writing %s..."%outfile)
    fields.write(outfile)
    return outfile

def header():    
    import ephem
    from obztak.utils.date import datestring
    now = ephem.now()
    header  = "# author: %s@%s\n"%(get_username(),get_hostname())
    header += "# date: %s UTC\n"%(datestring(ephem.now(),0))
    header += "# version: obztak v%s\n"%(__version__)
    return header


def rec_append_fields(rec, names, arrs, dtypes=None):
    """
    Return a new record array with field names populated with data
    from arrays in *arrs*.  If appending a single field, then *names*,
    *arrs* and *dtypes* do not have to be lists. They can just be the
    values themselves.
    """
    if (not isstring(names) and iterable(names) and len(names) and isstring(names[0])):
        if len(names) != len(arrs):
            raise ValueError("number of arrays do not match number of names")
    else:  # we have only 1 name and 1 array
        names = [names]
        arrs = [arrs]
    arrs = list(map(np.asarray, arrs))
    if dtypes is None:
        dtypes = [a.dtype for a in arrs]
    elif not iterable(dtypes):
        dtypes = [dtypes]
    if len(arrs) != len(dtypes):
        if len(dtypes) == 1:
            dtypes = dtypes * len(arrs)
        else:
            raise ValueError("dtypes must be None, a single dtype or a list")
    old_dtypes = rec.dtype.descr
    if six.PY2:
        old_dtypes = [(name.encode('utf-8'), dt) for name, dt in old_dtypes]
    newdtype = np.dtype(old_dtypes + list(zip(names, dtypes)))
    newrec = np.recarray(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    for name, arr in zip(names, arrs):
        newrec[name] = arr
    return newrec

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
