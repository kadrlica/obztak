#!/usr/bin/env python
"""
Module for working with survey fields.
"""
import os
import copy
from collections import OrderedDict as odict

import numpy as np

from maglites.utils import constants
from maglites.utils import fileio

DEFAULTS = odict([
    ('ID',        dict(dtype=int,value=0)),
    ('OBJECT',    dict(dtype='S80',value='')),
    ('SEQID',     dict(dtype='S80',value='')),
    ('SEQNUM',    dict(dtype=int,value=0)),
    ('RA',        dict(dtype=float,value=None)),
    ('DEC',       dict(dtype=float,value=None)),
    ('FILTER',    dict(dtype='S1',value='')),
    ('EXPTIME',   dict(dtype=float,value=90)),
    ('TILING',    dict(dtype=int,value=0)),
    ('PRIORITY',  dict(dtype=int,value=1)),
    ('DATE',      dict(dtype='S20',value='')),
    ('AIRMASS',   dict(dtype=float,value=None)),
    ('SLEW',      dict(dtype=float,value=None)),
    ('MOONANGLE', dict(dtype=float,value=None)),
    ('HOURANGLE', dict(dtype=float,value=None)),
])

DTYPES = odict([(k,v['dtype']) for k,v in DEFAULTS.items()])
VALUES = odict([(k,v['value']) for k,v in DEFAULTS.items()])

SEPARATOR = ' : '
OBJECT_PREFIX = 'MAGLITES field'
OBJECT_FMT = OBJECT_PREFIX + SEPARATOR + '%(ID)i,%(TILING)1i,%(FILTER)1s'
SEQID_PREFIX = 'MAGLITES sequence'
SEQID_FMT = OBJECT_PREFIX + SEPARATOR + '%(DATE)s,%(PRIORITY)1i'

SISPI_DICT = odict([
    ("seqtot",  2),
    ("seqnum",  None), # 1-indexed
    ("seqid",   None),
    ("object",  None),
    ("exptime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
])

SISPI_MAP = odict([ 
    ('seqnum','SEQNUM'),
    ('seqid','SEQID'),
    ('object','OBJECT'),
    ('exptime','EXPTIME'),
    ('RA','RA'),
    ('filter','FILTER'),
])

class Field(odict):
    def __init__(self, **kwargs):
        super(Field,self).__init__(**kwargs)
        self.update(VALUES)

    @classmethod
    def from_recarray(cls,recarray):
        field = cls()
        field.load_recarray(recarray)
        return field

    @classmethod
    def from_sispi(cls,sispi):
        field = cls()
        field.load_sispi(sispi)
        return field

    def finalize(self, names):
        if 'OBJECT' not in names: self.fill_object()
        elif 'ID' not in names: self.from_object()

        if 'SEQID' not in names: self.fill_seqid()
        elif 'DATE' not in names: self.from_seqid()

        if 'COUNT' not in names: self.fill_count()

    def load_recarray(self,recarray): 
        if len(recarray) != 1:
            msg = "Incorrect recarray length: %s"%len(recarray)
            raise Exception(msg)

        keys = dict([(n.upper(),n) for n in recarray.dtype.names])
        for k in self.keys():
            self[k] = recarray[keys[k]]

        self.finalize(keys.keys())
            
    def load_dict(self,d): 
        keys  = dict([(k.upper(),k) for k in d.keys()])
        names = map(str.upper,dict.keys())
        for k in self.keys():
            self[k] = d[keys[k]]

        self.finalize(keys.keys())

    def load_sispi(self,sispi): 
        for sispi_key,field_key in SISPI_MAP.items():
            self[field_key] = sispi[sispi_key]

        self.finalize(SISPI_MAP.values())

        
    def to_recarray(self): 
        recarray = np.recarray(1, dtype=DTYPES)
        for key in recarray.dtypes.names:
            recarray[key] = self[key]
        return recarray

    def to_dict(self): 
        return dict(copy.deepcopy(self))

    def to_sispi(self): 
        sispi = copy.deepcopy(SISPI_DICT)
        for sispi_key,field_key in SISPI_MAP.items():
            sispi[sispi_key] = self[field_key]
        return sispi

    def fill_object(self):
        self['OBJECT'] =  OBJECT_FMT%self

    def from_object(self,string):
        id,tiling,filter = string.split(SEPARATOR)[-1].split(',')
        self['ID'] = int(id)
        self['TILING'] = int(tiling)
        self['FILTER'] = str(filter)
        
    def fill_seqid(self):
        self['SEQID'] = SEQID_FMT%self

    def from_seqid(self, string):
        date,priority = string.split(SEPARATOR)[-1].split(',')
        self['DATE'] = str(date)
        self['PRIORITY'] = int(priority)

    def fill_seqnum(self):
        self['SEQNUM'] = constants.BANDS.index(self['FILTER'])

class FieldArray(np.recarray):

    def __new__(cls):
        # Need to do it this way so that array can be resized...
        dtype = DTYPES.items()
        self = np.recarray(0,dtype=dtype).view(cls)
        return self
    
    def __add__(self, other):
        return np.concatenate([self,other])

    def to_sispi(self):
        return [f.to_sispi() for f in self]

    def to_recarrayi(self):
        return self.view(np.recarray)

    @classmethod
    def from_sispi(cls, sispi):
        array = cls()
        for s in sispi:
            f = Field.from_sispi(s)
            array = np.append(array,f)
        return array

    @classmethod
    def from_recarray(cls, recarray):
        array = cls()
        for r in recarray:
            f = Field.from_recarray(r)
            array = np.append(array,f)
        return array

    @classmethod
    def from_file(cls, filename):
        base,ext = os.path.splitext(filename)
        if ext in ('.json'):
            sispi = fileio.read_json(filename)
            return cls.from_sispi(sispi)
        elif ext in ('.csv','.txt'):
            dtype = DTYPES.items()
            recarray = np.genfromtxt(filename,delimiter=',',names=True,dtype=dtype)
            return cls.from_recarray(recarray)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise Exception(msg)

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
