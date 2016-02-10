#!/usr/bin/env python
"""
Module for working with survey fields.
"""
import os
import copy
from collections import OrderedDict as odict
import logging

import numpy as np

from maglites.utils import constants
from maglites.utils import fileio

DEFAULTS = odict([
    ('ID',        dict(dtype=int,value=0)),
    ('SMASH_ID',  dict(dtype=int,value=0)),
    #('OBJECT',    dict(dtype='S80',value='')),
    #('SEQID',     dict(dtype='S80',value='')),
    #('SEQNUM',    dict(dtype=int,value=0)),
    ('RA',        dict(dtype=float,value=None)),
    ('DEC',       dict(dtype=float,value=None)),
    ('FILTER',    dict(dtype='S1',value='')),
    ('EXPTIME',   dict(dtype=float,value=90)),
    ('TILING',    dict(dtype=int,value=0)),
    ('PRIORITY',  dict(dtype=int,value=1)),
    ('DATE',      dict(dtype='S20',value='')),
    ('AIRMASS',   dict(dtype=float,value=-1.0)),
    ('SLEW',      dict(dtype=float,value=-1.0)),
    ('MOONANGLE', dict(dtype=float,value=-1.0)),
    ('HOURANGLE', dict(dtype=float,value=-1.0)),
])

DTYPES = odict([(k,v['dtype']) for k,v in DEFAULTS.items()])
VALUES = odict([(k,v['value']) for k,v in DEFAULTS.items()])

SEPARATOR = ' : '
OBJECT_PREFIX = 'MAGLITES field'
OBJECT_FMT = OBJECT_PREFIX + SEPARATOR + '%s'
SEQID_PREFIX = 'MAGLITES sequence'
SEQID_FMT = SEQID_PREFIX + SEPARATOR + '%(DATE)s'

SISPI_DICT = odict([
    ("seqtot",  2),
    ("seqnum",  None), # 1-indexed
    ("seqid",   None),
    ("object",  None),
    ("expTime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
    ("propid",  "2016A-0366"),
])

SISPI_MAP = odict([ 
    ('expTime','EXPTIME'),
    ('RA','RA'),
    ('dec','DEC'),
    ('filter','FILTER'),
])

class FieldArray(np.recarray):

    def __new__(cls,shape=0):
        # Need to do it this way so that array can be resized...
        dtype = DTYPES.items()
        self = np.recarray(shape,dtype=dtype).view(cls)
        values = VALUES.items()
        for k,v in values: self[k].fill(v)
        return self
    
    def __add__(self, other):
        return np.concatenate([self,other]).view(self.__class__)

    def append(self,other):
        return np.concatenate([self,other]).view(self.__class__)

    def keys(self):
        return self.dtype.names

    @property
    def index(self):
        return np.char.mod('%(SMASH_ID)i.%(TILING)02d',self)

    @property
    def object(self):
        return np.char.mod(OBJECT_FMT,self.index).astype('S80')

    @property
    def seqid(self):
        return np.char.mod(SEQID_FMT,self).astype('S80')

    @property
    def seqnum(self):
        return np.array([constants.BANDS.index(f)+1 for f in self['FILTER']],dtype=int)

    def from_index(self,string):
        smash_id,tiling = map(int,string.split('.'))
        self['SMASH_ID'] = smash_id
        self['TILING'] = tiling

    def from_object(self,string):
        self.from_index(string.split(SEPARATOR)[-1])

    def from_seqid(self, string):
        date = string.split(SEPARATOR)[-1].split(',')
        self['DATE'] = date

    def to_recarray(self):
        return self.view(np.recarray)

    def to_sispi(self):
        sispi = []
        objects = self.object
        seqnums = self.seqnum
        seqids = self.seqid
        for i,r in enumerate(self):
            sispi_dict = copy.deepcopy(SISPI_DICT)
            for sispi_key,field_key in SISPI_MAP.items():
                sispi_dict[sispi_key] = r[field_key]
            sispi_dict['object'] = objects[i]
            sispi_dict['seqnum'] = seqnums[i]
            sispi_dict['seqid']  = seqids[i]
            sispi.append(sispi_dict)
        return sispi

    @classmethod
    def load_sispi(cls,sispi):
        fields = cls()
        for i,s in enumerate(sispi):
            f = cls(1)
            for sispi_key,field_key in SISPI_MAP.items():
                f[field_key] = s[sispi_key]
            f.from_object(s['object'])
            # Only do this if DATE is not present
            f.from_seqid(s['seqid'])
            fields = fields + f
        return fields

    @classmethod
    def load_recarray(cls,recarray): 
        fields = cls(len(recarray))
        keys = dict([(n.upper(),n) for n in recarray.dtype.names])

        for k in fields.dtype.names:
            if k not in keys: 
                logging.warning('Key %s not found in input array'%k)
                continue
            fields[k] = recarray[keys[k]]
        return fields

    @classmethod
    def load_database(cls,database='db-fnal'):
        """
        Get the fields that have been observed from the telemetry DB.
        """
        from maglites.utils import Database

        if not isinstance(database,Database):
            database = Database(database)
            database.connect()

        defaults = dict(propid='2016A-0366', limit='', dbname='db-fnal')
        params = copy.deepcopy(defaults)
        params.update(kwargs)
            
        db = Database()
        db.connect()

        query ="""
        select object, seqid, seqnum, telra as RA, teldec as dec, 
        expTime, filter, to_char(utc_beg, 'YYYY/MM/DD HH24:MI:SS') AS DATE, 
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE, 
        COALESCE(ha,-1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW,
        from exposure where propid = '%(propid)s' %(limit)s
        """%params

        data = database.execute(query)
        names = map(str.upper,database.get_columns())
        
        if not len(data): return FieldArray(0)
        sispi = np.recarray(data,names=names)

        fields = self.load_sispi(sispi)
        return fields
        
    @classmethod
    def read(cls, filename):
        base,ext = os.path.splitext(filename)
        if ext in ('.json'):
            sispi = fileio.read_json(filename)
            return cls().load_sispi(sispi)
        elif ext in ('.csv','.txt'):
            dtype = DTYPES.items()
            recarray = fileio.csv2rec(filename,dtype=dtype)
            return cls().load_recarray(recarray)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)

    def write(self, filename, **kwargs):
        base,ext = os.path.splitext(filename)
        if ext in ('.json'):
            data = self.to_sispi()
            fileio.write_json(filename,data,**kwargs)
        elif ext in ('.csv','.txt'):
            data = self.to_recarray()
            fileio.rec2csv(filename,data,**kwargs)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)


def fields2sispi(infile,outfile=None,force=False):
    if not outfile: outfile = os.path.splitext(infile)[0]+'.json'
    fields = FieldArray.read(infile)
    if os.path.exists(outfile) and not force:
        msg = "Output file already exists: %s"%(outfile)
        raise IOError(msg)
    logging.debug("Writing %s..."%outfile)
    fields.write(outfile)
    return outfile

            
if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
