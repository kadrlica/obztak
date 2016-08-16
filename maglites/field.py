#!/usr/bin/env python
"""
Module for working with survey fields.
"""
import os
import copy
from collections import OrderedDict as odict
import logging

import numpy as np
from maglites import __version__
from maglites.utils import constants
from maglites.utils import fileio

DEFAULTS = odict([
    ('HEX',       dict(dtype=int,value=0)),
    ('RA',        dict(dtype=float,value=None)),
    ('DEC',       dict(dtype=float,value=None)),
    ('FILTER',    dict(dtype='S1',value='')),
    ('EXPTIME',   dict(dtype=float,value=90)),
    ('TILING',    dict(dtype=int,value=0)),
    ('PRIORITY',  dict(dtype=int,value=1)),
    ('DATE',      dict(dtype='S30',value='')),
    ('AIRMASS',   dict(dtype=float,value=-1.0)),
    ('SLEW',      dict(dtype=float,value=-1.0)),
    ('MOONANGLE', dict(dtype=float,value=-1.0)),
    ('HOURANGLE', dict(dtype=float,value=-1.0)),
])

DTYPES = odict([(k,v['dtype']) for k,v in DEFAULTS.items()])
VALUES = odict([(k,v['value']) for k,v in DEFAULTS.items()])

OBJECT_PREFIX = 'MAGLITES field: '
OBJECT_FMT = OBJECT_PREFIX + '%s'
SEQID_PREFIX = 'MAGLITES scheduled: '
SEQID_FMT = SEQID_PREFIX + '%(DATE)s'

SISPI_DICT = odict([
    ("object",  None),
    ("seqnum",  None), # 1-indexed
    ("seqtot",  2),
    ("seqid",   None),
    ("expTime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
    ("propid",  "2016A-0366"),
    ("comment", ""),
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

    def __getitem__(self,key):
        if key == 'ID':
            return self.unique_id
        else:
            return super(FieldArray,self).__getitem__(key)

    def append(self,other):
        return np.concatenate([self,other]).view(self.__class__)

    def keys(self):
        return self.dtype.names

    @property
    def unique_id(self):
        return np.char.mod('%(HEX)i-%(TILING)02d-%(FILTER)s',self)

    @property
    def object(self):
        return np.char.mod(OBJECT_FMT,self.unique_id).astype('S80')

    @property
    def seqid(self):
        return np.char.mod(SEQID_FMT,self).astype('S80')

    @property
    def seqnum(self):
        return np.array([constants.BANDS.index(f)+1 for f in self['FILTER']],dtype=int)

    @property
    def comment(self):
        comment = 'MAGLITES v%s: '%__version__
        comment += 'PRIORITY=%(PRIORITY)i, '

        fmt = '%s=%%(%s).4f'
        names = ['AIRMASS','SLEW','MOONANGLE','HOURANGLE']
        comment += ', '.join([fmt%(n,n) for n in names])
        return np.char.mod(comment,self)

    def from_unique_id(self,string):
        hex,tiling = map(int,string.split('-')[:2])
        self['HEX'] = hex
        self['TILING'] = tiling

    def from_object(self,string):
        self.from_unique_id(string.lstrip(OBJECT_PREFIX))

    def from_seqid(self, string):
        date = string.lstrip(SEQID_PREFIX)
        self['DATE'] = date

    def from_comment(self, string):
        integers = ['PRIORITY']
        floats   = ['AIRMASS','SLEW','MOONANGLE','HOURANGLE']
        values = dict([x.strip().split('=') for x in string.split(':')[-1].split(',')])
        for key,val in values.items():
            if key in integers:
                self[key] = int(val)
            elif key in floats:
                self[key] = float(val)
            elif key in strings:
                self[key] = str(val)
            else:
                msg = "Unrecognized comment field: %s"%key
                logging.warning(msg)

    def to_recarray(self):
        return self.view(np.recarray)

    def to_sispi(self):
        sispi = []
        objects = self.object
        seqnums = self.seqnum
        seqids = self.seqid
        comments = self.comment
        for i,r in enumerate(self):
            sispi_dict = copy.deepcopy(SISPI_DICT)
            for sispi_key,field_key in SISPI_MAP.items():
                sispi_dict[sispi_key] = r[field_key]
            sispi_dict['object'] = objects[i]
            sispi_dict['seqnum'] = seqnums[i]
            sispi_dict['seqid']  = seqids[i]
            sispi_dict['comment'] = comments[i]
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
            # Parse scheduled date if date is not present
            if 'date' in s: f['DATE'] = 'date'
            else: f.from_seqid(s['seqid'])
            f.from_comment(s['comment'])
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

        defaults = dict(propid='2016A-0366', limit='')
        params = copy.deepcopy(defaults)

        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec, 
        expTime, filter, 
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE, 
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE, 
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW 
        FROM exposure where propid = '%(propid)s' and exptime > 89 
        and discard = False and delivered = True and flavor = 'object'
        ORDER BY utc_beg %(limit)s 
        """%params

        data = database.execute(query)
        names = map(str.upper,database.get_columns())
        objidx = names.index('OBJECT')        
        if not len(data): return FieldArray(0)

        fields = cls()
        for d in data:
            f = cls(1)
            for i,key in enumerate(names):
                if key in f.dtype.names:
                    f[key] = d[i]
            f.from_object(d[objidx])
            fields = fields + f

        return fields

        
    @classmethod
    def read(cls, filename):
        base,ext = os.path.splitext(filename)
        if ext in ('.json'):
            sispi = fileio.read_json(filename)
            return cls().load_sispi(sispi)
        elif ext in ('.csv','.txt'):
            dtype = DTYPES.items()
            #recarray = fileio.csv2rec(filename,dtype=dtype)
            recarray = fileio.csv2rec(filename)
            return cls().load_recarray(recarray)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)

    def write(self, filename, **kwargs):
        base,ext = os.path.splitext(filename)
        logging.debug('Writing %s...'%filename)
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
