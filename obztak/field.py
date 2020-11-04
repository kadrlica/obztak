#!/usr/bin/env python
"""
Module for working with survey fields.
"""
import os
import copy
from collections import OrderedDict as odict
import logging

import numpy as np

from obztak import __version__
from obztak.utils import constants
from obztak.utils import fileio
from obztak.utils.date import setdefaults

# Default field array values
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
    ('PROGRAM',   dict(dtype='S30',value='')),
])
DTYPES = odict([(k,v['dtype']) for k,v in DEFAULTS.items()])
VALUES = odict([(k,v['value']) for k,v in DEFAULTS.items()])

# Separator for comments and sequences
SEP = ':'

# Default sispi dictionary
SISPI_DICT = odict([
    ("object",  None),
    ("seqnum",  None), # 1-indexed
    ("seqtot",  1),
    ("seqid",   ""),
    ("expTime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", None),
    ("wait",    "False"),
    ("propid",  None),
    ("comment", ""),
])

# Mapping between sispi dict keys and field array columns
SISPI_MAP = odict([
    ('expTime','EXPTIME'),
    ('RA','RA'),
    ('dec','DEC'),
    ('filter','FILTER'),
    ('program','PROGRAM'),
])

class FieldArray(np.recarray):
    """ Array for holding observation fields. """
    PROGRAM = 'maglites'
    PROPID  = '2016A-0366'

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT['program'] = PROGRAM
    SISPI_DICT['propid'] = PROPID

    OBJECT_FMT = 'MAGLITES field'+SEP+' %s'
    SEQID_FMT  = 'MAGLITES scheduled'+SEP+' %(DATE)s'

    BANDS = constants.BANDS

    def __new__(cls,shape=0):
        # Need to do it this way so that array can be resized...
        dtype = DTYPES.items()
        self = np.recarray(shape,dtype=dtype).view(cls)
        values = VALUES.items()
        for k,v in values: self[k].fill(v)
        return self
    
    #def __array_finalize__(self,obj):
    #    print('In array_finalize:')
    #    print('   self type is %s' % type(self))
    #    print('   obj type is %s' % type(obj))

    #def __array_wrap__(self, out_arr, context=None):
    #    print('In __array_wrap__:')
    #    print('   self is %s' % repr(self))
    #    print('   arr is %s' % repr(out_arr))

    def __add__(self, other):
        return np.concatenate([self,other]).view(self.__class__)

    def __getitem__(self,key):
        if isinstance(key,basestring) and key == 'ID':
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
    def field_id(self):
        return np.char.mod('%(HEX)i-%(TILING)02d',self)

    @property
    def object(self):
        return np.char.mod(self.OBJECT_FMT,self.unique_id).astype('S80')

    @property
    def seqid(self):
        return np.char.mod(self.SEQID_FMT,self).astype('S80')

    @property
    def seqnum(self):
        #return np.array([self.BANDS.index(f)+1 for f in self['FILTER']],dtype=int)
        return np.ones(len(self['FILTER']),dtype=int)

    @property
    def propid(self):
        return np.repeat(PROPID,len(self))

    @property
    def comment(self):
        comment = 'obztak v%s: '%__version__
        comment += 'PRIORITY=%(PRIORITY)i, '

        fmt = '%s=%%(%s).4f'
        names = ['AIRMASS','SLEW','MOONANGLE','HOURANGLE']
        comment += ', '.join([fmt%(n,n) for n in names])
        return np.char.mod(comment,self)

    def from_unique_id(self,string):
        try:
            hex,tiling = map(int,string.split('-')[:2])
            self['HEX'] = hex
            self['TILING'] = tiling
        except ValueError:
            logging.warn("Unparsed unique ID: '%s'"%string)
            self['HEX'] = -1
            self['TILING'] = -1

    def from_object(self,string):
        #self.from_unique_id(string.lstrip(OBJECT_PREFIX))
        self.from_unique_id(string.split(SEP,1)[-1].strip())

    def from_seqid(self, string):
        #date = string.lstrip(SEQID_PREFIX)
        if SEP not in string: return
        date = string.split(SEP,1)[-1].strip()
        self['DATE'] = date

    def from_comment(self, string):
        if SEP not in string: return
        integers = ['PRIORITY']
        floats   = ['AIRMASS','SLEW','MOONANGLE','HOURANGLE']
        values = dict([x.strip().split('=') for x in string.split(SEP,1)[-1].split(',')])
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
            sispi_dict = copy.deepcopy(self.SISPI_DICT)
            for sispi_key,field_key in SISPI_MAP.items():
                sispi_dict[sispi_key] = r[field_key]
            # Fill default program
            if not sispi_dict['program']:
                sispi_dict['program'] = self.SISPI_DICT['program']
            sispi_dict['object'] = objects[i]
            sispi_dict['seqnum'] = seqnums[i]
            sispi_dict['seqid']  = seqids[i]
            sispi_dict['comment'] = comments[i]
            sispi.append(sispi_dict)
        return sispi

    @classmethod
    def check_sispi(cls, sdict, check_propid=False):
        """Check that a dictionary is a valid sispi file

        Parameters
        ----------
        sdict        : sispi exposure dictionary
        check_propid : check that the propid matches this survey

        Returns
        -------
        check : boolean of whether check passed
        """
        # Ignore null exposures
        if sdict is None:
            logging.warn("Null exposure; skipping...")
            return False

        # Ignore exposures with the wrong propid
        propid = sdict.get('propid')
        if check_propid and propid != cls.PROPID: 
            logging.warn("Found exposure with propid=%s; skipping..."%propid)
            return False

        # Ignore exposures without RA,DEC columns
        if (sdict.get('RA') is None) or (sdict.get('dec') is None):
            logging.warn("RA,DEC not found; skipping...")
            return False

        return True
        
    @classmethod
    def load_sispi(cls,sispi,check_propid=False):
        fields = cls()
        # SISPI can do weird things...
        if (sispi is None) or (not len(sispi)): return fields
        for i,s in enumerate(sispi):
            # Check for some minimal exposure contents
            if not cls.check_sispi(s,check_propid): continue
            # SISPI can still do weird things...
            try:
                f = cls(1)
                for sispi_key,field_key in SISPI_MAP.items():
                    f[field_key] = s[sispi_key]
                f.from_object(s['object'])
                # Parse scheduled date if date is not present
                if 'date' in s: f['DATE'] = s['date']
                else: f.from_seqid(s['seqid'])
                f.from_comment(s['comment'])
                fields = fields + f
            #ADW: This is probably too inclusive...
            except (AttributeError,KeyError,ValueError,TypeError) as e: 
                logging.warn("Failed to load exposure\n%s"%s)
                logging.info(str(e))
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
    def query(cls, **kwargs):
        """ Generate the database query.

        Parameters:
        -----------
        kwargs : Keyword arguments to fill the query.

        Returns:
        --------
        query  : The query string.
        """
        defaults = dict(propid=cls.SISPI_DICT['propid'], limit='',
                        object_fmt = cls.OBJECT_FMT%'')
        kwargs = setdefaults(kwargs,copy.deepcopy(defaults))

        # Should pull this out to be accessible (self.query())?
        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW
        FROM exposure where propid = '%(propid)s' and exptime > 89
        and discard = False and delivered = True and flavor = 'object'
        and object like '%(object_fmt)s%%'
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query

    @classmethod
    def load_database(cls,database=None):
        """
        Load fields from the telemetry database.

        Parameters:
        -----------
        database : String or Database object to connect to.

        Returns:
        --------
        fields : A FieldArray filled from the database
        """
        try: from obztak.utils.database import Database
        except ImportError as e:
            logging.warn(e)
            return cls()

        try:
            database = Database(database)
        except IOError as e:
            logging.warn(e)
            return cls()

        database.connect()
        query = cls.query()

        logging.debug(query)
        data = database.execute(query)
        names = map(str.upper,database.get_columns())
        objidx = names.index('OBJECT')
        if not len(data):
            logging.warn("No fields found in database.")
            return cls()

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
    def load(cls, filename):
        base,ext = os.path.splitext(filename)
        if ext in ('.json'):
            sispi = fileio.read_json(filename)
            return cls().load_sispi(sispi)
        elif ext in ('.csv','.txt','.gz'):
            #dtype = DTYPES.items()
            #recarray = fileio.csv2rec(filename,dtype=dtype)
            recarray = fileio.csv2rec(filename)
            return cls().load_recarray(recarray)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)

    read = load

    def write(self, filename, **kwargs):
        base,ext = os.path.splitext(filename)
        logging.debug('Writing %s...'%filename)
        if ext in ('.json'):
            data = self.to_sispi()
            fileio.write_json(filename,data,**kwargs)
        elif ext in ('.csv','.txt','.gz'):
            data = self.to_recarray()
            fileio.rec2csv(filename,data,**kwargs)
        elif ext in ('.fits','.fz','.gz'):
            import fitsio
            data = self.to_recarray()
            kwargs.setdefault('clobber',True)
            fitsio.write(filename,data, **kwargs)
        else:
            msg = "Unrecognized file extension: %s"%ext
            raise IOError(msg)


class AllFieldArray(FieldArray):
    """ Array of all fields except DES and engineering. """
    PROGRAM  = 'all'
    PROPID   = 'none'
    PROPOSER = 'none'

    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID
    SISPI_DICT["proposer"] = PROPOSER

    OBJECT_FMT = '%s'
    SEQID_FMT = '%(DATE)s'

    BANDS = constants.BANDS

    @classmethod
    def query(cls, **kwargs):
        """ Generate the database query.

        Parameters:
        -----------
        kwargs : Keyword arguments to fill the query.

        Returns:
        --------
        query  : The query string.
        """
        from obztak.utils.date import setdefaults
        defaults = dict(propid=cls.SISPI_DICT['propid'], limit='')
        kwargs = setdefaults(kwargs,copy.deepcopy(defaults))

        query ="""
        SELECT object, seqid, seqnum, telra as RA, teldec as dec,
        expTime, filter,
        to_char(to_timestamp(utc_beg), 'YYYY/MM/DD HH24:MI:SS.MS') AS DATE,
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE,
        COALESCE(ha, -1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW
        FROM exposure where exptime > 59 and qc_teff > 0.1
        and propid != '2012B-0001' and propid not like '%%-9999'
        and discard = False and delivered = True and flavor = 'object'
        ORDER BY utc_beg %(limit)s
        """%kwargs
        return query


#def fields2sispi(infile,outfile=None,force=False):
#    if not outfile: outfile = os.path.splitext(infile)[0]+'.json'
#    fields = FieldArray.read(infile)
#    if os.path.exists(outfile) and not force:
#        msg = "Output file already exists: %s"%(outfile)
#        raise IOError(msg)
#    logging.debug("Writing %s..."%outfile)
#    fields.write(outfile)
#    return outfile

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
