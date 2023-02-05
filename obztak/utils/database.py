#!/usr/bin/env python
"""Simple interface interface to a postgres database taking
connection infromation from '.desservices.ini' file.

For more documentation on desservices, see here:
https://opensource.ncsa.illinois.edu/confluence/x/lwCsAw

"""

import os
import logging
import datetime

import psycopg2
import pandas as pd
import numpy as np

try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

def desc2dtype(desc):
    """ Covert from postgres type description to numpy dtype.
    Tries to conform to the type mapping on the psycopg2 documentation:
    https://pythonhosted.org/psycopg2/usage.html"""
    # A list of all string types can be found in:
    # from psycopg2.extensions import string_types
    dtype = []
    for d in desc:
        name = d.name
        code = d.type_code
        size = d.internal_size
        # Booleans
        if code == psycopg2.extensions.BOOLEAN:
            dt = (name, bool)
        # Numeric types
        elif code == psycopg2.extensions.LONGINTEGER:
            dt = (name, long)
        elif code == psycopg2.extensions.INTEGER:
            dt = (name, 'i%i'%size)
        elif code == psycopg2.extensions.FLOAT:
            dt = (name, 'f%i'%size)
        elif code == psycopg2.NUMBER:
            # Catchall for other numbers
            dt = (name, float)
        # Character strings
        elif code == psycopg2.STRING:
            if size > 0:
                dt = (name, 'S%i'%size)
            else:
                # These are TEXT objects of undefined length
                dt = (name, object)
        elif code == psycopg2.extensions.UNICODE:
            # Probably doesn't get called because STRING is a catchall
            dt = (name, 'U%i'%size)
        # Dates and times (should eventually move to np.datetime64)
        elif code == psycopg2.extensions.DATE:
            dt = (name, datetime.date)
        elif code == psycopg2.extensions.TIME:
            dt = (name, datetime.time)
        elif code == psycopg2.extensions.INTERVAL:
            dt = (name, datetime.timedelta)
        elif code in (psycopg2.DATETIME, 1184):
            dt = (name, datetime.datetime)
        elif code == psycopg2._psycopg.UNKNOWN:
            dt = (name, object)
        # Binary stuff
        elif code == psycopg2.BINARY:
            dt = (name, bytearray)
        elif code == psycopg2.ROWID:
            dt = (name, bytearray)
        else: # Ignore other types for now.
            msg = "Unrecognized type code: "+str(d)
            raise TypeError(msg)
        dtype.append(dt)
    return dtype


class Database(object):
    
    def __init__(self,dbname=None):
        self.dbname = self.parse_dbname(dbname)
        self.conninfo = self.parse_config(section=self.dbname)
        self.connection = None
        self.cursor = None

    def __str__(self):
        ret = str(self.connection)
        return ret

    def parse_dbname(self, dbname):
        import platform
        if dbname is None:
            hostname = platform.node()
            # Only from ctio machines
            if hostname in ('observer2.ctio.noao.edu','observer3.ctio.noao.edu'):
            # Any machine on ctio network
            #if hostname.endswith('.ctio.noao.edu'):
                return 'db-ctio'
            else:
                return 'db-fnal'
        return dbname

    def parse_config(self, filename=None, section='db-fnal'):
        if filename: pass
        elif os.getenv("DES_SERVICES"):
            filename=os.getenv("DES_SERVICES")
        elif os.path.exists(".desservices.ini"):
            filename=os.path.expandvars("$PWD/.desservices.ini")
        else:
            filename=os.path.expandvars("$HOME/.desservices.ini")
        logging.debug('.desservices.ini: %s'%filename)
        if not os.path.exists(filename): raise IOError("%s does not exist"%filename)

        # ConfigParser throws "no section error" if file does not exist...
        # That's confusing, so 'open' to get a more understandable error
        #open(filename)
        c = RawConfigParser()
        c.read(filename)

        d={}
        d['host']     = c.get(section,'server')
        d['dbname']   = c.get(section,'name')
        d['user']     = c.get(section,'user')
        d['password'] = c.get(section,'passwd')
        d['port']     = c.get(section,'port')
        return d

    def connect(self):
        logging.debug("Connecting to: %s"%self.dbname)
        self.connection = psycopg2.connect(**self.conninfo)
        self.cursor = self.connection.cursor()
        logging.debug("Connection: %s"%str(self.connection).split("'")[1])

    def disconnect(self):
        self.cursor.close()
        self.cursor = None
        self.connection.close()
        self.connection = None

    def execute(self,query):
        self.cursor.execute(query)      
        try: 
            return self.cursor.fetchall()
        except Exception as e:
            self.reset()
            raise(e)
        
    def reset(self):
        self.connection.reset()

    def get_description(self,query=None):
        if query: self.select(query)
        return self.cursor.description

    def get_columns(self,query=None):
        return [d[0] for d in self.get_description(query)]

    def get_dtypes(self,query=None):
        desc = self.get_description(query)
        return desc2dtype(desc)

    def query2recarray(self,query):
        # Doesn't work for all data types
        data = self.execute(query)
        names = self.get_columns()
        dtypes = self.get_dtypes()
        if not len(data):
            msg = "No data returned by query"
            #raise ValueError(msg)
            return np.recarray(0,dtype=dtypes)

        #return np.rec.array(data,names=names)
        return np.rec.array(data,dtype=dtypes)

    query2rec = query2recarray


    def qcInv(self, timedelta=None, propid=None):
        """Get qc information.
     
        Parameters:
        -----------
        timedelta : time interval to query.
        propid    : proposal id to select on.
     
        Returns:
        --------
        df        : pd.DataFrame with query results
        """
        if timedelta is None: timedelta = '12h'
        propid = '' if propid is None else "and propid = '%s'"%propid
        
        query = """
        SELECT 
        id as expnum, telra as ra, teldec as dec, 
        to_char(date, 'HH24:MI') AS utc,
        filter as fil, CAST(exptime AS INT) as time, airmass as secz,
        qc_fwhm as psf, qc_sky as sky, qc_cloud as cloud, qc_teff as teff,
        object
        FROM exposure
        WHERE
        flavor = 'object' and
        date > (now() - interval '{timedelta}') {propid}
        ORDER BY expnum ASC
        """.format(timedelta=timedelta, propid=propid)
         
        return pd.DataFrame(self.query2rec(query))

    def duplicates(self, timedelta=None, propid=None):
        """Get qc information for duplicate exposures.

        Parameters:
        -----------
        timedelta : time interval to query.
        propid    : proposal id to select on.

        Returns:
        --------
        df        : pd.DataFrame with query results
        """
        if timedelta is None: timedelta = '12h'
        propid = '' if propid is None else "and propid = '%s'"%propid

        query = """
        SELECT id as expnum, telra as ra, teldec as dec,
        to_char(date, 'HH24:MI') AS utc,
        filter as fil, CAST(exptime AS INT) as time, airmass as secz,
        qc_fwhm as psf, qc_sky as sky, qc_cloud as cloud, qc_teff as teff,
        e.object
        FROM exposure e,
          (SELECT object from exposure where
           ( date > (now() - interval '{timedelta}') )
           GROUP BY object having count(object)>1) AS t
        WHERE t.object = e.object and e.flavor = 'object'
        {propid}
        AND ( e.date > (now() - interval '{timedelta}') )
        ORDER BY object, id
        """.format(timedelta=timedelta, propid=propid)

        return pd.DataFrame(self.query2rec(query))


if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    db = Database()
    db.connect()
    print(db)
