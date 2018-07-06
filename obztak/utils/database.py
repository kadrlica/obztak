#!/usr/bin/env python
"""Simple interface interface to a postgres database taking
connection infromation from '.desservices.ini' file.

For more documentation on desservices, see here:
https://opensource.ncsa.illinois.edu/confluence/x/lwCsAw

"""

import os
import logging

import psycopg2
import numpy as np

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
            if hostname in ('observer2.ctio.noao.edu','observer3.ctio.noao.edu'):
                return 'db-ctio'
            else:
                return 'db-fnal'
        return dbname

    def parse_config(self, filename=None, section='db-fnal'):
        #if not filename: filename=os.getenv("DES_SERVICES")
        if os.path.exists(".desservices.ini"):
            filename=os.path.expandvars("$PWD/.desservices.ini")
        else:
            filename=os.path.expandvars("$HOME/.desservices.ini")
        logging.debug('.desservices.ini: %s'%filename)

        # ConfigParser throws "no section error" if file does not exist...
        # That's confusing, so 'open' to get a more understandable error
        open(filename) 
        import ConfigParser
        c = ConfigParser.RawConfigParser()
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
        except Exception, e:
            self.reset()
            raise(e)
        
    def reset(self):
        self.connection.reset()

    def get_columns(self,query=None):
        if query: self.select(query)
        return [d[0] for d in self.cursor.description]

    def query2recarray(self,query):
        # Doesn't work for all data types
        data = self.execute(query)
        names = self.get_columns()
        if not len(data):
            msg = "No data returned by query"
            raise ValueError(msg)
        return np.rec.array(data,names=names)

    query2rec = query2recarray

if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    db = Database()
    db.connect()
    print db
