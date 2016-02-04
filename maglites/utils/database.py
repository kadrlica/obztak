#!/usr/bin/env python
"""
Simple database interface taking connection
infromation from .desservices.ini
"""

import os
import psycopg2
import psycopg2.extras
import numpy as np
import logging

class Database(object):
    
    def __init__(self,dbname='db-fnal'):
        self.dbname = dbname
        self.conninfo = self.parse_config(section=dbname)
        self.connection = None
        self.cursor = None

    def __str__(self):
        ret = str(self.connection)
        return ret

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
        d['host'] = c.get(section,'server')
        d['dbname'] = c.get(section,'name')
        d['user'] = c.get(section,'user')
        d['password'] = c.get(section,'passwd')
        d['port'] = c.get(section,'port')
        return d

    def connect(self):
        #self.connection = psycopg2.connect(**CONNECTIONS[self.name])
        self.connection = psycopg2.connect(**self.conninfo)
        #self.cursor = self.connection.cursor(cursor_factory=psycopg2.extras.DictCursor)
        self.cursor = self.connection.cursor()

    def disconnect(self):
        self.cursor.close()
        self.cursor = None
        self.connection.close()
        self.connection = None

    def execute(self,query):
        self.cursor.execute(query)      
        return self.cursor.fetchall()

    def reset(self):
        self.connection.reset()
        
if __name__ == "__main__":
    import argparse
    description = "python script"
    parser = argparse.ArgumentParser(description=description)
    opts = parser.parse_args()

    db = Database()
    db.connect()
    print db
