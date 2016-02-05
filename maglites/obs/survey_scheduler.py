#!/usr/bin/env python
"""
Module for running survey operations.
"""
import os,sys
import logging
import copy
import numpy as np
import ephem

from maglites.sim import Simulator
from maglites.utils import Database, datestring
import maglites.utils.constants as constants

class Scheduler(Simulator):

    def __init__(self, infile_target_fields):
        super(Scheduler,self).__init__(infile_target_fields)
        observed_fields = self.getObservedFields()
        if observed_fields:
            self.accomplished_fields = np.append(self.accomplished_fields,observed_fields)
        self.scheduled_fields = self.createFieldArray()
        
    def run(self, tstart=None, tstop=None, plot=True):
        # If no tstop, run for 90 minutes
        timedelta = 90*ephem.minute
        if tstart is None: tstart = ephem.now()
        if tstop is None: tstop = tstart + timedelta

        # Convert strings into dates
        if isinstance(tstart,basestring):
            tstart = ephem.Date(tstart)
        if isinstance(tstop,basestring):
            tstop = ephem.Date(tstop)

        msg = "Previously accomplished fields: %i"%len(self.accomplished_fields)
        logging.info(msg)

        date = tstart
        latch = True
        while latch:
            # Check to see if in valid observation window
            if self.observation_windows is not None:
                inside = False
                for window in self.observation_windows:
                    if date >= window[0] and date <= window[-1]: 
                        inside = True 

                if not inside:
                    msg = 'Date is outside  nominal observing window: %s'%datestring(date)
                    logging.warning(msg)


            if plot:
                print '  %s'%(datestring(date)),
            else:
                print '  %s'%(datestring(date))
            
            # Check 
            compute_slew = True
            if len(self.accomplished_fields['ID']) == 0:
                compute_slew = False
            else:
                if (date - ephem.Date(self.accomplished_fields['DATE'][-1])) > (30. * ephem.minute):
                    compute_slew = False
            if compute_slew:
                field_select = self.selectField(date, ra_previous=self.accomplished_fields['RA'][-1], dec_previous=self.accomplished_fields['DEC'][-1], plot=plot)
            else:
                field_select = self.selectField(date, plot=plot)
                
            id_select = field_select['ID']

            date = date + constants.FIELDTIME

            self.accomplished_field_ids.append(id_select)
            new_field = self.createFieldArray(1)
            for key in field_select.keys():
                new_field[-1][key] = field_select[key]

            self.accomplished_fields = np.append(self.accomplished_fields,new_field)
            self.scheduled_fields = np.append(self.scheduled_fields,new_field)

            if date > tstop: break
        print "Newly scheduled fields: %i"%len(self.scheduled_fields)
        return self.scheduled_fields

        self.accomplished_field_ids = []

    def getObservedFields(self, **kwargs):
        """
        Get the fields that have been observed from the telemetry DB.
        """
        defaults = dict(propid='2016A-0366', limit='', dbname='db-fnal')
        
        params = dict(defaults)
        #if opts is not None: params.update(vars(opts))
        
        db = Database()
        db.connect()

        #From survey_simulator, the expected columns are:
        # ['ID', 'RA', 'DEC', 'TILING', 'PRIORITY', 'DATE', 'AIRMASS', 'SLEW', 'MOONANGLE', 'HOURANGLE']
        query ="""
        select object as ID, telra as RA, teldec as DEC, 
        1 as TILING, 1 as PRIORITY, 
        to_char(utc_beg, 'YYYY/MM/DD HH24:MI:SS') AS DATE, 
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE, 
        COALESCE(ha,-1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW 
        from exposure where propid = '%(propid)s' %(limit)s
        """%params

        data = db.execute(query)
        names = map(str.upper,db.get_columns())

        if len(data): ret = np.rec.array(data,names=names)
        else:         ret = None
        return ret

    @classmethod
    def parser(cls):
        parser = cls.common_parser()
        parser.add_argument('--tstart',
                            help="Start time for observation.")
        parser.add_argument('--tstop' ,
                            help="Stop time for observation.")

        return parser
    


def main():
    args = Scheduler.parser().parse_args()

    obs = Scheduler(args.fields)
    obs.loadObservationWindows(args.windows)
    obs.loadAccomplishedFields(args.done)
    obs.run(args.tstart,args.tstop,plot=args.plot)
    if args.outfile: 
        obs.saveFields(args.outfile,obs.scheduled_fields)

    if not sys.flags.interactive and args.plot:
        raw_input(' ...finish...')

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
    
    main()
    
