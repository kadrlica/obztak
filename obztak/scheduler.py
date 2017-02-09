"""
Module providing the survey scheduler.
"""

import os,sys
import copy
import numpy as np
import time
import ephem
import matplotlib.pyplot as plt
import logging
from collections import OrderedDict as odict

import obztak.utils.projector
import obztak.utils.constants
import obztak.utils.ortho

from obztak.utils import constants
from obztak.utils import ortho
from obztak.utils import fileio

from obztak.ctio import CTIO
from obztak.field import FieldArray
from obztak.tactician import CoverageTactician
from obztak.utils.date import get_nite, datestr, datestring, nitestring, utc2nite
from obztak.factory import tactician_factory

# For debugging (can also use the verbose command line argument)
#logging.basicConfig(level=20) # KCB

############################################################

class Scheduler(object):
    """
    Deal with survey scheduling.
    """
    _defaults = odict([
        ('tactician','coverage'),
        ('windows',os.path.join(fileio.get_datadir(),"maglites-windows.csv")),
        ('targets',os.path.join(fileio.get_datadir(),"maglites-target-fields.csv")),
    ])
    FieldType = FieldArray

    def __init__(self,target_fields=None,windows=None,completed_fields=None):
        self.load_target_fields(target_fields)
        self.load_windows(windows)
        self.load_observed_fields()
        self.load_completed_fields(completed_fields)

        self.scheduled_fields = self.FieldType()
        self.observatory = CTIO()

    def load_target_fields(self, target_fields=None):
        if target_fields is None:
            target_fields = self._defaults['targets']

        if isinstance(target_fields,basestring):
            self.target_fields = self.FieldType.read(target_fields)
        else:
            self.target_fields = self.FieldType(target_fields)
        return self.target_fields

    def load_windows(self, windows=None):
        """
        Load the set of start and stop times for the observation windows.
        """
        if windows is None:
            windows = self._defaults['windows']
            logging.info("Setting default observing windows:\n %s"%windows)

        if isinstance(windows,basestring):
            windows = fileio.csv2rec(windows)

        self.windows = []
        for start,end in windows:
            self.windows.append([ephem.Date(start), ephem.Date(end)])

        # Sanity check that observation windows are properly sorted
        for ii,(start,end) in enumerate(self.windows):
            msg = 'Observation windows are not properly sorted\n'
            msg+= '%s -- %s'%(datestr(start),datestr(end))
            if (end < start):
                logging.warn(msg)
            if ii > 0 and (start < self.windows[ii-1][1]):
                logging.warn(msg)

        logging.info('Observation Windows:')
        for start,end in self.windows:
            logging.info(' %s UTC -- %s UTC'%(datestr(start),datestr(end)))
        logging.info(30*'-')

    def load_observed_fields(self):
        """
        Load fields from the telemetry database that were already observed.
        """
        try:
            fields = self.FieldType.load_database()
        except Exception as e:
            logging.warn("Failed to load completed exposures from database")
            logging.info(e)
            fields = self.FieldType()
        self.observed_fields = fields
        return self.observed_fields

    def load_completed_fields(self, completed_fields=None):
        """Load completed fields. The default behavior is to load the
        observed_fields as completed_fields. However, if the string
        'None' is passed then return an empty FieldArray.

        Parameters:
        -----------
        completed_fields : Filename, list of filenames, or FieldArray-type object.

        Returns:
        --------
        fields           : FieldArray of the completed fields
        """
        # Deal with 'None' string
        if isinstance(completed_fields,list):
            if completed_fields[0].lower()=='none':
                self.completed_fields = self.FieldType()
                return self.completed_fields
        elif isinstance(completed_fields,basestring):
            if completed_fields.lower()=='none':
                self.completed_fields = self.FieldType()
                return self.completed_fields

        self.completed_fields = copy.deepcopy(self.observed_fields)

        if not completed_fields:
            return self.completed_fields

        if isinstance(completed_fields,basestring):
            completed_fields = [completed_fields]

        if isinstance(completed_fields,list):
            fields = self.FieldType()
            for filename in completed_fields:
                fields = fields + self.FieldType.read(filename)

            completed_fields = fields

        new=~np.in1d(completed_fields.unique_id,self.completed_fields.unique_id)
        new_fields = completed_fields[new]
        self.completed_fields = self.completed_fields + new_fields
        return self.completed_fields

    def create_tactician(self,tactician=None):
        if tactician is None: tactician = self._defaults['tactician']
        return tactician_factory(tactician,mode=tactician)

    def select_field(self, date, mode='coverage'):
        """
        Select field(s) using the survey tactician.

        Parameters:
        -----------
        date       : ephem.Date object
        mode       : Type of tactician to use for selecting field

        Returns:
        --------
        field      : selected field(s) from tactician
        """
        sel = ~np.in1d(self.target_fields['ID'],self.completed_fields['ID'])

        self.tactician = self.create_tactician(mode)
        self.tactician.set_date(date)
        self.tactician.set_target_fields(self.target_fields[sel])
        self.tactician.set_completed_fields(self.completed_fields)

        field_select = self.tactician.select_fields()

        logging.debug(str(field_select))

        # For diagnostic purposes
        if False and len(self.scheduled_fields) % 10 == 0:
            weight = self.tactician.weight
            ortho.plotWeight(field_select[-1], self.target_fields, self.tactician.weight)
            raw_input('WAIT')

        if len(field_select) == 0:
            logging.error("No field selected... we've got problems.")
            msg  = "date=%s\n"%(datestr(date))
            msg += "index_select=%s, index=%s\n"%(index_select,index)
            msg += "nselected=%s, selection=%s\n"%(cut.sum(),cut[index_select])
            msg += "weights=%s"%weight
            logging.info(msg)
            #ortho.plotWeight(self.scheduled_fields[-1], self.target_fields, self.tactician.weight)
            ortho.plotField(self.scheduled_fields[-1],self.scheduled_fields,options_basemap=dict(date='2017/02/20 05:00:00'))
            raw_input('WAIT')
            import pdb; pdb.set_trace()
            raise Exception()

        return field_select


    def run(self, tstart=None, tstop=None, clip=False, plot=False, mode='coverage'):
        """
        Schedule a chunk of exposures. This is the loop where date is incremented

        Parameters:
        -----------
        tstart : Chunk start time
        tstop  : Chunk end time (may be replace with chunk length)
        plot   : Plot the chunk (may be removed)

        Returns:
        --------
        fields : Scheduled fields
        """
        # Reset the scheduled fields
        self.scheduled_fields = self.FieldType()

        # If no tstop, run for 90 minutes
        timedelta = 90*ephem.minute
        if tstart is None: tstart = ephem.now()
        if tstop is None: tstop = tstart + timedelta

        # Convert strings into dates
        if isinstance(tstart,basestring):
            tstart = ephem.Date(tstart)
        if isinstance(tstop,basestring):
            tstop = ephem.Date(tstop)

        msg  = "Run start: %s\n"%datestr(tstart,4)
        msg += "Run end: %s\n"%datestr(tstop,4)
        msg += "Run time: %s minutes"%(timedelta/ephem.minute)
        logging.debug(msg)

        msg = "Previously completed fields: %i"%len(self.completed_fields)
        logging.info(msg)

        msg = "Scheduling with tactician: %s"%mode
        logging.info(msg)

        date = tstart
        latch = True
        while latch:
            logging.debug(' '+datestr(date,4))

            # Check to see if in valid observation window
            if self.windows is not None:
                inside = False
                for window in self.windows:
                    if date >= window[0] and date < window[-1]:
                        inside = True
                        break

                if not inside:
                    if clip:
                        break
                    else:
                        msg = 'Date outside of nominal observing windows'
                        logging.warning(msg)

            # Select one (or more) fields from the tactician
            field_select = self.select_field(date, mode)

            # Now update the time from the selected field
            date = ephem.Date(field_select[-1]['DATE']) + constants.FIELDTIME

            self.completed_fields = self.completed_fields + field_select
            self.scheduled_fields = self.scheduled_fields + field_select

            msg=" %(DATE).19s: id=%(ID)10s, secz=%(AIRMASS).2f, slew=%(SLEW).2f"
            msg+=", moon=%(PHASE).0f%%,%(ALT).0fdeg"
            for i,f in zip(field_select.unique_id,field_select):
                params = dict([('ID',i)]+[(k,f[k]) for k in f.dtype.names])
                params.update({'PHASE':self.tactician.moon.phase,"ALT":np.degrees(self.tactician.moon.alt)})
                logging.info(msg%params)

            #if plot: self.plotField(date, field_select)
            if plot:
                ortho.plotField(field_select[:-1],self.target_fields,self.completed_fields)
            if date >= tstop: break

        msg = "Newly scheduled fields: %i"%len(self.scheduled_fields)
        logging.info(msg)

        return self.scheduled_fields

    def schedule_field(self, hex, tiling, band=None, date=None, plot=False, mode=None):
        """
        Schedule a single filed at a given time.

        Parameters:
        -----------
        hexid  : the hex ID of the field
        tiling : the tiling number of the field
        band   : The band of the field
        date   : The date/time for observation
        plot   : Plot the output
        mode   : Mode for scheduler tactician

        Returns:
        --------
        field : The scheduled field
        """
        # Probably cleaner to make this it's own tactician
        date = ephem.Date(date) if date else ephem.now()

        select  = (self.target_fields['HEX']==hex)
        select &= (self.target_fields['TILING']==tiling)
        if band is not None:
            select &= (self.target_fields['FILTER']==band)
        index = np.nonzero(select)[0]

        field = self.target_fields[select]
        nfields = select.sum()
        field['DATE'] = map(datestring,nfields*[date])
        return field

    def schedule_chunk(self,tstart=None,chunk=60,clip=False,plot=False,mode=None):
        """
        Schedule a chunk of exposures.

        Parameters:
        -----------
        tstart : Start time (UTC); in `None` use `ephem.now()`
        chunk  : Chunk of time to schedule.
        plot   : Dynamically plot each scheduled exposure
        mode   : Mode for scheduler tactician

        Returns:
        --------
        fields : Scheduled fields
        """
        # If no tstop, run for 90 minutes
        if tstart is None: tstart = ephem.now()
        tstop = tstart + chunk*ephem.minute

        return self.run(tstart,tstop,clip,plot,mode)

    def schedule_nite(self,date=None,chunk=60,clip=False,plot=False,mode=None):
        """
        Schedule a night of observing.

        A `nite` is defined by the day (UTC) at noon local time before
        observing started.

        Parameters:
        -----------
        date  : The date of the nite to schedule
        chunk : The duration of a chunk of exposures (minutes)
        plot  : Dynamically plot the progress after each chunk
        mode  : Mode for scheduler tactician

        Returns:
        --------
        chunks : A list of the chunks generated for the scheduled nite.
        """

        # Create the nite
        nite = get_nite(date)

        # Convert chunk to MJD
        if chunk > 1: chunk = chunk*ephem.minute

        try:
            nites = [get_nite(w[0]) for w in self.windows]
            idx = nites.index(nite)
            start,finish = self.windows[idx]
        except (TypeError, ValueError):
            msg = "Requested nite (%s) not found in windows:\n"%nite
            msg += '['+', '.join([n for n in nites])+']'
            logging.warning(msg)

            start = date
            self.observatory.date = date
            self.observatory.horizon = self.observatory.twilight
            finish = self.observatory.next_rising(ephem.Sun(), use_center=True)
            self.observatory.horizon = '0'

            logging.info("Night start (UTC):  %s"%datestr(start))
            logging.info("Night finish (UTC): %s"%datestr(finish))

        chunks = []
        i = 0
        while start < finish:
            i+=1
            msg = "Scheduling %s -- Chunk %i"%(start,i)
            logging.debug(msg)
            end = start+chunk
            scheduled_fields = self.run(start,end,clip=clip,plot=False,mode=mode)

            if plot:
                field_select = scheduled_fields[-1:]
                ortho.plotField(field_select,self.target_fields,self.completed_fields)
                if (raw_input(' ...continue ([y]/n)').lower()=='n'):
                    break

            chunks.append(scheduled_fields)
            start = ephem.Date(chunks[-1]['DATE'][-1]) + constants.FIELDTIME
            #start = end

        if plot: raw_input(' ...finish... ')

        return chunks

    def schedule_survey(self,start=None,end=None,chunk=60,plot=False,mode=None):
        """
        Schedule the entire survey.

        Parameters:
        -----------
        start : Start of survey (int or str)
        end   : End of survey (int or str)
        chunk : The duration of a chunk of exposures (minutes)
        plot  : Dynamically plot the progress after each night
        mode  : Mode of scheduler tactician

        Returns:
        --------
        scheduled_nites : An ordered dictionary of scheduled nites
        """

        self.scheduled_nites = odict()

        for tstart,tend in self.windows:
            if start is not None and ephem.Date(tstart) < ephem.Date(start):
                continue
            if end is not None and ephem.Date(tend) > ephem.Date(end):
                continue

            #nite = nitestring(tstart)
            nite = get_nite(tstart)

            try:
                chunks = self.schedule_nite(tstart,chunk,clip=True,plot=False,mode=mode)
            except ValueError as error:
                ortho.plotField(self.completed_fields[-1:],self.target_fields,
                                self.completed_fields)
                raise(error)

            self.scheduled_nites[nite] = chunks

            if plot:
                ortho.plotField(self.completed_fields[-1:],self.target_fields,self.completed_fields)#,options_basemap=dict(date='2017/02/21 05:00:00'))

                if (raw_input(' ...continue ([y]/n)').lower()=='n'):
                    break

        if plot: raw_input(' ...finish... ')
        return self.scheduled_nites

    def write(self,filename):
        self.scheduled_fields.write(filename)

    @classmethod
    def common_parser(cls):
        """
        Comman argument parser for scheduler tools.
        """
        from obztak.utils.parser import Parser, DatetimeAction

        description = __doc__
        parser = Parser(description=description)
        #parser.add_argument('--survey',choices=['obztak','maglites','bliss'],
        #                    default = None, help='choose survey to schedule.')
        parser.add_argument('-p','--plot',action='store_true',
                            help='create visual output.')
        parser.add_argument('--utc','--utc-start',dest='utc_start',action=DatetimeAction,
                            help="start time for observation.")
        parser.add_argument('--utc-end',action=DatetimeAction,
                            help="end time for observation.")
        parser.add_argument('-k','--chunk', default=60., type=float,
                            help = 'time chunk')
        parser.add_argument('-f','--fields',default=None,
                            help='all target fields.')
        parser.add_argument('-m','--mode',default='coverage',
                            help='Mode for scheduler tactician.')
        parser.add_argument('-w','--windows',default=None,
                            help='observation windows.')
        parser.add_argument('-c','--complete',nargs='?',action='append',
                            help="fields that have been completed.")
        parser.add_argument('-o','--outfile',default=None,
                            help='save output file of scheduled fields.')
        parser.add_argument('--write-protect',action='store_true',
                            help='write-protect output files')
        return parser

    @classmethod
    def parser(cls):
        return cls.common_parser()

    @classmethod
    def main(cls):
        args = cls.parser().parse_args()
        scheduler = cls(args.fields,args.windows,args.complete)
        scheduler.run(tstart=args.utc_start,tstop=args.utc_end,plot=args.plot)
        if args.outfile: 
            scheduler.scheduled_fields.write(args.outfile)
         
        return scheduler

############################################################

if __name__ == '__main__':
    scheduler = Scheduler.main()

############################################################
