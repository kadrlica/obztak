"""
Module providing the survey scheduler.
"""

import os,sys
import copy
import numpy as np
import time
import ephem
import logging
from collections import OrderedDict as odict

from obztak.utils import constants
from obztak.utils import ortho
from obztak.utils import fileio

from obztak.ctio import CTIO
from obztak.field import FieldArray
from obztak.tactician import CoverageTactician
from obztak.utils.date import get_nite, isstring
from obztak.utils.date import datestr, datestring, nitestring, utc2nite
from obztak.factory import tactician_factory

# For debugging (use the verbose command line argument)
#logging.basicConfig(level=20) # KCB

############################################################

class Scheduler(object):
    """
    Deal with survey scheduling.
    """
    _defaults = odict([
        #('tactician','coverage'),
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

        self.create_seeing()

    def create_seeing(self,filename=None,mode='qc'):
        import obztak.seeing
        #dirname ='/Users/kadrlica/delve/observing/data/'
        #basename = 'delve_sim_01.csv.gz'
        #filename = os.path.join(dirname,basename)
        if mode == 'dimm':
            self.seeing = obztak.seeing.DimmSeeing(filename=filename)
        elif mode == 'qc':
            self.seeing = obztak.seeing.QcSeeing(filename=filename)
        else:
            self.seeing = obztak.seeing.QcSeeing(filename=filename)

        return self.seeing

    def load_target_fields(self, target_fields=None):
        if target_fields is None:
            target_fields = self._defaults['targets']

        if isstring(target_fields):
            self.target_fields = self.FieldType.read(target_fields)
            logging.info("Loading target fields...\n %s"%target_fields)
        else:
            self.target_fields = self.FieldType(target_fields)
        return self.target_fields

    def load_windows(self, windows=None):
        """
        Load the set of start and stop times for the observation windows.
        """
        if windows is None:
            windows = self._defaults['windows']

        if isstring(windows):
            logging.info("Loading observing windows...\n %s"%windows)
            windows = fileio.csv2rec(windows)

        self.windows = []
        for start,end in windows:
            self.windows.append([ephem.Date(start), ephem.Date(end)])

        # Sanity check that observation windows are properly sorted
        for ii,(start,end) in enumerate(self.windows):
            msg = 'Observation windows are not properly sorted\n'
            msg+= '%s: %s -- %s'%(get_nite(start),datestr(start),datestr(end))
            if (end < start):
                logging.warn(msg)
            if ii > 0 and (start < self.windows[ii-1][1]):
                logging.warn(msg)

        logging.debug('Observation Windows:')
        for start,end in self.windows:
            logging.debug('  %s: %s UTC -- %s UTC'%(get_nite(start),datestr(start),datestr(end)))
        logging.debug(30*'-')

    def load_observed_fields(self):
        """
        Load fields from the telemetry database that were already observed.
        """
        logging.info("Loading observed fields...")
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
        logging.info("Loading completed fields...")
        if isinstance(completed_fields,list):
            if completed_fields[0].lower()=='none':
                self.completed_fields = self.FieldType()
                return self.completed_fields
        elif isstring(completed_fields):
            if completed_fields.lower()=='none':
                self.completed_fields = self.FieldType()
                return self.completed_fields

        self.completed_fields = copy.deepcopy(self.observed_fields)

        if not completed_fields:
            return self.completed_fields

        if isstring(completed_fields):
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

    def create_tactician(self, cls=None, mode=None):
        """ Create a tactician in the given mode.

        Parameters:
        -----------
        cls : the tactician class [defaults to survey]
        mode: the tactician mode

        Returns:
        --------
        tac : the tactician
        """
        return tactician_factory(cls=cls, mode=mode)

    def select_field(self, date, mode=None):
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

        # ADW: Why do we create the tactician each time?
        self.tactician = self.create_tactician(mode=mode)
        self.tactician.set_date(date)
        self.tactician.set_target_fields(self.target_fields[sel])
        self.tactician.set_completed_fields(self.completed_fields)
        self.tactician.fwhm = self.fwhm

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
            #ortho.plotField(self.scheduled_fields[-1],self.scheduled_fields,options_basemap=dict(date='2017/02/20 05:00:00'))
            raw_input('WAIT')
            import pdb; pdb.set_trace()
            raise Exception()

        return field_select


    def run(self, tstart=None, tstop=None, clip=False, plot=False, mode=None):
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
        if tstart is None: tstart = ephem.now()
        if tstop is None:
            timedelta = 90*ephem.minute
            tstop = tstart + timedelta

        # Convert strings into dates
        if isstring(tstart):
            tstart = ephem.Date(tstart)
        if isstring(tstop):
            tstop = ephem.Date(tstop)

        msg  = "\nRun start: %s\n"%datestr(tstart,4)
        msg +=   "Run stop : %s\n"%datestr(tstop,4)
        logging.debug(msg)

        msg = "Previously completed fields: %i"%len(self.completed_fields)
        logging.info(msg)

        # This is not safe since tactician is re-created in select_field
        self.tactician = self.create_tactician(mode=mode)
        msg = "Scheduling with '%s' in mode '%s'"%(self.tactician.__class__.__name__,self.tactician.mode)
        logging.info(msg)

        self.seeing.set_date(datestr(tstart))
        self.fwhm = self.seeing.get_fwhm(band='i',airmass=1.0)
        logging.info("Predicted i-band zenith fwhm: %.2f arcsec"%self.fwhm)
        logging.debug(self.seeing.raw)

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
            try:
                field_select = self.select_field(date, mode)
            except Exception as e:
                # Only write if error occurred outside observing window
                if not inside:
                    logging.warning(str(e))
                    break
                else:
                    raise(e)

            # Now update the time from the last selected field (note duplication in tactician.select_field)
            fieldtime = field_select[-1]['EXPTIME']*ephem.second + constants.OVERHEAD
            date = ephem.Date(field_select[-1]['DATE'].astype(str)) + fieldtime

            self.completed_fields = self.completed_fields + field_select
            self.scheduled_fields = self.scheduled_fields + field_select

            #msg=" %(DATE).19s: id=%(ID)10s, secz=%(AIRMASS).2f, slew=%(SLEW).2f"
            msg=" %(DATE).19s: id=%(ID)10s, ra,dec=%(RA).2f,%(DEC).2f"
            msg+=", secz=%(AIRMASS).2f, moon=%(PHASE).0f%%,%(ALT).0fdeg"
            for i,f in zip(field_select.unique_id,field_select):
                params = dict([('ID',i)]+[(k,f[k]) for k in f.dtype.names])
                params.update({
                    'PHASE':self.tactician.moon.phase,
                    "ALT":np.degrees(self.tactician.moon.alt),
                    'DATE':params['DATE'].astype(str),
                    'FILTER':params['FILTER'].astype(str)
                })
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
        field['DATE'] = list(map(datestring,nfields*[date]))
        return field

    def schedule_chunk(self, tstart=None, chunk=60, clip=False,
                       plot=False, mode=None):
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

    def schedule_nite(self, date=None, start=None, end=None,
                      chunk=60, clip=False, plot=False, mode=None):
        """
        Schedule a night of observing.

        A `nite` is defined by the day (UTC) at noon local time before
        observing started.

        Parameters:
        -----------
        date  : The date of the nite to schedule
        start : When to start the observations (UTC)
        end   : When to end the observations (UTC)
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
            window_start,window_end = self.windows[idx]

            if start is None:
                start = window_start
            else:
                logging.warn("Over-writing nite start time.")

            if end is None:
                end = window_end
            else:
                logging.warn("Over-writing nite end time.")

        except (TypeError, ValueError):
            msg = "Requested nite (%s) not found in windows"%nite
            logging.warning(msg)
            msg = '['+', '.join([n for n in nites])+']'
            logging.debug(msg)

            self.observatory.date = date
            self.observatory.horizon = self.observatory.twilight
            if start is None:
                start = date
            if end is None:
                end = self.observatory.next_rising(ephem.Sun(), use_center=True)
            self.observatory.horizon = '0'

        logging.info("Night start (UTC):  %s"%datestr(start))
        logging.info("Night finish (UTC): %s"%datestr(end))

        chunk_start, chunks = start, []
        while chunk_start < end:
            msg = "Scheduling %s -- Chunk %i"%(chunk_start, len(chunks)+1)
            logging.debug(msg)
            chunk_end = chunk_start+chunk

            try:
                scheduled_fields = self.run(chunk_start,chunk_end,
                                            clip=clip,plot=False,mode=mode)
            except ValueError as e:
                # Write fields even if there is an error
                chunks.append(self.scheduled_fields)
                logging.warning(str(e))
                break

            if len(scheduled_fields) == 0:
                # No new fields scheduled (probably error)
                logging.warning("No new fields scheduled.")
                break

            chunks.append(scheduled_fields)
            fieldtime = chunks[-1]['EXPTIME'][-1]*ephem.second + constants.OVERHEAD
            chunk_start = ephem.Date(chunks[-1]['DATE'][-1].astype(str)) + fieldtime
            #start = end

        if plot: raw_input(' ...enter to finish... ')

        return chunks

    def schedule_survey(self, start=None, end=None, chunk=60,
                        plot=False, mode=None, write=False, dirname=None):
        """
        Schedule the entire survey.

        Parameters:
        -----------
        start : Start date of survey (int or str)
        end   : End date of survey (int or str)
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

            if write:
                self.write_nite(nite,chunks,dirname=dirname)

            if plot:
                ortho.plotField(self.completed_fields[-1:],self.target_fields,
                                self.completed_fields)
                if (raw_input(' ...continue ([y]/n)').lower()=='n'):
                    import pdb; pdb.set_trace()

        if plot: raw_input(' ...finish... ')
        return self.scheduled_nites

    def write_nite(self,nite,chunks,dirname=None):
        if dirname:
            outdir = os.path.join(dirname,nite)
        else:
            outdir = os.path.join(nite)
        if not os.path.exists(outdir): os.makedirs(outdir)
        outfile = os.path.join(outdir,nite+'.json')
        base,ext = os.path.splitext(outfile)

        for i,chunk in enumerate(chunks):
            if len(chunks) > 1:
                outfile = base+'_%02d'%(i+1)+ext
            logging.debug("Writing %s..."%outfile)
            chunk.write(outfile)


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
                            help = 'time chunk (minutes)')
        parser.add_argument('-f','--fields',default=None,
                            help='all target fields.')
        #parser.add_argument('-m','--mode',default='coverage',
        #                    help='Mode for scheduler tactician.')
        parser.add_argument('-m','--mode',default=None,
                            help='mode for scheduler tactician.')
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
