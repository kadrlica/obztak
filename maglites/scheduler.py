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

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.constants as constants
import maglites.utils.ortho
import maglites.utils.ortho as ortho
import maglites.utils.fileio as fileio

from maglites.field import FieldArray
from maglites.utils.ortho import get_nite, datestring

CONDITIONS = odict([
    ('great', [1.4, 2.0]),
    ('good',  [0.0, 2.0]),
    ('fine',  [0.0, 1.9]),
    ('ok',    [0.0, 1.6]),
    ('poor',  [0.0, 1.5]),
    ('bad',   [0.0, 1.4]),
])

# For debugging (can also use the verbose command line argument)
#logging.basicConfig(level=20) # KCB

############################################################

class Scheduler(object):
    """
    Deal with survey scheduling.
    """

    def __init__(self,target_fields=None,observation_windows=None,completed_fields=None):
        self.loadTargetFields(target_fields)
        self.loadObservationWindows(observation_windows)
        self.loadObservedFields()
        self.loadCompletedFields(completed_fields)
        
        self.scheduled_fields = FieldArray()

        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        self.loadBlancoConstraints()

    def loadTargetFields(self, target_fields=None):
        if target_fields is None:
            datadir = fileio.get_datadir()
            target_fields = os.path.join(datadir,"maglites-target-fields.csv")
        
        if isinstance(target_fields,basestring):
            self.target_fields = FieldArray.read(target_fields)
        else:
            self.target_fields = target_fields
        
    def loadObservationWindows(self, observation_windows=None):
        """
        Load the set of start and stop times for the observation windows.
        """
        if observation_windows is None: 
            datadir = fileio.get_datadir()
            observation_windows = os.path.join(datadir,"maglites-windows.csv")
            logging.info("Setting default observing windows: %s"%observation_windows)
            

        if isinstance(observation_windows,basestring):
            observation_windows = fileio.csv2rec(observation_windows)
            
        self.observation_windows = []
        for start,end in observation_windows:
            self.observation_windows.append([ephem.Date(start), ephem.Date(end)])

        # Do a sanity check to make sure that observation windows are properly sorted
        for ii in range(0, len(self.observation_windows)):
            if self.observation_windows[ii][1] < self.observation_windows[ii][0]:
                logging.warning('Observation windows are not properly sorted')
                logging.info('%s -- %s'%(self.observation_windows[ii][0],self.observation_windows[ii][1]))
            if ii > 0:
                if self.observation_windows[ii][0] < self.observation_windows[ii - 1][1]:
                    logging.warning('Observation windows are not properly sorted')
                    logging.info('%s -- %s'%(self.observation_windows[ii][0],self.observation_windows[ii][1]))

        logging.info('Observation Windows:')
        for start,end in self.observation_windows:
            logging.info('  %s UTC -- %s UTC'%(start,end))
        logging.info(30*'-')

    def loadObservedFields(self, **kwargs):
        """
        Load observed fields from the telemetry database.
        """
        try: 
            fields = FieldArray.load_database()
        except Exception as e: 
            logging.warning("Failed to load completed exposures from database")
            logging.info(e)
            fields = FieldArray()
        self.observed_fields = fields
        return self.observed_fields


    def loadCompletedFields(self, completed_fields=None):
        """Load completed fields. The default behavior is to load the
        observed fields as completed fields. However, if the string
        'None' is passed then return an empty FieldArray.

        Parameters:
        -----------
        completed_fields : Filename, list of filenames, or FieldArray object.

        Returns:
        --------
        FieldArray of the completed fields

        """
        # Deal with 'None' string
        if isinstance(completed_fields,list):
            if completed_fields[0].lower()=='none':
                self.completed_fields = FieldArray()
                return self.completed_fields
        elif isinstance(completed_fields,basestring):
            if completed_fields.lower()=='none':
                self.completed_fields = FieldArray()
                return self.completed_fields

        self.completed_fields = copy.deepcopy(self.observed_fields)

        if not completed_fields:
            return self.completed_fields

        if isinstance(completed_fields,basestring):
            completed_fields = [completed_fields]

        if isinstance(completed_fields,list):
            fields = FieldArray()
            for filename in completed_fields:
                fields = fields + FieldArray.read(filename)
        
            completed_fields = fields

        new=~np.in1d(completed_fields.unique_id,self.completed_fields.unique_id)
        new_fields = completed_fields[new]
        self.completed_fields = self.completed_fields + new_fields
        return self.completed_fields

    def loadBlancoConstraints(self):
        """
        Load telescope pointing constraints
        """
        # Updated to remove the dependence on scipy (which is broken on the mountain)
        datadir = fileio.get_datadir()
        data = np.recfromtxt(os.path.join(datadir,'blanco_hour_angle_limits.dat'), names=True)

        self.blanco_constraints = data
        ha_degrees = np.tile(0., len(self.blanco_constraints['HA']))
        for ii in range(0, len(self.blanco_constraints['HA'])):
            ha_degrees[ii] = maglites.utils.projector.hms2dec(self.blanco_constraints['HA'][ii])
        
        ha_degrees -= 1.25 # Buffer to protect us from the chicken

        self.f_hour_angle_limit = lambda dec: np.interp(dec,self.blanco_constraints['Dec'], ha_degrees, left=-1, right=-1)
        self.f_airmass_limit = lambda dec: np.interp(dec,self.blanco_constraints['Dec'], self.blanco_constraints['AirmassLimit'], left=-1, right=-1)

        return self.f_hour_angle_limit,self.f_airmass_limit

    def selectField(self, date, ra_previous=None, dec_previous=None, plot=False, mode='coverage'):
        """
        Select the `best` field to observe at a given time.

        A single field can contain multiple exposures (for example g- and r-band).
        
        Available modes:
        `balance`  : 
        `balance2` : 
        `balance3` : 

        Parameters:
        -----------
        date         : The time to schedule the exposure
        ra_previous  : The ra of the previous exposure
        dec_previous : The dec of the previous exposure
        plot         : Plot the output
        mode         : Algorithm used to select the exposure

        Returns:
        --------
        field        :  The selected exposures as a FieldArray object
        """

        self.observatory.date = ephem.Date(date)

        ra_zenith, dec_zenith = self.observatory.radec_of(0, '90') # RA and Dec of zenith
        ra_zenith = np.degrees(ra_zenith)
        dec_zenith = np.degrees(dec_zenith)
        airmass = maglites.utils.projector.airmass(ra_zenith, dec_zenith, self.target_fields['RA'], self.target_fields['DEC'])
        airmass_next = maglites.utils.projector.airmass(ra_zenith + 15., dec_zenith, self.target_fields['RA'], self.target_fields['DEC'])

        # Include moon angle
        moon = ephem.Moon()
        moon.compute(date)
        ra_moon = np.degrees(moon.ra)
        dec_moon = np.degrees(moon.dec)
        moon_angle = maglites.utils.projector.angsep(ra_moon, dec_moon, self.target_fields['RA'], self.target_fields['DEC'])

        # Slew from the previous pointing
        if ra_previous is not None and dec_previous is not None:
            slew = maglites.utils.projector.angsep(ra_previous, dec_previous, self.target_fields['RA'], self.target_fields['DEC'])
            slew_ra = np.fabs(ra_previous - self.target_fields['RA'])
            slew_dec = np.fabs(dec_previous - self.target_fields['DEC'])
        else:
            slew = np.tile(0., len(self.target_fields['RA']))
            slew_ra = np.tile(0., len(self.target_fields['RA']))
            slew_dec = np.tile(0., len(self.target_fields['RA']))

        # Hour angle restrictions
        #hour_angle_degree = copy.copy(self.target_fields['RA']) - ra_zenith # BUG
        #hour_angle_degree[hour_angle_degree > 180.] = hour_angle_degree[hour_angle_degree > 180.] - 360. # BUG
        hour_angle_degree = copy.copy(self.target_fields['RA']) - ra_zenith
        hour_angle_degree[hour_angle_degree < -180.] += 360.
        hour_angle_degree[hour_angle_degree > 180.] -= 360.
        cut_hour_angle = np.fabs(hour_angle_degree) < self.f_hour_angle_limit(self.target_fields['DEC']) # Check the hour angle restrictions at south pole
        
        # Airmass restrictions
        cut_airmass = airmass < self.f_airmass_limit(self.target_fields['DEC'])

        # Declination restrictions
        cut_declination = self.target_fields['DEC'] > maglites.utils.constants.SOUTHERN_REACH

        # Don't consider fields which have already been observed
        cut_todo = np.logical_not(np.in1d(self.target_fields['ID'], self.completed_fields['ID']))
        cut = cut_todo & cut_hour_angle & cut_airmass & cut_declination & (airmass < 2.) # Now with Blanco telescope constraints
        #cut = cut_todo & (airmass < 2.) # Original

        # Exclude special fields unless using special tacticians
        if mode not in ['smcnod']:
            cut = cut & (self.target_fields['PRIORITY'] < 90)

        # Need to figure out what to do if there are no available fields...

        # Now apply some kind of selection criteria, e.g., 
        # select the field with the lowest airmass
        #airmass[np.logical_not(cut)] = 999.

        if mode == 'airmass':
            airmass_effective = copy.copy(airmass)
            airmass_effective[np.logical_not(cut)] = np.inf # Do not observe fields that are unavailable
            airmass_effective += self.target_fields['TILING'] # Priorize coverage over multiple tilings
            index_select = np.argmin(airmass_effective)
        elif mode == 'ra':
            # Different selection
            #ra_effective = copy.copy(self.target_fields['RA'])
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = np.inf
            ra_effective += 360. * self.target_fields['TILING']
            index_select = np.argmin(ra_effective)
        elif mode == 'slew':
            #ra_effective = copy.copy(self.target_fields['RA'])
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = np.inf
            ra_effective += 360. * self.target_fields['TILING']
            ra_effective += slew**2
            #ra_effective += 2. * slew
            index_select = np.argmin(ra_effective)
        elif mode == 'balance':
            """
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = np.inf
            ra_effective += 360. * self.target_fields['TILING']
            #ra_effective += 720. * self.target_fields['TILING']
            ra_effective += slew**2
            ra_effective += 100. * (airmass - 1.)**3
            weight = ra_effective
            index_select = np.argmin(weight)
            weight = hour_angle_degree
            """
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = np.inf
            weight += 3. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'balance2':
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = np.inf
            weight += 360. * self.target_fields['TILING']
            weight += slew_ra**2
            weight += slew_dec
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'balance3':
            logging.debug("Slew: %s"%slew)
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = np.inf
            weight += 3. * 360. * self.target_fields['TILING']
            """
            x_slew, y_slew = zip(*[[0., 0.],
                                   [2.5, 10.],
                                   [5., 30.],
                                   [10., 150.],
                                   [20., 250.],
                                   [50., 500.],
                                   [180., 5000.]])
            """
            x_slew, y_slew = zip(*[[0., 0.],
                                   [2.5, 10.],
                                   [5., 30.],
                                   [10., 500.], # 
                                   [20., 1000.], # 500
                                   [50., 5000.], # 1000
                                   [180., 5000.]])
            weight += np.interp(slew, x_slew, y_slew, left=np.inf, right=np.inf)
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'airmass2':
            weight = 200. * (airmass - airmass_next)
            weight[np.logical_not(cut)] = np.inf
            weight += 360. * self.target_fields['TILING']
            weight += 100. * (airmass - 1.)**3
            weight += slew**2
            index_select = np.argmin(weight)
        elif mode in ('coverage','good'):
            weight = copy.copy(hour_angle_degree)
            #weight[np.logical_not(cut)] = 9999.
            weight[np.logical_not(cut)] = np.inf
            weight += 6. * 360. * self.target_fields['TILING'] # Was 6, 60
            weight += slew**3 # slew**2
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'coverage2':
            weight = copy.copy(hour_angle_degree)
            weight *= 2.
            weight[np.logical_not(cut)] = np.inf
            weight += 6. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'coverage3':
            weight = copy.copy(hour_angle_degree)
            weight *= 0.5
            weight[np.logical_not(cut)] = np.inf
            weight += 6. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'lowairmass':
            weight = 2.0 * copy.copy(hour_angle_degree)
            #if len(self.scheduled_fields) == 0:
            #    weight += 200. * maglites.utils.projector.angsep(self.target_fields['RA'], 
            #                                                     self.target_fields['DEC'],
            #                                                     90., -70.)
            weight[np.logical_not(cut)] = np.inf
            weight += 3. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            #weight += 2000. * (airmass - 1.)**3 # 200
            weight += 5000. * (airmass > 1.5)
            index_select = np.argmin(weight)
            
            """
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = np.inf
            weight += 3. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            weight += 1000. * (airmass - 1.)**3
            index_select = np.argmin(weight)
            """
        elif mode in CONDITIONS.keys():
            weight = 2.0 * copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = np.inf
            weight += 3. * 360. * self.target_fields['TILING']
            weight += slew**3
            airmass_min, airmass_max = CONDITIONS[mode]
            airmass_sel = ((airmass < airmass_min) | (airmass > airmass_max))
            # ADW: This should probably also be in there
            weight += 100. * (airmass - 1.)**3
            weight += 5000. * airmass_sel
            index_select = np.argmin(weight)
        elif mode == 'smcnod':
            weight = 10000. * np.logical_not(np.in1d(self.target_fields['HEX'], maglites.utils.constants.HEX_SMCNOD)).astype(float)
            weight[np.logical_not(cut)] = np.inf
            weight += 360. * self.target_fields['TILING']
            weight += slew
            index_select = np.argmin(weight)
        else:
            msg = "Unrecognized mode: %s"%mode
            raise Exception(msg)

        # Search for other exposures in the same field
        field_id = self.target_fields['HEX'][index_select]
        tiling = self.target_fields['TILING'][index_select]        

        index = np.nonzero( (self.target_fields['HEX']==field_id) & \
                                   (self.target_fields['TILING']==tiling) & cut)[0]

        
        timedelta = constants.FIELDTIME*np.arange(len(index))
        if np.any(slew[index] > 5.):
            # Apply a 30 second penalty for slews over 5 deg.
            # This is not completely realistic, but better than nothing
            # This is also broken when selecting two fields at once
            timedelta += 30*ephem.second
        field_select = self.target_fields[index]
        #print 'CHECKPOINT 1', field_select, index, field_id, tiling # KCB
        #print 'CHECKPOINT 2', np.any((self.target_fields['HEX']==field_id) & (self.target_fields['TILING']==tiling)), np.sum(cut), np.min(weight)
        #print 'CHECKPOINT 3', np.nonzero( (self.target_fields['HEX']==field_id) & \
        #                           (self.target_fields['TILING']==tiling) & cut)
        field_select['AIRMASS'] = airmass[index]
        field_select['DATE'] = map(datestring,date+timedelta)
        field_select['SLEW'] = slew[index]
        field_select['MOONANGLE'] = moon_angle[index]
        field_select['HOURANGLE'] = hour_angle_degree[index]

        msg = str(field_select)
        logging.debug(msg)

        # For diagnostic purposes
        if False and len(self.scheduled_fields) % 10 == 0:
            ortho.plotWeight(field_select[-1], self.target_fields, weight)
            raw_input('WAIT')

        if len(field_select) == 0:
            msg = 'No field selected... now we\'ve got problems'
            logging.warning(msg)
            print field_id, tiling
            print index_select
            print cut[index_select]
            print index
            print cut.sum()
            print weight
            #ortho.plotWeight(field_select, self.target_fields, weight)
            raw_input('WAIT')
            

        return field_select


    def run(self, tstart=None, tstop=None, clip=False, plot=True, mode=None):
        """
        Schedule a chunk of exposures.
        
        Parameters:
        -----------
        tstart : Chunk start time
        tstop  : Chunk end time (may be replace with chunk length)
        plot   : Plot the chunk (may be removed)
        
        Returns:
        --------
        fields : Scheduled fields
        """
        if mode is None: mode='coverage'

        # Reset the scheduled fields
        self.scheduled_fields = FieldArray(0)

        # If no tstop, run for 90 minutes
        timedelta = 90*ephem.minute
        if tstart is None: tstart = ephem.now()
        if tstop is None: tstop = tstart + timedelta
        msg  = "Run start: %s\n"%datestring(tstart)
        msg += "Run end: %s\n"%datestring(tstop)
        msg += "Run time: %s minutes"%(timedelta/ephem.minute)
        logging.debug(msg)

        # Convert strings into dates
        if isinstance(tstart,basestring):
            tstart = ephem.Date(tstart)
        if isinstance(tstop,basestring):
            tstop = ephem.Date(tstop)

        msg = "Previously completed fields: %i"%len(self.completed_fields)
        logging.info(msg)

        msg = "Scheduling with tactician mode: %s"%mode
        logging.info(msg)

        date = tstart
        latch = True
        while latch:
            logging.debug('  '+datestring(date))

            # Check to see if in valid observation window
            if self.observation_windows is not None:
                inside = False
                for window in self.observation_windows:
                    if date >= window[0] and date < window[-1]: 
                        inside = True 

                if not inside:
                    if clip: 
                        break
                    else:
                        msg = 'Date outside of nominal observing windows'
                        logging.warning(msg)

                
            # FIXME: I think that ra_previous and dec_previous don't need to be passed
            compute_slew = True
            if len(self.completed_fields) == 0:
                compute_slew = False
            else:
                if (date - ephem.Date(self.completed_fields['DATE'][-1])) > (30. * ephem.minute):
                    compute_slew = False

                
            if compute_slew:
                field_select = self.selectField(date, ra_previous=self.completed_fields['RA'][-1], dec_previous=self.completed_fields['DEC'][-1], plot=plot,mode=mode)
            else:
                field_select = self.selectField(date, plot=plot, mode=mode)

            id_select = field_select['ID']
            # Previously, increment time by a constant
            #date = date + len(field_select)*constants.FIELDTIME
            # Now update the time from the selected field
            date = ephem.Date(field_select[-1]['DATE']) + constants.FIELDTIME

            self.completed_fields = self.completed_fields + field_select
            self.scheduled_fields = self.scheduled_fields + field_select

            msg = "  %(DATE).19s: id=%(ID)s, airmass=%(AIRMASS).2f, slew=%(SLEW).2f"
            for i,f in zip(field_select.unique_id,field_select):
                params = dict([('ID',i)]+[(k,f[k]) for k in f.dtype.names])
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
        date = ephem.Date(date) if date else ephem.now()

        select  = (self.target_fields['HEX']==hex)
        select &= (self.target_fields['TILING']==tiling)
        if band is not None:
            select &= (self.target_fields['FILTER']==band)
        index = np.nonzero(select)[0]

        field = self.target_fields[select]

        field['DATE'] = map(datestring,select.sum()*[date])
        #field['AIRMASS'] = 
        #field['DATE'] = 
        #field['SLEW'] = 
        #field['MOONANGLE'] = 
        #field['HOURANGLE'] = 
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

        A `nite` is defined by the day (UTC) at noon local time before observing started.

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
        nite_tuple = nite.tuple()[:3]

        # Convert chunk to MJD
        if chunk > 1: chunk = chunk*ephem.minute

        try:
            nites = [get_nite(w[0]) for w in self.observation_windows]
            nite_tuples = [n.tuple()[:3] for n in nites]
            idx = nite_tuples.index(nite_tuple)
            start,finish = self.observation_windows[idx]
        except (TypeError, ValueError):
            msg = "Requested nite not found in windows:\n"
            msg += "%s/%s/%s : "%nite_tuple
            msg += '['+', '.join(['%s/%s/%s'%t for t in nite_tuples])+']'
            logging.warning(msg)

            # WARNING: copy.deepcopy doesn't work for ephem.Observer
            start = date
            self.observatory.date = date
            self.observatory.horizon = '-14'
            finish = self.observatory.next_rising(ephem.Sun(), use_center=True)
            self.observatory.horizon = '0'

            logging.info("Night start time: %s"%datestring(start))
            logging.info("Night finish time: %s"%datestring(finish))

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
        nites : A list of the nightly schedule
        """

        nites = odict()
        
        for tstart,tend in self.observation_windows:
            if start is not None and ephem.Date(tstart) < ephem.Date(start): continue
            if end is not None and ephem.Date(tend) > ephem.Date(end): continue

            chunks = self.schedule_nite(tstart,chunk,clip=True,plot=False,mode=mode)
            nite_name = '%d%02d%02d'%tstart.tuple()[:3]
            nites[nite_name] = chunks

            if plot:
                field_select = self.completed_fields[-1:]
                ortho.plotField(field_select,self.target_fields,self.completed_fields)

                #self.plotField(end,field_select)
                if (raw_input(' ...continue ([y]/n)').lower()=='n'): 
                    break

        if plot: raw_input(' ...finish... ')
        return nites

    def write(self,filename):
        self.scheduled_fields.write(filename)

    @classmethod
    def common_parser(cls):
        from maglites.utils.parser import Parser, DatetimeAction

        description = __doc__
        parser = Parser(description=description)
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
