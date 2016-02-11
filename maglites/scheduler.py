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
import maglites.utils.fileio as fileio

from maglites.utils.ortho import datestring
from maglites.field import FieldArray
############################################################

class Scheduler(object):

    def __init__(self,target_fields=None,observation_windows=None,completed_fields=None):
        self.loadTargetFields(target_fields)
        self.loadObservationWindows(observation_windows)
        self.loadCompletedFields(completed_fields)
        
        self.scheduled_fields = FieldArray()
        self.observed_fields  = self.loadObservedFields()
        self.completed_fields = self.completed_fields + self.observed_fields

        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        self.loadBlancoConstraints()

    def loadTargetFields(self, target_fields=None):
        if target_fields is None:
            target_fields = os.path.expandvars("$MAGLITESDIR/maglites/data/maglites-target-fields.csv")
        
        if isinstance(target_fields,basestring):
            self.target_fields = FieldArray.read(target_fields)
        else:
            self.target_fields = target_fields
        
    def loadObservationWindows(self, observation_windows=None):
        if observation_windows is None: 
            observation_windows = os.path.expandvars("$MAGLITESDIR/maglites/data/maglites-windows.csv")
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
            logging.info('  %s -- %s'%(start,end))
        logging.info(30*'-')

    def loadCompletedFields(self, completed_fields=None):
        if not completed_fields:
            self.completed_fields = FieldArray()
            return

        if isinstance(completed_fields,basestring):
            completed_fields = [completed_fields]

        if isinstance(completed_fields,list):
            fields = FieldArray()
            for filename in completed_fields:
                fields = fields + FieldArray.read(filename)
            self.completed_fields = fields

    def loadObservedFields(self, **kwargs):
        """
        Get observed fields from the telemetry database.
        """

        from maglites.field import FieldArray
        try: 
            fields = FieldArray.load_database()
        except: 
            fields = FieldArray()
        return fields

    def loadBlancoConstraints(self):
        """
        Load telescope pointing constraints
        """
        # Updated to remove the dependence on scipy (which is broken on the mountain)
        data = np.recfromtxt('%s/maglites/data/blanco_hour_angle_limits.dat'%(os.environ['MAGLITESDIR']), names=True)
        self.blanco_constraints = data
        ha_degrees = np.tile(0., len(self.blanco_constraints['HA']))
        for ii in range(0, len(self.blanco_constraints['HA'])):
            ha_degrees[ii] = maglites.utils.projector.hms2dec(self.blanco_constraints['HA'][ii])
        
        self.f_hour_angle_limit = lambda dec: np.interp(dec,self.blanco_constraints['Dec'], ha_degrees, left=-1, right=-1)
        self.f_airmass_limit = lambda dec: np.interp(dec,self.blanco_constraints['Dec'], self.blanco_constraints['AirmassLimit'], left=-1, right=-1)

        return self.f_hour_angle_limit,self.f_airmass_limit

    def selectField(self, date, ra_previous=None, dec_previous=None, plot=False, mode='balance'):
        """
        Input is pyephem date object
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

        # Need to figure out what to do if there are no available fields

        # Now apply some kind of selection criteria, e.g., select the field with the lowest airmass
        #airmass[np.logical_not(cut)] = 999.
        
        if mode == 'airmass':
            airmass_effective = copy.copy(airmass)
            airmass_effective[np.logical_not(cut)] = 999. # Do not observe fields that are unavailable
            airmass_effective += self.target_fields['TILING'] # Priorize coverage over multiple tilings
            index_select = np.argmin(airmass_effective)
        elif mode == 'ra':
            # Different selection
            #ra_effective = copy.copy(self.target_fields['RA'])
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = 9999.
            ra_effective += 360. * self.target_fields['TILING']
            index_select = np.argmin(ra_effective)
        elif mode == 'slew':
            #ra_effective = copy.copy(self.target_fields['RA'])
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = 9999.
            ra_effective += 360. * self.target_fields['TILING']
            ra_effective += slew**2
            #ra_effective += 2. * slew
            index_select = np.argmin(ra_effective)
        elif mode == 'balance':
            """
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = 9999.
            ra_effective += 360. * self.target_fields['TILING']
            #ra_effective += 720. * self.target_fields['TILING']
            ra_effective += slew**2
            ra_effective += 100. * (airmass - 1.)**3
            weight = ra_effective
            index_select = np.argmin(weight)
            weight = hour_angle_degree
            """
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = 9999.
            weight += 3. * 360. * self.target_fields['TILING']
            weight += slew**3 # slew**2
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'balance2':
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = 9999.
            weight += 360. * self.target_fields['TILING']
            weight += slew_ra**2
            weight += slew_dec
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'balance3':
            print slew
            weight = copy.copy(hour_angle_degree)
            weight[np.logical_not(cut)] = 9999.
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
            #plt.figure()
            #plt.plot(x_slew, y_slew)
            #raw_input('WAIT')
            weight += np.interp(slew, x_slew, y_slew, left=9999., right=9999.)
            weight += 100. * (airmass - 1.)**3
            index_select = np.argmin(weight)
        elif mode == 'airmass2':
            weight = 200. * (airmass - airmass_next)
            weight[np.logical_not(cut)] = 9999.
            weight += 360. * self.target_fields['TILING']
            weight += 100. * (airmass - 1.)**3
            weight += slew**2
            index_select = np.argmin(weight)
        
        # Search for other exposures in the same field
        field_id = self.target_fields['HEX'][index_select]
        tiling = self.target_fields['TILING'][index_select]        
        index_select = np.nonzero( (self.target_fields['HEX']==field_id) & \
                                   (self.target_fields['TILING']==tiling) & cut)[0]

        timedelta = constants.FIELDTIME*np.arange(len(index_select))
        if np.any(slew[index_select] > 5.):
            # Apply a 30 second penalty for longer slews for more accurate predictions
            timedelta += 30*ephem.second
        field_select = self.target_fields[index_select]
        field_select['AIRMASS'] = airmass[index_select]
        field_select['DATE'] = map(datestring,date+timedelta)
        field_select['SLEW'] = slew[index_select]
        field_select['MOONANGLE'] = moon_angle[index_select]
        field_select['HOURANGLE'] = hour_angle_degree[index_select]

        print np.sum(cut), len(field_select), field_select['FILTER'], 
        print np.unique(field_select['AIRMASS']), np.unique(field_select['SLEW'])

        # For diagnostic purposes
        #if len(self.accomplished_fields) % 10 == 0:
        #    self.plotWeight(date, field_select, weight)
        #    raw_input('WAIT')

        return field_select

    def plotWeight(self, date, field_select, weight):
        if plt.get_fignums(): plt.cla()
        fig, basemap = maglites.utils.ortho.makePlot(date,name='weight')
        
        index_sort = np.argsort(weight)[::-1]
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][index_sort], self.target_fields['DEC'][index_sort])
        weight_min = np.min(weight)
        basemap.scatter(*proj, c=weight[index_sort], edgecolor='none', s=50, vmin=weight_min, vmax=weight_min + 300., cmap='Spectral')
        #colorbar = plt.colorbar(label='Weight')

        #cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
        #proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        #basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        
        """
        cut_accomplished = np.in1d(self.target_fields['ID'],self.accomplished_fields['ID'])
        proj = maglites.utils.ortho.safeProj(basemap, 
                                             self.target_fields['RA'][~cut_accomplished], 
                                             self.target_fields['DEC'][~cut_accomplished])
        basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
        """

        # Draw colorbar in existing axis
        if len(fig.axes) == 2:
            colorbar = plt.colorbar(label='Weight',cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar(label='Weight')
            
        # Show the selected field
        proj = maglites.utils.ortho.safeProj(basemap, [field_select['RA']], [field_select['DEC']])
        basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

        plt.draw()
        time.sleep(0) # 0.1

    def plotField(self, date, field_select):
        if plt.get_fignums(): plt.cla()
        fig, basemap = maglites.utils.ortho.makePlot(date,name='ortho')

        """
        # Plot airmass
        cut_completed = np.in1d(self.target_fields['ID'], self.completed_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_completed], self.target_fields['DEC'][cut_completed])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo], self.target_fields['DEC'][cut_todo])
        basemap.scatter(*proj, c=airmass[cut_todo], edgecolor='none', s=50, vmin=1., vmax=2., cmap='summer_r')
        #basemap.scatter(*proj, c=cut_airmass.astype(float)[cut_todo], edgecolor='none', s=50, vmin=0., vmax=1., cmap='summer_r')
        colorbar = plt.colorbar(label='Airmass')
        """
        """
        # Plot hour angle
        cut_completed = np.in1d(self.target_fields['ID'], self.completed_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_completed], self.target_fields['DEC'][cut_completed])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo], self.target_fields['DEC'][cut_todo])
        basemap.scatter(*proj, c=np.fabs(hour_angle_degree[cut_todo]), edgecolor='none', s=50, vmin=0., vmax=78.75, cmap='summer_r')
        #basemap.scatter(*proj, c=cut_hour_angle.astype(float)[cut_todo], edgecolor='none', s=50, vmin=0., vmax=1., cmap='summer_r')
        colorbar = plt.colorbar(label='Hour Angle')
        """
        """
        # Plot RA
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo], self.target_fields['DEC'][cut_todo])
        ra_effective = self.target_fields['RA'][cut_todo] - ra_zenith
        ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
        basemap.scatter(*proj, c=ra_effective, edgecolor='none', s=50, cmap='summer_r')
        colorbar = plt.colorbar(label='RA')

        cut_completed = np.in1d(self.target_fields['ID'], self.completed_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_completed], self.target_fields['DEC'][cut_completed])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        """
        """
        # Plot weight
        index_sort = np.argsort(ra_effective[cut_todo])[::-1]
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo][index_sort], self.target_fields['DEC'][cut_todo][index_sort])
        weight_min = np.min(ra_effective[cut_todo])
        basemap.scatter(*proj, c=ra_effective[cut_todo][index_sort], edgecolor='none', s=50, vmin=weight_min, vmax=weight_min + 100., cmap='summer_r')
        colorbar = plt.colorbar(label='Weight')

        cut_completed = np.in1d(self.target_fields['ID'], self.completed_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_completed], self.target_fields['DEC'][cut_completed])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        """

        # ADW: Need to be careful about the size of the marker. It
        # does not change with the size of the frame so it is
        # really safest to scale to the size of the zenith circle
        # (see PlotPointings). That said, s=50 is probably roughly ok.
        
        # Plot number of tilings 
        cut_completed = np.in1d(self.target_fields['ID'],
                                self.completed_fields['ID'])
        proj = maglites.utils.ortho.safeProj(basemap, 
                                             self.target_fields['RA'][~cut_completed], 
                                             self.target_fields['DEC'][~cut_completed])
        basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_completed))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_completed], self.target_fields['DEC'][cut_completed])
        basemap.scatter(*proj, c=self.target_fields['TILING'][cut_completed], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')

        # Draw colorbar in existing axis
        if len(fig.axes) == 2:
            colorbar = plt.colorbar(label='Tiling',cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar(label='Tiling')
            
        # Show the selected field
        proj = maglites.utils.ortho.safeProj(basemap, [field_select['RA']], [field_select['DEC']])
        basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

        plt.draw()
        time.sleep(0.1)


    def run(self, tstart=None, tstop=None, plot=True):
        # Reset the scheduled fields
        self.scheduled_fields = FieldArray(0)

        # If no tstop, run for 90 minutes
        timedelta = 90*ephem.minute
        if tstart is None: tstart = ephem.now()
        if tstop is None: tstop = tstart + timedelta
        print 'start:',datestring(tstart), 'stop:',datestring(tstop)

        # Convert strings into dates
        if isinstance(tstart,basestring):
            tstart = ephem.Date(tstart)
        if isinstance(tstop,basestring):
            tstop = ephem.Date(tstop)

        msg = "Previously completed fields: %i"%len(self.completed_fields)
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
                    msg = 'Date outside of nominal observing windows'
                    logging.warning(msg)
            logging.info('  '+datestring(date))

            # Check 
            compute_slew = True
            if len(self.completed_fields) == 0:
                compute_slew = False
            else:
                if (date - ephem.Date(self.completed_fields['DATE'][-1])) > (30. * ephem.minute):
                    compute_slew = False
            if compute_slew:
                field_select = self.selectField(date, ra_previous=self.completed_fields['RA'][-1], dec_previous=self.completed_fields['DEC'][-1], plot=plot)
            else:
                field_select = self.selectField(date, plot=plot)

            id_select = field_select['ID']
            # Previously, increment time by a constant
            #date = date + len(field_select)*constants.FIELDTIME
            # Now update the time from the selected field
            date = ephem.Date(field_select[-1]['DATE']) + constants.FIELDTIME

            self.completed_fields = self.completed_fields + field_select
            self.scheduled_fields = self.scheduled_fields + field_select

            if plot: self.plotField(date, field_select)

            if date > tstop: break
        print "Newly scheduled fields: %i"%len(self.scheduled_fields)
        return self.scheduled_fields

    schedule_chunk = run

    def schedule_nite(self,nite=None,chunk=60.,plot=False):
        # Create the nite
        nite = ephem.Date(nite) if nite else ephem.now()
        nite_tuple = nite.tuple()[:3]

        # Convert chunk to MJD
        if chunk > 1: chunk = chunk*ephem.minute

        try:
            nites = [w[0] for w in self.observation_windows]
            nite_tuples = [n.tuple()[:3] for n in nites]
            idx = nite_tuples.index(nite_tuple)
            start,finish = self.observation_windows[idx]
        except (TypeError, ValueError):
            msg = "Requested nite not found in windows:\n"
            msg += "%s/%s/%s : "%nite_tuple
            msg += '['+', '.join(['%s/%s/%s'%t for t in nite_tuples])+']'
            logging.warning(msg)
            
            sun = ephem.Sun()
            obs = self.observatory
            if obs.previous_setting(sun) > obs.previous_rising(sun):
                # We are at night
                start = obs.previous_setting(sun) + 1*ephem.hour
            else:
                start = obs.next_setting(sun) + 1*ephem.hour
            finish = obs.next_rising(sun) - 1*ephem.hour

            logging.debug("Night start time: %s"%datestring(start))
            logging.debug("Night finish time: %s"%datestring(finish))

        chunks = []

        i = 0
        while start < finish:
            i+=1
            msg = "Scheduling %s -- Chunk %i"%(start,i)
            logging.debug(msg)
            end = start+chunk
            scheduled_fields = self.run(start, end, plot=False)
            if plot:
                field_select = scheduled_fields[-1]
                self.plotField(end,field_select)
                if (raw_input(' ...continue ([y]/n)').lower()=='n'): 
                    break
            chunks.append(scheduled_fields)
            start = end
        return chunks

    def schedule_survey(self, chunk=60., plot=False):
        nites = odict()

        for start,end in self.observation_windows:
            chunks = self.schedule_nite(start,chunk,plot=False)
            nite_name = '%d%02d%02d'%start.tuple()[:3]
            nites[nite_name] = chunks

            field_select = self.completed_fields[-1]
            self.plotField(end,field_select)
            if (raw_input(' ...continue ([y]/n)').lower()=='n'): 
                break
            
        return nites

    def write(self,filename):
        self.scheduled_fields.write(filename)

    @classmethod
    def common_parser(cls):
        from maglites.utils.parser import Parser

        description = __doc__
        parser = Parser(description=description)
        parser.add_argument('-p','--plot',action='store_true',
                            help='create visual output.')
        parser.add_argument('--utc-start',
                            help="start time for observation.")
        parser.add_argument('--utc-end',
                            help="end time for observation.")
        parser.add_argument('-k','--chunk', default=60., type=float,
                            help = 'time chunk')
        parser.add_argument('-f','--fields',default=None,
                            help='all target fields.')
        parser.add_argument('-w','--windows',default=None,
                            help='observation windows.')
        parser.add_argument('-c','--complete',nargs='?',action='append',
                            help="fields that have been completed.")
        parser.add_argument('-o','--outfile',default=None,
                            help='save output file of scheduled fields.')

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
