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

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.constants as constants
import maglites.utils.ortho
import maglites.utils.fileio as fileio

from maglites.utils.ortho import datestring
from maglites.field import FieldArray
from maglites.utils import Database
############################################################

class Scheduler(object):

    def __init__(self, target_fields, observation_windows=None, accomplished_fields=None):
        if isinstance(target_fields,basestring):
            self.target_fields = FieldArray.read(target_fields)
        else:
            self.target_fields = target_fields

        self.loadObservationWindows(observation_windows)
        self.loadAccomplishedFields(accomplished_fields)
        
        self.scheduled_fields = FieldArray()
        self.observed_fields  = self.loadObservedFields()
        self.accomplished_fields = self.accomplished_fields + self.observed_fields

        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        self.loadBlancoConstraints()

    def loadObservationWindows(self, observation_windows = None):
        if not observation_windows: 
            self.observation_windows = None            
            return
        
        if isinstance(observation_windows,basestring):
            observation_windows = np.recfromcsv(observation_windows)
            
        self.observation_windows = []
        for start,end in observation_windows:
            self.observation_windows.append([ephem.Date(start), ephem.Date(end)])

        # Do a sanity check to make sure that observation windows are properly sorted
        for ii in range(0, len(self.observation_windows)):
            if self.observation_windows[ii][1] < self.observation_windows[ii][0]:
                logging.warning('Observation windows are not properly sorted')
            if ii > 0:
                if self.observation_windows[ii][0] < self.observation_windows[ii - 1][1]:
                    logging.warning('Observation windows are not properly sorted')

        logging.info('Observation Windows:')
        for start,end in self.observation_windows:
            logging.info('  %s -- %s'%(start,end))
        logging.info(30*'-')

    def loadAccomplishedFields(self, accomplished_fields = None):
        if not accomplished_fields:
            self.accomplished_fields = FieldArray()
        if isinstance(accomplished_fields,basestring):
            accomplished_fields = [accomplished_fields]

        if isinstance(accomplished_fields,list):
            fields = FieldArray()
            for filename in accomplished_fields:
                fields = fields + FieldArray.read(filename)
            self.accomplished_fields = fields

    def loadObservedFields(self, **kwargs):
        """
        Get the fields that have been observed from the telemetry DB.
        """
        defaults = dict(propid='2016A-0366', limit='', dbname='db-fnal')
        params = copy.deepcopy(defaults)
        params.update(kwargs)

        db = Database()
        db.connect()

        query ="""
        select object as OBJECT, seqid as SEQID, seqnum as SEQNUM,
        telra as RA, teldec as DEC, 
        to_char(utc_beg, 'YYYY/MM/DD HH24:MI:SS') AS DATE, 
        COALESCE(airmass,-1) as AIRMASS, COALESCE(moonangl,-1) as MOONANGLE, 
        COALESCE(ha,-1) as HOURANGLE, COALESCE(slewangl,-1) as SLEW 
        from exposure where propid = '%(propid)s' %(limit)s
        """%params

        data = db.execute(query)
        names = map(str.upper,db.get_columns())
        
        if not len(data): return FieldArray(0)

        fields = FieldArray.load_recarray(recarray)
        fields.from_object(recarray['OBJECT'])
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

        # Include moon angle
        moon = ephem.Moon()
        moon.compute(date)
        ra_moon = np.degrees(moon.ra)
        dec_moon = np.degrees(moon.dec)
        moon_angle = maglites.utils.projector.angsep(ra_moon, dec_moon, self.target_fields['RA'], self.target_fields['DEC'])

        # Slew from the previous pointing
        if ra_previous is not None and dec_previous is not None:
            slew = maglites.utils.projector.angsep(ra_previous, dec_previous, self.target_fields['RA'], self.target_fields['DEC'])
        else:
            slew = np.tile(0., len(self.target_fields['RA']))

        # Hour angle restrictions
        hour_angle_degree = copy.copy(self.target_fields['RA']) - ra_zenith
        hour_angle_degree[hour_angle_degree > 180.] = hour_angle_degree[hour_angle_degree > 180.] - 360.
        cut_hour_angle = np.fabs(hour_angle_degree) < self.f_hour_angle_limit(self.target_fields['DEC']) # Check the hour angle restrictions at south pole
        
        # Airmass restrictions
        cut_airmass = airmass < self.f_airmass_limit(self.target_fields['DEC'])

        # Declination restrictions
        cut_declination = self.target_fields['DEC'] > -89.

        # Don't consider fields which have already been observed
        cut_todo = np.logical_not(np.in1d(self.target_fields['ID'], self.accomplished_fields['ID']))
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
            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
            ra_effective[np.logical_not(cut)] = 9999.
            ra_effective += 360. * self.target_fields['TILING']
            #ra_effective += 720. * self.target_fields['TILING']
            ra_effective += slew**2
            ra_effective += 100. * (airmass - 1.)**3
            index_select = np.argmin(ra_effective)

        # Search for other exposures in the same field
        field_id = self.target_fields['SMASH_ID'][index_select]
        tiling = self.target_fields['TILING'][index_select]        
        index_select = np.nonzero( (self.target_fields['SMASH_ID']==field_id) & \
                                   (self.target_fields['TILING']==tiling) & cut)[0]



        field_select = self.target_fields[index_select]
        field_select['AIRMASS'] = airmass[index_select]
        field_select['DATE'] = maglites.utils.ortho.datestring(date)
        field_select['SLEW'] = slew[index_select]
        field_select['MOONANGLE'] = moon_angle[index_select]
        field_select['HOURANGLE'] = hour_angle_degree[index_select]

        print np.sum(cut), len(field_select), field_select['FILTER'], 
        print np.unique(field_select['AIRMASS']), np.unique(field_select['SLEW'])

        return field_select


    def plotField(self, date, field_select):
        if plt.get_fignums(): plt.cla()
        fig, basemap = maglites.utils.ortho.makePlot(date,name='ortho')

        """
        # Plot airmass
        cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo], self.target_fields['DEC'][cut_todo])
        basemap.scatter(*proj, c=airmass[cut_todo], edgecolor='none', s=50, vmin=1., vmax=2., cmap='summer_r')
        #basemap.scatter(*proj, c=cut_airmass.astype(float)[cut_todo], edgecolor='none', s=50, vmin=0., vmax=1., cmap='summer_r')
        colorbar = plt.colorbar(label='Airmass')
        """
        """
        # Plot hour angle
        cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
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

        cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        """
        """
        # Plot weight
        index_sort = np.argsort(ra_effective[cut_todo])[::-1]
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo][index_sort], self.target_fields['DEC'][cut_todo][index_sort])
        weight_min = np.min(ra_effective[cut_todo])
        basemap.scatter(*proj, c=ra_effective[cut_todo][index_sort], edgecolor='none', s=50, vmin=weight_min, vmax=weight_min + 100., cmap='summer_r')
        colorbar = plt.colorbar(label='Weight')

        cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
        """

        # ADW: Need to be careful about the size of the marker. It
        # does not change with the size of the frame so it is
        # really safest to scale to the size of the zenith circle
        # (see PlotPointings). That said, s=50 is probably roughly ok.
        
        # Plot number of tilings 
        cut_accomplished = np.in1d(self.target_fields['ID'],self.accomplished_fields['ID'])
        proj = maglites.utils.ortho.safeProj(basemap, 
                                             self.target_fields['RA'][~cut_accomplished], 
                                             self.target_fields['DEC'][~cut_accomplished])
        basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
        
        proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
        basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')

        # Draw colorbar in existing axis
        if len(fig.axes) == 2:
            colorbar = plt.colorbar(label='Tiling',cax=fig.axes[-1])
        else:
            colorbar = plt.colorbar(label='Tiling')
            
        # Show the selected field
        proj = maglites.utils.ortho.safeProj(basemap, [field_select['RA']], [field_select['DEC']])
        basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

        plt.draw()
        time.sleep(0.2)


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
                    msg = 'Date outside of nominal observing windows'
                    logging.warning(msg)
            logging.info('  '+datestring(date))

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
            date = date + len(field_select)*constants.FIELDTIME

            self.accomplished_fields = self.accomplished_fields + field_select
            self.scheduled_fields    = self.scheduled_fields + field_select

            if plot: self.plotField(date, field_select)

            if date > tstop: break
        print "Newly scheduled fields: %i"%len(self.scheduled_fields)
        return self.scheduled_fields

    def write(self,filename):
        self.scheduled_fields.write(filename)

    @classmethod
    def common_parser(cls):
        import logging
        import argparse

        class VerboseAction(argparse._StoreTrueAction):

            def __call__(self, parser, namespace, values, option_string=None):
                super(VerboseAction,self).__call__(parser, namespace, values, option_string)
                #setattr(namespace, self.dest, self.const)
                if self.const: logging.getLogger().setLevel(logging.DEBUG)

        class SpecialFormatter(logging.Formatter):
            """
            Class for overloading log formatting based on level.
            """
            FORMATS = {'DEFAULT'       : "%(message)s",
                       logging.WARNING : "WARNING: %(message)s",
                       logging.ERROR   : "ERROR: %(message)s",
                       logging.DEBUG   : "DEBUG: %(message)s"}
         
            def format(self, record):
                self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
                return logging.Formatter.format(self, record)
         
        logger = logging.getLogger()
        handler = logging.StreamHandler()
        handler.setFormatter(SpecialFormatter())
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

        description = __doc__
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument('-v','--verbose',action=VerboseAction,
                            help='Output verbosity.')
        parser.add_argument('-p','--plot',action='store_true',
                            help='Plot output.')
        parser.add_argument('--utc-start',
                            help="Start time for observation.")
        parser.add_argument('--utc-end',
                            help="End time for observation.")
        parser.add_argument('-f','--fields',default='target_fields.csv',
                            help='List of all target fields.')
        parser.add_argument('-w','--windows',default='observation_windows.csv',
                            help='List of observation windows.')
        parser.add_argument('-a','--accomplish',action='append',
                            help="List of fields that have been accomplished.")
        parser.add_argument('-o','--outfile',default='scheduled_fields.csv',
                            help='Save output file of scheduled fields.')

        return parser

    @classmethod
    def parser(cls):
        return cls.common_parser()

    @classmethod
    def main(cls):
        args = cls.parser().parse_args()
        scheduler = cls(args.fields,args.windows,args.accomplish)
        scheduler.run(plot=args.plot)
        if args.outfile: 
            scheduler.scheduled_fields.write(args.outfile)
         
        return scheduler


############################################################



if __name__ == '__main__':
    scheduler = Scheduler.main()
############################################################
