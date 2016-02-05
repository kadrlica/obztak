# Make a Simulator class
# Keeps track of time, what exposures have been taken
# What exposures need to be taken
# Selects exposure
# Visualize

# Inputs: 
# * Listing of all desired pointings (RA, Dec, tiling, priority)
# * Listing of all the exposures taken (RA, Dec, tiling, airmass, time)
# * Listing of 

# Functions:
# * Run full set of observations
# * Select next sequence of exposures from a given start time
# * Visualize (afterburner?)

import os,sys
import copy
import numpy as np
import time
import ephem
import matplotlib.pyplot as plt

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.constants as constants
import maglites.utils.ortho
import maglites.utils.fileio as fileio

############################################################

class Simulator(object):

    def __init__(self, infile_target_fields, infile_observation_windows=None, infile_accomplished_fields=None):
        self.target_fields = np.recfromtxt(infile_target_fields, names=True)

        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        self.loadObservationWindows(infile_observation_windows)
        self.loadAccomplishedFields(infile_accomplished_fields)

        self.loadBlancoConstraints()

    def loadObservationWindows(self, infile_observation_windows = None):
        if not infile_observation_windows: 
            self.observation_windows = None            
            return

        reader = open(infile_observation_windows)
        lines = reader.readlines()
        reader.close()

        self.observation_windows = []
        for line in lines:
            parts = line.strip().split(',')
            if len(parts) != 2:
                print 'WARNING! Unable to parse line: %s'%(line)
            self.observation_windows.append([ephem.Date(parts[0]), ephem.Date(parts[1])])

        # Do a sanity check to make sure that observation windows are properly sorted
        for ii in range(0, len(self.observation_windows)):
            if self.observation_windows[ii][1] < self.observation_windows[ii][0]:
                print 'WARNING! Observation windows are not properly sorted'
            if ii > 0:
                if self.observation_windows[ii][0] < self.observation_windows[ii - 1][1]:
                    print 'WARNING! Observation windows are not properly sorted'

        print 'Observation Windows'
        for observation_window in self.observation_windows:
            print '  %s -- %s'%(observation_window[0], observation_window[1])

    def loadAccomplishedFields(self, infile_accomplished_fields = None):
        if not infile_accomplished_fields:
            self.accomplished_fields = self.createFieldArray()
        else:
            self.accomplished_fields = self.loadFields(infile_accomplished_fields)
        self.accomplished_field_ids = self.accomplished_fields['ID'].tolist()


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

    def loadBlancoConstraints2(self):
        """
        Load telescope pointing constraints
        """
        import scipy.interpolate
        data = np.recfromtxt('%s/maglites/data/blanco_hour_angle_limits.dat'%(os.environ['MAGLITESDIR']), names=True)
        self.blanco_constraints = data
        ha_degrees = np.tile(0., len(self.blanco_constraints['HA']))
        for ii in range(0, len(self.blanco_constraints['HA'])):
            ha_degrees[ii] = maglites.utils.projector.hms2dec(self.blanco_constraints['HA'][ii])

        self.f_hour_angle_limit = scipy.interpolate.interp1d(self.blanco_constraints['Dec'], ha_degrees, 
                                                             bounds_error=False, fill_value=-1.)
        
        self.f_airmass_limit = scipy.interpolate.interp1d(self.blanco_constraints['Dec'], self.blanco_constraints['AirmassLimit'], 
                                                          bounds_error=False, fill_value=-1.)

        return self.f_hour_angle_limit,self.f_airmass_limit


    def selectField(self, date, ra_previous=None, dec_previous=None, plot=False, mode='balance'):
        """
        Input is pyephem date object
        """

        if type(date) == ephem.Date:
            self.observatory.date = date
        else:
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
        cut_todo = np.logical_not(np.in1d(self.target_fields['ID'], self.accomplished_field_ids))
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

        #plt.figure()
        #plt.scatter(np.arange(len(ra_effective)), ra_effective)
        #raw_input()

        print np.sum(cut), self.target_fields['TILING'][index_select], airmass[index_select], slew[index_select]

        if plot:
            #fig, ax, basemap = maglites.utils.ortho.makePlot(date)

            if plt.get_fignums(): plt.cla()
            fig, basemap = maglites.utils.ortho.makePlot(date)

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
            
            # Plot number of tilings 
            # ADW: Need to be careful about the size of the marker. It
            # does not change with the size of the frame so it is
            # really safest to scale to the size of the zenith circle
            # (see PlotPointings). That said, s=50 is probably roughly ok.
            cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
            cut_accomplished[index_select] = True
            proj = maglites.utils.ortho.safeProj(basemap, 
                                                 self.target_fields['RA'][np.logical_not(cut_accomplished)], 
                                                 self.target_fields['DEC'][np.logical_not(cut_accomplished)])
            basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
            
            proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
            basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')

            # Draw colorbar in existing axis
            if len(fig.axes) == 2:
                colorbar = plt.colorbar(label='Tiling',cax=fig.axes[-1])
            else:
                colorbar = plt.colorbar(label='Tiling')
                
            # Show the selected field
            proj = maglites.utils.ortho.safeProj(basemap, [self.target_fields['RA'][index_select]], [self.target_fields['DEC'][index_select]])
            basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

            plt.draw()
            time.sleep(0.2)

        field_select_dict = {}
        for name in self.target_fields.dtype.names:
            field_select_dict[name] = self.target_fields[name][index_select]
        field_select_dict['AIRMASS'] = airmass[index_select]
        field_select_dict['DATE'] = maglites.utils.ortho.datestring(date)
        field_select_dict['SLEW'] = slew[index_select]
        field_select_dict['MOONANGLE'] = moon_angle[index_select]
        field_select_dict['HOURANGLE'] = hour_angle_degree[index_select]

        #return self.target_fields['ID'][index_select]
        return field_select_dict

    def run(self, date=None, plot=True):

        if date is None:
            if self.observation_windows is not None:
                date = self.observation_windows[0][0]
            else:
                print 'WARNING! Please supply a date to begin observations'
        else:
            if type(date) != ephem.Date:
                date = ephem.Date(date)
        
        latch = True
        while latch:
            # Check to see if in valid observation window

            ### ADW: It would be good to pad these windows by a bit in
            ### case we start early or end late
            if self.observation_windows is not None:
                if date < self.observation_windows[0][0]:
                    date = self.observation_windows[0][0]
                    print 'Advance to start of first night'
                if date > self.observation_windows[-1][1]:
                    print 'Reached end of observation windows'
                    latch = False
                    plot = True # Kludge to make end of night summary plots
                else:
                    for ii in range(1, len(self.observation_windows)):
                        if date > self.observation_windows[ii - 1][1] and date < self.observation_windows[ii][0]:
                            date = self.observation_windows[ii][0]
                            print 'Advance to start of next night'
                            plot = True # Kludge to make end of night summary plots

            if plot:
                print '  %s'%(maglites.utils.ortho.datestring(date)),
            else:
                print '  %s'%(maglites.utils.ortho.datestring(date))
            
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
            if plot:
                if latch:
                    latch = (raw_input(' ...continue ([y]/n)').lower() != 'n')
                else:
                    raw_input(' ...finish...')
                plt.clf()

            date = date + (4. * ephem.minute) # 4 minutes = 90 sec exposures in g and r with 30 sec between exposures
            self.accomplished_field_ids.append(id_select)
            plot = False # Kludge to make end of night summary plots

            new_field = self.createFieldArray(1)
            for key in field_select.keys():
                new_field[-1][key] = field_select[key]

        print len(self.accomplished_field_ids)

        # Clean up
        self.accomplished_field_ids = []

    @classmethod
    def createFieldArray(cls,size=0):
        dtype = [('ID', int),
                 ('RA', float),
                 ('DEC', float),
                 ('TILING', int),
                 ('PRIORITY', int),
                 ('DATE', 'a20'),
                 ('AIRMASS', float),
                 ('SLEW', float),
                 ('MOONANGLE', float),
                 ('HOURANGLE', float)]
        return np.recarray(size, dtype=dtype)
    
    @classmethod
    def saveFields(cls, filename, fields):
        base,ext = os.path.splitext(filename)
        if ext in ('.csv','.txt','.dat'):
            outdata = fields
            fileio.rec2csv(filename,outdata)
        elif ext in ('.json'):
            outdata = cls.fields2sispi(fields)
            fileio.write_json(filename,outdata)
        else:
            msg = 'Unrecognized file type'
            raise Exception(msg)

    @classmethod
    def loadFields(cls, filename):
        base,ext = os.path.splitext(filename)

        orig = cls.createFieldArray()
        dtype = copy.deepcopy(orig.dtype)

        if ext in ('.csv','.txt','.dat'):
            load = np.genfromtxt(filename,delimiter=',',names=True,dtype=dtype)
        elif ext in ('.json'):
            load = cls.sispi2fields(fileio.read_json(filename))

        orig_names = np.char.array(orig.dtype.names)
        load_names = np.char.array(load.dtype.names)
        if np.any(load_names != orig_names):
            msg =  "Unexpected input name:\n"
            msg += str(load_names[load_names != orig_names])
            raise Exception(msg)

        return load

    @classmethod
    def fields2sispi(cls,fields):
        out_dicts = []
        for field in fields:
            object_name = constants.OBJECT_FMT%(field)
            seqid = constants.SEQID_FMT%(field)
            for i,band in enumerate(constants.BANDS):
                sispi_dict = copy.deepcopy(constants.SISPI_DICT)
                sispi_dict['seqnum'] = i+1
                sispi_dict['seqid'] = seqid
                sispi_dict['object'] = object_name
                sispi_dict['RA'] = field['RA']
                sispi_dict['dec'] = field['DEC']
                sispi_dict['filter'] = band
                out_dicts.append(sispi_dict)
        return out_dicts

    @classmethod
    def sispi2fields(cls,sispi):
        fields = cls.createFieldArray()
        for row in sispi:
            if row['seqnum'] > 1: continue

            row.update(constants.OBJECT2FIELD(row['object']))
            row.update(constants.SEQID2FIELD(row['seqid']))
            
            row = dict([(str(k.upper()),v) for k,v in row.items()])
            row['AIRMASS']   = None
            row['SLEW']      = None
            row['MOONANGLE'] = None
            row['HOURANGLE'] = None

            field = cls.createFieldArray(1)
            for key in field.dtype.names:
                print key, row[key]
                field[key] = row[key]
            fields = np.append(fields,field)

        return fields


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
        parser.add_argument('-f','--fields',default='target_fields.txt',
                            help='List of all target fields.')
        parser.add_argument('-w','--windows',default='observation_windows.txt',
                            help='List of observation windows.')
        parser.add_argument('-d','--done',
                            help="List of fields that have been observed.")
        parser.add_argument('-o','--outfile',default='accomplished_fields.txt',
                            help='Save output fields surveyed.')

        return parser

    @classmethod
    def parser(cls):
        return cls.common_parser()

############################################################


def main():
    args = Simulator.parser().parse_args()
    sim = Simulator(args.fields,args.windows,args.done)
    sim.run(plot=args.plot)
    if args.outfile: 
        sim.saveFields(args.outfile,sim.accomplished_fields)
    
    return sim

if __name__ == '__main__':
    main()

############################################################
