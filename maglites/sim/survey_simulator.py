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

import os
import copy
import numpy as np
import scipy.interpolate
import ephem

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.ortho


class FormatFloatForce(mlab.FormatFormatStr): 
    """
    mlab not doing well...
    """
    def __init__(self,fmt="%.4f"): 
        mlab.FormatFormatStr.__init__(self,fmt) 
    def toval(self, x): 
        return x 
    def fromstr(self, s): 
        return float(s) 

############################################################

class Simulator(object):

    def __init__(self, infile_target_fields, infile_observation_windows=None, infile_accomplished_fields=None):
        self.target_fields = np.recfromtxt(infile_target_fields, names=True)
        if not infile_accomplished_fields:
            self.accomplished_field_ids = []
            self.accomplished_fields = {}
            for key in ['ID', 'RA', 'DEC', 'TILING', 'PRIORITY', 'DATE', 'AIRMASS', 'SLEW', 'MOONANGLE', 'HOURANGLE']:
                self.accomplished_fields[key] = []
        else:
            self.accomplished_fields = np.recfromtxt(infile_accomplished_fields, delimiter=',', names=True)
            
        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        if infile_observation_windows is not None:
            self.loadObservationWindows(infile_observation_windows)
        else:
            self.observation_windows = None

        self.loadBlancoConstraints()

    def loadObservationWindows(self, infile_observation_windows):
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

    def loadBlancoConstraints(self):
        """
        Load telescope pointing constraints
        """
        data = np.recfromtxt('%s/maglites/data/blanco_hour_angle_limits.dat'%(os.environ['MAGLITESDIR']), names=True)
        self.blanco_constraints = data
        ha_degrees = np.tile(0., len(self.blanco_constraints['HA']))
        for ii in range(0, len(self.blanco_constraints['HA'])):
            ha_degrees[ii] = maglites.utils.projector.hms2dec(self.blanco_constraints['HA'][ii])

        self.f_hour_angle_limit = scipy.interpolate.interp1d(self.blanco_constraints['Dec'], ha_degrees, 
                                                             bounds_error=False, fill_value=-1.)
        self.f_airmass_limit = scipy.interpolate.interp1d(self.blanco_constraints['Dec'], self.blanco_constraints['AirmassLimit'], 
                                                          bounds_error=False, fill_value=-1.)

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
            ra_effective += slew**2
            ra_effective += 100. * (airmass - 1.)**3
            index_select = np.argmin(ra_effective)

        #plt.figure()
        #plt.scatter(np.arange(len(ra_effective)), ra_effective)
        #raw_input()

        print np.sum(cut), self.target_fields['TILING'][index_select], airmass[index_select], slew[index_select]

        if plot:
            #fig, ax, basemap = maglites.utils.ortho.makePlot(date)
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
            cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
            cut_accomplished[index_select] = True
            proj = maglites.utils.ortho.safeProj(basemap, 
                                                 self.target_fields['RA'][np.logical_not(cut_accomplished)], 
                                                 self.target_fields['DEC'][np.logical_not(cut_accomplished)])
            basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
            
            proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
            basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
            colorbar = plt.colorbar(label='Tiling')
            

            # Show the selected field
            proj = maglites.utils.ortho.safeProj(basemap, [self.target_fields['RA'][index_select]], [self.target_fields['DEC'][index_select]])
            basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)


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
                    raw_input(' ...continue')
                plt.clf()

            date = date + (4. * ephem.minute) # 4 minutes = 90 sec exposures in g and r with 30 sec between exposures
            self.accomplished_field_ids.append(id_select)
            plot = False # Kludge to make end of night summary plots

            for key in field_select.keys():
                self.accomplished_fields[key].append(field_select[key])
                #if key not in self.accomplished_fields.keys():
                #    self.accomplished_fields[key] = [field_select[key]]
                #else:
                #    self.accomplished_fields[key].append(field_select[key])

        print len(self.accomplished_field_ids)

        # Clean up
        self.accomplished_field_ids = []

    def consolidateAccomplishedFields(self):
        """
        slew = np.tile(0., len(self.accomplished_fields['RA']))
        for ii in range(1, len(self.accomplished_fields['RA'])):
            time_difference = ephem.Date(self.accomplished_fields['DATE'][ii]) - ephem.Date(self.accomplished_fields['DATE'][ii - 1]) 
            angsep = maglites.utils.projector.angsep(self.accomplished_fields['RA'][ii], self.accomplished_fields['DEC'][ii], 
                                                     self.accomplished_fields['RA'][ii - 1], self.accomplished_fields['DEC'][ii - 1])
            if time_difference > (30. * ephem.minute):
                angsep = 0.
            print time_difference, angsep
            slew[ii] = angsep
        """

        data = zip(self.accomplished_fields['ID'],
                   self.accomplished_fields['RA'],
                   self.accomplished_fields['DEC'],
                   self.accomplished_fields['TILING'],
                   self.accomplished_fields['PRIORITY'],
                   self.accomplished_fields['DATE'],
                   self.accomplished_fields['AIRMASS'],
                   self.accomplished_fields['SLEW'],
                   self.accomplished_fields['MOONANGLE'],
                   self.accomplished_fields['HOURANGLE'])
        #slew]
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
        accomplished_fields = np.array(data, dtype=dtype)

        return accomplished_fields

    def saveAccomplishedFields(self, outfile):
        ### ADW: It would probably be better to use pylab.rec2csv or
        ### some other standard csv creation routine. The spaces after
        ### the commas make me worried...
        data = self.consolidateAccomplishedFields()
        #np.savetxt(outfile, data, fmt='%i, %.4f, %.4f, %i, %i, %s, %.4f, %.4f, %.4f, %.4f', header='ID, RA, DEC, TILING, PRIORITY, DATE, AIRMASS, SLEW, MOONANGLE, HOURANGLE')
        formatd = dict()
        for name,(dtype,size) in data.dtype.fields.items():
            if dtype.kind == 'f': formatd[name] = FormatFloatForce()

        print formatd
        mlab.rec2csv(data,outfile,formatd=formatd)

############################################################

def main():
    my_simulator = Simulator('target_fields.txt')
    my_simulator.loadObservationWindows('observation_windows.txt')
    my_simulator.run(plot=False)
    my_simulator.saveAccomplishedFields('accomplished_fields_2.txt')

    # Example diagnostics
    """
    d = np.recfromtxt('accomplished_fields_2.txt', delimiter=',', names=True)
    plt.figure()
    plt.scatter(np.arange(len(d)), d['SLEW'])

    plt.figure()
    #plt.hist(d['AIRMASS'], bins=31)
    plt.scatter(np.arange(len(d)), d['AIRMASS'])
    """

if __name__ == '__main__':
    main()

############################################################
