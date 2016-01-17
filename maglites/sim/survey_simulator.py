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

import numpy as np
import ephem
import matplotlib.pyplot as plt

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.ortho

############################################################

class Simulator:

    def __init__(self, infile_target_fields, infile_observation_windows=None, infile_acomplished=None):
        self.target_fields = np.recfromtxt(infile_target_fields, names=True)
        if not infile_acomplished:
            self.accomplished_field_ids = [] 
        else:
            # Load the IDs of accomplished fields
            pass
            
        self.observatory = ephem.Observer()
        self.observatory.lon = maglites.utils.constants.LON_CTIO
        self.observatory.lat = maglites.utils.constants.LAT_CTIO
        self.observatory.elevation = maglites.utils.constants.ELEVATION_CTIO

        if infile_observation_windows is not None:
            self.loadObservationWindows(infile_observation_windows)
        else:
            self.observation_windows = None

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

    def selectField(self, date, plot=False):
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

        airmass = 1. / np.cos(np.radians(maglites.utils.projector.angsep(self.target_fields['RA'],
                                                                         self.target_fields['DEC'],
                                                                         ra_zenith,
                                                                         dec_zenith)))

        # Don't consider fields which have already been observed
        cut_todo = np.logical_not(np.in1d(self.target_fields['ID'], self.accomplished_field_ids))
        cut = cut_todo & (airmass < 2.)

        print np.sum(cut)

        # Now apply some kind of selection criteria, e.g., select the field with the lowest airmass
        airmass[np.logical_not(cut)] = 999.
        airmass_effective = airmass + self.target_fields['TILING']
        index_select = np.argmin(airmass_effective)

        if plot:
            #fig, ax, basemap = maglites.utils.ortho.makePlot(date)
            fig, basemap = maglites.utils.ortho.makePlot(date)
            
            """
            cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
            proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
            basemap.scatter(*proj, c='0.75', edgecolor='none', s=50)
            
            proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_todo], self.target_fields['DEC'][cut_todo])
            basemap.scatter(*proj, c=airmass[cut_todo], edgecolor='none', s=50, vmin=1., vmax=2., cmap='summer_r')
            colorbar = plt.colorbar(label='Airmass')
            """
            
            cut_accomplished = np.in1d(self.target_fields['ID'], self.accomplished_field_ids)
            cut_accomplished[index_select] = True

            proj = maglites.utils.ortho.safeProj(basemap, 
                                                 self.target_fields['RA'][np.logical_not(cut_accomplished)], 
                                                 self.target_fields['DEC'][np.logical_not(cut_accomplished)])
            basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=3, cmap='summer_r')
            
            proj = maglites.utils.ortho.safeProj(basemap, self.target_fields['RA'][cut_accomplished], self.target_fields['DEC'][cut_accomplished])
            basemap.scatter(*proj, c=self.target_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=3, cmap='summer_r')
            colorbar = plt.colorbar(label='Tiling')

            proj = maglites.utils.ortho.safeProj(basemap, [self.target_fields['RA'][index_select]], [self.target_fields['DEC'][index_select]])
            basemap.scatter(*proj, c='magenta', edgecolor='none', s=50)

        return self.target_fields['ID'][index_select]

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
            id_select = self.selectField(date, plot=plot)
            if plot:
                if latch:
                    latch = (raw_input(' ...continue ([y]/n)').lower() != 'n')
                else:
                    raw_input(' ...continue')
                plt.clf()

            date = date + (4. * ephem.minute) # 4 minutes = 90 sec exposures in g and r with 30 sec between exposures
            self.accomplished_field_ids.append(id_select)
            plot = False # Kludge to make end of night summary plots

        print len(self.accomplished_field_ids)

        # Clean up
        self.accomplished_field_ids = [] 

############################################################

if __name__ == '__main__':

    my_simulator = Simulator('target_fields.txt')
    #my_simulator.selectField('2016/2/11 03:00', plot=True)
    
    my_simulator.loadObservationWindows('observation_windows.txt')
    #my_simulator.run('2016/2/11 03:00')
    #my_simulator.run('2016/2/11 08:00')
    my_simulator.run(plot=False)

############################################################
