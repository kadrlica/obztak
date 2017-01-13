#!/usr/bin/env python
"""
Scheduling tactician.
"""
import os

import numpy as np
import ephem
from collections import OrderedDict as odict

from obztak.utils.projector import angsep

from obztak.utils import projector as proj
from obztak.ctio import CTIO
from obztak.utils import constants
from obztak.utils.date import datestring

CONDITIONS = odict([
    ('great', [1.4, 2.0]),
    ('good',  [0.0, 2.0]),
    ('fine',  [0.0, 1.9]),
    ('ok',    [0.0, 1.6]),
    ('poor',  [0.0, 1.5]),
    ('bad',   [0.0, 1.4]),
])

class Tactician(object):
    name = 'tactician'

    def __init__(self,fields,observatory=None):
        """ Initialize the scheduling tactician.

        Parameters:
        -----------
        fields      : The available fields.
        observatory : The observatory (defaults to CTIO)
        
        Returns:
        --------
        Tactician   : The Tactician object
        """
        self.fields = fields.copy()
        if not observatory: observatory = CTIO()
        self.observatory = observatory
        self.moon = ephem.Moon()
        self.previous_field = None

    def set_date(self,date):
        if date is not None:
            self.observatory.date = ephem.Date(date)
            self.moon.compute(self.observatory.date)
        return self.observatory.date

    def set_previous_field(self,field):
        if field is not None:
            self.previous_field = np.copy(field)
        else:
            self.previous_field = None

    @property
    def date(self):
        return self.observatory.date

    @property
    def hour_angle_limit(self):
        return self.observatory.hour_angle_limit(self.fields['DEC']) 

    @property
    def airmass_limit(self):
        return self.observatory.airmass_limit(self.fields['DEC'])

    @property
    def zenith_angle(self):
        # RA and Dec of zenith
        return np.degrees(self.observatory.radec_of(0,'90'))

    @property
    def airmass(self):
        """ Calculate the airmass of each field. """
        ra_zenith,dec_zenith = self.zenith_angle
        return proj.airmass(ra_zenith, dec_zenith, 
                            self.fields['RA'], self.fields['DEC'])
    @property
    def moon_angle(self):
        # Include moon angle
        ra_moon,dec_moon = np.degrees([self.moon.ra,self.moon.dec])
        return proj.angsep(ra_moon, dec_moon, 
                           self.fields['RA'], self.fields['DEC'])
    @property
    def moon_phase(self):
        return self.moon.phase

    @property
    def slew(self):
        if self.previous_field:
            return angsep(self.previous_field['RA'],self.previous_field['DEC'],
                          self.fields['RA'], self.fields['DEC'])
        else:
            return np.zeros(len(self.fields))

    @property
    def hour_angle(self):
        ra_zenith,dec_zenith = self.zenith_angle
        hour_angle = np.copy(self.fields['RA']) - ra_zenith
        hour_angle[hour_angle < -180.] += 360.
        hour_angle[hour_angle > 180.] -= 360.
        return hour_angle

    @property
    def viable_fields(self):
        # Check the hour angle restrictions at south pole
        sel_hour_angle = np.fabs(self.hour_angle) < self.hour_angle_limit

        # Blanco airmass restrictions
        sel_airmass = self.airmass < self.airmass_limit

        # Declination restrictions
        sel_declination = self.fields['DEC'] > constants.SOUTHERN_REACH

        # Exclude special fields (unless using special tacticians)
        sel_special = self.fields['PRIORITY'] < 90

        viable = sel_hour_angle & sel_airmass & sel_declination & sel_special
        return viable

    @property
    def weight(self):
        weight = self.hour_angle
        sel = self.viable_fields
        weight[~sel] = np.inf
        weight += 6. * 360. * self.fields['TILING'] # Was 6, 60
        weight += self.slew**3 # slew**2
        weight += 100. * (self.airmass - 1.)**3
        return weight

    def select_index(self):
        index_select = np.argmin(self.weight)

        # Search for other exposures in the same field
        field_id = self.fields['HEX'][index_select]
        tiling   = self.fields['TILING'][index_select]

        index = np.nonzero( (self.fields['HEX'] == field_id) & 
                            (self.fields['TILING'] == tiling))[0]
        return index

    def select_fields(self):
        index = self.select_index()

        timedelta = constants.FIELDTIME*np.arange(len(index))
        if np.any(self.slew[index] > 5.):
            # Apply a 30 second penalty for slews over 5 deg.
            # This is not completely realistic, but better than nothing
            # This is also broken when selecting two fields at once
            timedelta += 30*ephem.second

        fields              = self.fields[index]
        fields['AIRMASS']   = self.airmass[index]
        fields['DATE']      = map(datestring,self.date+timedelta)
        fields['SLEW']      = self.slew[index]
        fields['MOONANGLE'] = self.moon_angle[index]
        fields['HOURANGLE'] = self.hour_angle[index]
        return fields

class CoverageTactician(Tactician):
    name = 'coverage'

    @property
    def weight(self):
        sel = self.viable_fields
        weight = self.hour_angle
        weight[~sel] = np.inf
        weight += 6. * 360. * self.fields['TILING'] # Was 6, 60
        weight += self.slew**3 # slew**2
        weight += 100. * (self.airmass - 1.)**3
        return weight

class ConditionTactician(Tactician):
    name = 'condition'

    @property
    def weight(self):
        airmass = self.airmass
        sel = self.viable_fields
        weight = 2.0 * self.hour_angle
        weight[~sel] = np.inf
        weight += 3. * 360. * self.fields['TILING']
        weight += self.slew**3
        airmass_min, airmass_max = CONDITIONS[mode]
        airmass_cut = ((airmass < airmass_min) | (airmass > airmass_max))

        # ADW: This should probably also be in there
        weight += 100. * (airmass - 1.)**3
        weight += 5000. * airmass_cut

        return weight


### class AirmassTactician(Tactician):
###     name = 'airmass'
###  
###     @property
###     def weight(self,airmass):
###         airmass_effective = copy.copy(airmass)
###         # Do not observe fields that are unavailable
###         airmass_effective[np.logical_not(cut)] = np.inf 
###         # Priorize coverage over multiple tilings
###         airmass_effective += self.target_fields['TILING'] 
###         index_select = np.argmin(airmass_effective)
###         return index_select
###  
### class RATactician(Tactician):
###     name = 'ra'
###     def calculate_weight(self,ra_zenith):
###         ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
###         ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
###         ra_effective[np.logical_not(cut)] = np.inf
###         ra_effective += 360. * self.target_fields['TILING']
###         index_select = np.argmin(ra_effective)
###         return index_select
###  
### class SlewTactician(Tactician):
###     name = 'slew'
###     def calculate_weights(self,fields,ra_zenith,slew):
###         ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
###         ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
###         ra_effective[np.logical_not(cut)] = np.inf
###         ra_effective += 360. * self.target_fields['TILING']
###         ra_effective += slew**2
###         #ra_effective += 2. * slew
###         index_select = np.argmin(ra_effective)
###         return index_select
###  
### class BalanceTactician(Tactician):
###     name = 'balance'
###     def calculate_weight(self,airmass, hour_angle_degree):
###             """
###             ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
###             ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
###             ra_effective[np.logical_not(cut)] = np.inf
###             ra_effective += 360. * self.target_fields['TILING']
###             #ra_effective += 720. * self.target_fields['TILING']
###             ra_effective += slew**2
###             ra_effective += 100. * (airmass - 1.)**3
###             weight = ra_effective
###             index_select = np.argmin(weight)
###             weight = hour_angle_degree
###             """
###             weight = copy.copy(hour_angle_degree)
###             weight[np.logical_not(cut)] = np.inf
###             weight += 3. * 360. * self.target_fields['TILING']
###             weight += slew**3 # slew**2
###             weight += 100. * (airmass - 1.)**3
###             index_select = np.argmin(weight)
###             return index_select
###  
### class BalanceTactician2(Tactician):
###     name = 'balance2'
###     def calculate_weight(self,hour_angle_degree,slew_ra,slew_dec):
###         weight = copy.copy(hour_angle_degree)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 360. * self.target_fields['TILING']
###         weight += slew_ra**2
###         weight += slew_dec
###         weight += 100. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###         return index_select
###  
### class BalanceTactician2(Tactician):
###     name = 'balance3'
###     def calculate_weight(self,hour_angle_degree,slew_ra,slew_dec):
###         logging.debug("Slew: %s"%slew)
###         weight = copy.copy(hour_angle_degree)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 3. * 360. * self.target_fields['TILING']
###         """
###         x_slew, y_slew = zip(*[[0., 0.],
###                                [2.5, 10.],
###                                [5., 30.],
###                                [10., 150.],
###                                [20., 250.],
###                                [50., 500.],
###                                [180., 5000.]])
###         """
###         x_slew, y_slew = zip(*[[0., 0.],
###                                [2.5, 10.],
###                                [5., 30.],
###                                [10., 500.], #
###                                [20., 1000.], # 500
###                                [50., 5000.], # 1000
###                                [180., 5000.]])
###         weight += np.interp(slew, x_slew, y_slew, left=np.inf, right=np.inf)
###         weight += 100. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###         return index_select
###  
### class AirmassTactician2(Tactician):
###     name = 'airmass2'
###     def calculate_weight(self,hour_angle_degree,slew_ra,slew_dec):
###         weight = 200. * (airmass - airmass_next)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 360. * self.target_fields['TILING']
###         weight += 100. * (airmass - 1.)**3
###         weight += slew**2
###         index_select = np.argmin(weight)
###         return index_select
###  
### class CoverageTactician(Tactician):
###     name = 'coverage'
###     def calculate_weight(self,airmass,hour_angle_degree,slew):
###         weight = copy.copy(hour_angle_degree)
###         #weight[np.logical_not(cut)] = 9999.
###         weight[np.logical_not(cut)] = np.inf
###         weight += 6. * 360. * self.target_fields['TILING'] # Was 6, 60
###         weight += slew**3 # slew**2
###         weight += 100. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###         return index_select
###  
### class CoverageTactician2(Tactician):
###     name = 'converage2'
###     def calculate_weight(self,hour_angle_degree,slew_ra,slew_dec):
###         weight = copy.copy(hour_angle_degree)
###         weight *= 2.
###         weight[np.logical_not(cut)] = np.inf
###         weight += 6. * 360. * self.target_fields['TILING']
###         weight += slew**3 # slew**2
###         weight += 100. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###         return index_select
###  
### class CoverageTactician3(Tactician):
###     name = 'coverage3'
###     def calculate_weight(self,airmass,hour_angle_degree,slew):
###         weight = copy.copy(hour_angle_degree)
###         weight *= 0.5
###         weight[np.logical_not(cut)] = np.inf
###         weight += 6. * 360. * self.target_fields['TILING']
###         weight += slew**3 # slew**2
###         weight += 100. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###  
### class LowAirmassTactician(Tactician):
###     name = 'lowairmass'
###     def calculate_weight(self):
###         weight = 2.0 * copy.copy(hour_angle_degree)
###         #if len(self.scheduled_fields) == 0:
###         #    weight += 200. * obztak.utils.projector.angsep(self.target_fields['RA'],
###         #                                                     self.target_fields['DEC'],
###         #                                                     90., -70.)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 3. * 360. * self.target_fields['TILING']
###         weight += slew**3 # slew**2
###         #weight += 2000. * (airmass - 1.)**3 # 200
###         weight += 5000. * (airmass > 1.5)
###         index_select = np.argmin(weight)
###  
###         """
###         weight = copy.copy(hour_angle_degree)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 3. * 360. * self.target_fields['TILING']
###         weight += slew**3 # slew**2
###         weight += 1000. * (airmass - 1.)**3
###         index_select = np.argmin(weight)
###         """
###  
### class ConditionTactician(Tactician):
###     name = 'condition'
###     def calculate_weight(self):
###         weight = 2.0 * copy.copy(hour_angle_degree)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 3. * 360. * self.target_fields['TILING']
###         weight += slew**3
###         airmass_min, airmass_max = CONDITIONS[mode]
###         airmass_sel = ((airmass < airmass_min) | (airmass > airmass_max))
###         # ADW: This should probably also be in there
###         weight += 100. * (airmass - 1.)**3
###         weight += 5000. * airmass_sel
###         index_select = np.argmin(weight)
###  
### class SmcnodTactician(Tactician):
###     name = 'smcnod'
###     def calculate_weight(self):
###         weight = 10000. * np.logical_not(np.in1d(self.target_fields['HEX'], obztak.utils.constants.HEX_SMCNOD)).astype(float)
###         weight[np.logical_not(cut)] = np.inf
###         weight += 360. * self.target_fields['TILING']
###         weight += slew
###         index_select = np.argmin(weight)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
