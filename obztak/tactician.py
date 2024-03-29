#!/usr/bin/env python
"""
Scheduling tactician.
"""
import os
import logging

import numpy as np
import ephem
from collections import OrderedDict as odict

from obztak.utils import projector as proj
from obztak.ctio import CTIO
from obztak.utils import constants
from obztak.utils.date import datestring

CONDITIONS = odict([
    (None,      [0.0, 2.0]), #default
    ('great',   [1.6, 2.0]),
    ('good',    [0.0, 2.0]),
    ('complete',[0.0, 2.0]),
    ('coverage',[0.0, 2.0]),
    ('maglites',[0.0, 2.0]),
    ('fine',    [0.0, 1.9]),
    ('ok',      [0.0, 1.6]),
    ('poor',    [0.0, 1.5]),
    ('bad',     [0.0, 1.4]),
])

class Tactician(object):
    name = 'tactician'

    def __init__(self, fields=None, observatory=None, **kwargs):
        """ Initialize the survey scheduling tactician.

        Parameters:
        -----------
        fields      : The available fields.
        observatory : The observatory (defaults to CTIO)
        
        Returns:
        --------
        Tactician   : The Tactician object
        """
        if not observatory: observatory = CTIO()
        self.observatory = observatory
        self.moon = ephem.Moon()
        self.sun = ephem.Sun()

        self.set_target_fields(fields)
        self.set_completed_fields(None)
        self.set_date(None)

    def set_date(self,date):
        if date is not None:
            self.observatory.date = ephem.Date(date)
            self.moon.compute(self.observatory)
            self.sun.compute(self.observatory)

    def set_target_fields(self,fields):
        if fields is not None:
            self.fields = fields.copy()
        else:
            self.fields = None

    def set_completed_fields(self,fields):
        if fields is not None:
            self.completed_fields = fields.copy()
        else:
            self.completed_fields = None

    def set_previous_field(self,field):
        #if field is not None:
        #    self.previous_field = field.copy()
        #else:
        #    self.previous_field = None
        pass

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
        # See here for ra,dec details: http://rhodesmill.org/pyephem/radec
        ra_moon,dec_moon = np.degrees([self.moon.ra,self.moon.dec])
        return proj.angsep(ra_moon, dec_moon, 
                           self.fields['RA'], self.fields['DEC'])
    @property
    def moon_phase(self):
        return self.moon.phase

    @property
    def slew(self):
        """Angular separation to previous field."""
        # Set previous field as last completed field
        previous_field = None
        if (self.completed_fields is not None) and len(self.completed_fields):
            previous_field = self.completed_fields[-1]

            # Ignore if more than 30 minutes has elapsed
            if (self.date-ephem.Date(previous_field['DATE'].astype(str))) > 30*ephem.minute:
                previous_field = None

            # Exposure being exposed has RA,DEC = 0,0; problem!
            # Could also require previous_field['PROGRAM'] == 'None'
            elif (previous_field['RA'] == 0) and (previous_field['DEC'] == 0):
                previous_field = None

        if previous_field:
            return proj.angsep(previous_field['RA'],previous_field['DEC'],
                               self.fields['RA'], self.fields['DEC'])
        else:
            return np.zeros(len(self.fields))

    @property
    def slew_time(self):
        """Estimate of the slew time (Alt/Az telescope)."""
        DeprecationWarning("'slew_time' has been deprecated; use 'slew'")
        # Set previous field as last completed field
        previous_field = None
        if (self.completed_fields is not None) and len(self.completed_fields):
            previous_field = self.completed_fields[-1]

            # Ignore if more than 30 minutes has elapsed
            if (self.date-ephem.Date(previous_field['DATE'].astype(str))) > 30*ephem.minute:
                previous_field = None

        if previous_field:
            return np.sqrt((previous_field['RA']-self.fields['RA'])**2 +
                           (previous_field['DEC']-self.fields['DEC'])**2)

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

        # Moon angle constraints
        moon_limit = 20.
        sel_moon = (self.moon_angle > moon_limit)

        viable = sel_hour_angle & sel_airmass & sel_declination & sel_special & sel_moon
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
        fields = self.fields[index]

        # This is janky since MagLiteS was originally scheduling pairs of exposures
        timedelta = np.arange(len(index),dtype=float)
        timedelta *= (fields['EXPTIME']*ephem.second + constants.OVERHEAD)
        timedelta[self.slew[index] > 5.] += 30*ephem.second

        #timedelta = constants.FIELDTIME*np.arange(len(index))
        #if np.any(self.slew[index] > 5.):
        #    # Apply a 30 second penalty for slews over 5 deg.
        #    # This is not completely realistic, but better than nothing
        #    # WARNING: This is broken when selecting two fields at once
        #    timedelta += 30*ephem.second

        fields['AIRMASS']   = self.airmass[index]
        fields['DATE']      = list(map(datestring,self.date+timedelta))
        fields['SLEW']      = self.slew[index]
        fields['MOONANGLE'] = self.moon_angle[index]
        fields['HOURANGLE'] = self.hour_angle[index]
        return fields

# DEPRECATED: ADW 2019-09-09
class CoverageTactician(Tactician):
    name = 'coverage'
    mode = None

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

    def __init__(self, *args, **kwargs):
        super(ConditionTactician,self).__init__(*args,**kwargs)
        self.mode = kwargs.get('mode',None)

    @property
    def coverage_weight(self):
        """ Hack to subsume CoverageTactician into ConditionTactician """
        sel = self.viable_fields
        weight = self.hour_angle
        weight[~sel] = np.inf
        weight += 6. * 360. * self.fields['TILING'] # Was 6, 60
        weight += self.slew**3 # slew**2
        weight += 100. * (self.airmass - 1.)**3
        return weight

    @property
    def weight(self):
        if self.mode == 'coverage':
            return self.coverage_weight

        airmass = self.airmass
        sel = self.viable_fields
        weight = 2.0 * self.hour_angle
        weight[~sel] = np.inf
        weight += 3. * 360. * self.fields['TILING']
        if self.mode == 'complete':
            weight += 100. * 360. * self.fields['TILING']

        weight += self.slew**3
        airmass_min, airmass_max = CONDITIONS[self.mode]
        airmass_cut = ((airmass < airmass_min) | (airmass > airmass_max))

        # ADW: This should probably also be in there
        weight += 100. * (airmass - 1.)**3
        weight += 5000. * airmass_cut

        if self.mode == 'great':
            weight += 5000. * (self.fields['DEC'] > -80)

        return weight


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
