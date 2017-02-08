#!/usr/bin/env python
"""
Scheduling tactician.
"""
import os
import logging

import numpy as np
import ephem
from collections import OrderedDict as odict

from obztak.utils.projector import angsep

from obztak.utils import projector as proj
from obztak.ctio import CTIO
from obztak.utils import constants
from obztak.utils.date import datestring

CONDITIONS = odict([
    ('great',   [1.4, 2.0]),
    ('good',    [0.0, 2.0]),
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
        # Set previous field as last completed field
        previous_field = None
        if (self.completed_fields is not None) and len(self.completed_fields):
            previous_field = self.completed_fields[-1]

            # Ignore if more than 30 minutes has elapsed
            if (self.date-ephem.Date(previous_field['DATE'])) > 30*ephem.minute:
                previous_field = None

        if previous_field:
            return angsep(previous_field['RA'],previous_field['DEC'],
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
            # WARNING: This is broken when selecting two fields at once
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

    def __init__(self, *args, **kwargs):
        super(ConditionTactician,self).__init__(*args,**kwargs)
        self.mode = kwargs.get('mode',None)

    @property
    def weight(self):
        airmass = self.airmass
        sel = self.viable_fields
        weight = 2.0 * self.hour_angle
        weight[~sel] = np.inf
        weight += 3. * 360. * self.fields['TILING']
        weight += self.slew**3
        airmass_min, airmass_max = CONDITIONS[self.mode]
        airmass_cut = ((airmass < airmass_min) | (airmass > airmass_max))

        # ADW: This should probably also be in there
        weight += 100. * (airmass - 1.)**3
        weight += 5000. * airmass_cut

        return weight

class SMCNODTactician(Tactician):

    @property
    def weight(self):
        sel = self.viable_fields
        weight = 10000. * np.logical_not(np.in1d(self.fields['HEX'], obztak.utils.constants.HEX_SMCNOD)).astype(float)
        weight[~sel] = np.inf
        weight += 360. * self.fields['TILING']
        weight += slew
        return weight


class BlissTactician(Tactician):
    CONDITIONS = odict([
        (None,    [0.0, 1.4]),
        ('bliss', [0.0, 1.4]),
        ('good',  [0.0, 1.4]),
        ('poor',  [0.0, 1.2]),
    ])

    def __init__(self, *args, **kwargs):
        super(BlissTactician,self).__init__(*args,**kwargs)
        self.mode = kwargs.get('mode',None)

    @property
    def weight(self):
        airmass = self.airmass
        moon_angle = self.moon_angle

        sel = self.viable_fields
        weight = np.zeros(len(sel))

        # Moon angle constraints
        moon_limit = 35.
        sel &= (moon_angle > moon_limit)

        # Moon band constraints
        if (self.moon.phase >= 80) and (self.moon.alt > -0.04):
            # Moon is very bright; only do z
            sel &= (np.char.count('z',self.fields['FILTER']) > 0)
        elif (self.moon.phase >= 45) and (self.moon.alt > -0.04):
            # Moon is more than half full; do i,z
            sel &= (np.char.count('iz',self.fields['FILTER']) > 0)
        else:
            # Moon is faint or down; do g,r (unless none available)
            sel &= (np.char.count('gr',self.fields['FILTER']) > 0)
            #weight += 1e8 * (np.char.count('iz',self.fields['FILTER']) > 0)
        if (self.sun.alt > -0.28):
            # No g-band if Sun altitude > -16 deg
            sel &= ~(np.char.count('g',self.fields['FILTER']) > 0)

        # Airmass cut
        airmass_min, airmass_max = self.CONDITIONS[self.mode]
        sel &= ((airmass > airmass_min) & (airmass < airmass_max))

        # Don't allow the same field to be scheduled in different bands
        # less than 8 hours apart
        if len(self.completed_fields):
            dates = np.array(map(ephem.Date,self.completed_fields['DATE']))
            recent = self.completed_fields[(self.date - dates) < 10*ephem.hour]
            sel &= ~np.in1d(self.fields.field_id,recent.field_id)
            #cut = np.in1d(self.fields.field_id,recent.field_id)
            #sel &= ~cut

        # Set the weights for each field. Lower weight means more favorable.

        # Higher weight for rising fields (higher hour angle)
        # HA [min,max] = [-53,54] (for airmass 1.4)
        weight = 5.0 * self.hour_angle
        #weight = 1.0 * self.hour_angle

        # Higher weight for larger slews
        # slew = 10 deg -> weight = 1e2
        weight += self.slew**2

        # Higher weight for higher airmass
        # airmass = 1.4 -> weight = 6.4
        weight += 100. * (airmass - 1.)**3

        # Higher weight for fields close to the moon
        # angle = 50 -> weight = 6.4
        weight += 100 * (35./moon_angle)**3

        # Try hard to do the first tiling
        weight += 1e6 * (self.fields['TILING'] - 1)

        # Prioritize Planet 9 Region late in the survey/night
        # Allow i,z exposures at high penalty
        ra_zenith, dec_zenith = np.degrees(self.observatory.radec_of(0,'90'))
        if ra_zenith > 270:
            weight += 1e6 * (self.fields['PRIORITY'] - 1)
            #sel &= (np.char.count('iz',self.fields['FILTER']) > 0)
            #weight += 1e8 * (np.char.count('iz',self.fields['FILTER']) > 0)

        # Set infinite weight to all disallowed fields
        weight[~sel] = np.inf

        return weight

    def select_index(self):
        weight = self.weight
        index = np.array([np.argmin(weight)],dtype=int)
        if np.any(~np.isfinite(weight[index])):
            msg = "Infinite weight selected"
            print(msg)
            import obztak.utils.ortho, pylab as plt
            airmass_min, airmass_max = self.CONDITIONS[self.mode]
            obztak.utils.ortho.plotFields(self.completed_fields[-1],self.fields,self.completed_fields,options_basemap=dict(airmass=airmass_max))
            import pdb; pdb.set_trace()
            raw_input()
            raise ValueError(msg)
        return index

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()



#        if mode == 'airmass':
#            airmass_effective = copy.copy(airmass)
#            # Do not observe fields that are unavailable
#            airmass_effective[np.logical_not(cut)] = np.inf
#            # Priorize coverage over multiple tilings
#            airmass_effective += self.target_fields['TILING']
#            index_select = np.argmin(airmass_effective)
#        elif mode == 'ra':
#            # Different selection
#            #ra_effective = copy.copy(self.target_fields['RA'])
#            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
#            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
#            ra_effective[np.logical_not(cut)] = np.inf
#            ra_effective += 360. * self.target_fields['TILING']
#            index_select = np.argmin(ra_effective)
#        elif mode == 'slew':
#            #ra_effective = copy.copy(self.target_fields['RA'])
#            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
#            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
#            ra_effective[np.logical_not(cut)] = np.inf
#            ra_effective += 360. * self.target_fields['TILING']
#            ra_effective += slew**2
#            #ra_effective += 2. * slew
#            index_select = np.argmin(ra_effective)
#        elif mode == 'balance':
#            """
#            ra_effective = copy.copy(self.target_fields['RA']) - ra_zenith
#            ra_effective[ra_effective > 180.] = ra_effective[ra_effective > 180.] - 360.
#            ra_effective[np.logical_not(cut)] = np.inf
#            ra_effective += 360. * self.target_fields['TILING']
#            #ra_effective += 720. * self.target_fields['TILING']
#            ra_effective += slew**2
#            ra_effective += 100. * (airmass - 1.)**3
#            weight = ra_effective
#            index_select = np.argmin(weight)
#            weight = hour_angle_degree
#            """
#            weight = copy.copy(hour_angle_degree)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 3. * 360. * self.target_fields['TILING']
#            weight += slew**3 # slew**2
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'balance2':
#            weight = copy.copy(hour_angle_degree)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 360. * self.target_fields['TILING']
#            weight += slew_ra**2
#            weight += slew_dec
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'balance3':
#            logging.debug("Slew: %s"%slew)
#            weight = copy.copy(hour_angle_degree)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 3. * 360. * self.target_fields['TILING']
#            """
#            x_slew, y_slew = zip(*[[0., 0.],
#                                   [2.5, 10.],
#                                   [5., 30.],
#                                   [10., 150.],
#                                   [20., 250.],
#                                   [50., 500.],
#                                   [180., 5000.]])
#            """
#            x_slew, y_slew = zip(*[[0., 0.],
#                                   [2.5, 10.],
#                                   [5., 30.],
#                                   [10., 500.], #
#                                   [20., 1000.], # 500
#                                   [50., 5000.], # 1000
#                                   [180., 5000.]])
#            weight += np.interp(slew, x_slew, y_slew, left=np.inf, right=np.inf)
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'airmass2':
#            weight = 200. * (airmass - airmass_next)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 360. * self.target_fields['TILING']
#            weight += 100. * (airmass - 1.)**3
#            weight += slew**2
#            index_select = np.argmin(weight)
#        elif mode in ('coverage','good'):
#            weight = copy.copy(hour_angle_degree)
#            #weight[np.logical_not(cut)] = 9999.
#            weight[np.logical_not(cut)] = np.inf
#            weight += 6. * 360. * self.target_fields['TILING'] # Was 6, 60
#            weight += slew**3 # slew**2
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'coverage2':
#            weight = copy.copy(hour_angle_degree)
#            weight *= 2.
#            weight[np.logical_not(cut)] = np.inf
#            weight += 6. * 360. * self.target_fields['TILING']
#            weight += slew**3 # slew**2
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'coverage3':
#            weight = copy.copy(hour_angle_degree)
#            weight *= 0.5
#            weight[np.logical_not(cut)] = np.inf
#            weight += 6. * 360. * self.target_fields['TILING']
#            weight += slew**3 # slew**2
#            weight += 100. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#        elif mode == 'lowairmass':
#            weight = 2.0 * copy.copy(hour_angle_degree)
#            #if len(self.scheduled_fields) == 0:
#            #    weight += 200. * obztak.utils.projector.angsep(self.target_fields['RA'],
#            #                                                     self.target_fields['DEC'],
#            #                                                     90., -70.)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 3. * 360. * self.target_fields['TILING']
#            weight += slew**3 # slew**2
#            #weight += 2000. * (airmass - 1.)**3 # 200
#            weight += 5000. * (airmass > 1.5)
#            index_select = np.argmin(weight)
#
#            """
#            weight = copy.copy(hour_angle_degree)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 3. * 360. * self.target_fields['TILING']
#            weight += slew**3 # slew**2
#            weight += 1000. * (airmass - 1.)**3
#            index_select = np.argmin(weight)
#            """
#        elif mode in CONDITIONS.keys():
#            weight = 2.0 * copy.copy(hour_angle_degree)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 3. * 360. * self.target_fields['TILING']
#            weight += slew**3
#            airmass_min, airmass_max = CONDITIONS[mode]
#            airmass_sel = ((airmass < airmass_min) | (airmass > airmass_max))
#            # ADW: This should probably also be in there
#            weight += 100. * (airmass - 1.)**3
#            weight += 5000. * airmass_sel
#            index_select = np.argmin(weight)
#        elif mode == 'smcnod':
#            weight = 10000. * np.logical_not(np.in1d(self.target_fields['HEX'], obztak.utils.constants.HEX_SMCNOD)).astype(float)
#            weight[np.logical_not(cut)] = np.inf
#            weight += 360. * self.target_fields['TILING']
#            weight += slew
#            index_select = np.argmin(weight)
#        else:
#            msg = "Unrecognized mode: %s"%mode
#            raise Exception(msg)
