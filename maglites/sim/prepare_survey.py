"""
Decide which fields to observe and time windows to observe.
"""

import os
import numpy as np
import ephem

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.ortho

############################################################

def prepareObservationWindows(nights, horizon=-14., outfile=None):
    """
    Use -14 deg twilight as default for start and end of observation windows.
    """

    observatory = ephem.Observer()
    observatory.lon = maglites.utils.constants.LON_CTIO
    observatory.lat = maglites.utils.constants.LAT_CTIO
    observatory.elevation = maglites.utils.constants.ELEVATION_CTIO
    observatory.horizon = '%.2f'%(horizon)

    sun = ephem.Sun()

    observation_windows = []
    for date, mode in nights:
        observatory.date = '%s 03:00'%(date)
        observatory.date = observatory.date + 24. * ephem.hour 

        time_setting = ephem.Date(observatory.previous_setting(sun))
        time_rising = ephem.Date(observatory.next_rising(sun))
        time_midpoint = ephem.Date(0.5 * (time_setting + time_rising))

        if mode == 'full':
            observation_windows.append([time_setting, time_rising])
        elif mode == 'first':
            observation_windows.append([time_setting, time_midpoint])
        elif mode == 'second':
            observation_windows.append([time_midpoint, time_rising])

        #print date
        #print observation_windows[-1]
        #print type(observation_windows[-1][0])

    if outfile:
        writer = open(outfile, 'w')
        for observation_window in observation_windows:
            writer.write('%s, %s\n'%(observation_window[0], observation_window[1]))
        writer.close()

    return observation_windows

############################################################

def prepareTargetList(infile, outfile=None):
    # Consider to remove SMASH fields?
    # How to create a dither??

    #data = np.recfromtxt('smash_fields_alltiles.txt', names=['RA', 'DEC'])
    data = np.recfromtxt('%s/maglites/data/smash_fields_alltiles.txt'%(os.environ['MAGLITESDIR']), names=True)
    ra = data['RA']
    dec = data['DEC']
    l, b = maglites.utils.projector.celToGal(ra, dec)

    angsep_lmc = maglites.utils.projector.angsep(maglites.utils.constants.RA_LMC, maglites.utils.constants.DEC_LMC, ra, dec)
    angsep_smc = maglites.utils.projector.angsep(maglites.utils.constants.RA_SMC, maglites.utils.constants.DEC_SMC, ra, dec)
    cut = (np.fabs(b) > 10.) \
          & ((angsep_lmc < 30.) | (angsep_smc < 30.)) \
          & (dec < -55.) & (ra > 100.) & (ra < 300.)
    #cut = cut | ((dec < -65.) & (angsep_lmc > 5.) & (angsep_smc > 5.))
    cut = cut | ((dec < -65.) & (ra > 300.) & (ra < 360.)) # SMC
    cut = cut | (dec < -80.)

    print np.sum(cut)

    ra_select = ra[cut]
    dec_select = dec[cut]
    tiling = np.tile(1, np.sum(cut))
    priority = np.tile(1, np.sum(cut))
    field_id = np.arange(1, np.sum(cut) + 1)

    # Kludge to make 3 effective tilings
    ra_select = np.tile(ra_select, 3)
    dec_select = np.tile(dec_select, 3)
    tiling = np.tile(np.arange(1, 4), np.sum(cut)).reshape(np.sum(cut), 3).transpose().flatten()
    priority = np.tile(1, 3 * np.sum(cut))
    field_id = np.arange(1, (3 * np.sum(cut)) + 1)

    #fig, ax, basemap = maglites.utils.ortho.makePlot('2016/2/11 03:00')
    fig, basemap = maglites.utils.ortho.makePlot('2016/2/11 03:00')

    proj = maglites.utils.ortho.safeProj(basemap, ra_select, dec_select)
    basemap.scatter(*proj, color='orange', edgecolor='none', s=50)

    if outfile:
        np.savetxt(outfile, zip(field_id, ra_select, dec_select, tiling, priority), 
                   fmt='%12i%12.4f%12.4f%12i%12i', 
                   header='%10s%12s%12s%12s%12s'%('ID', 'RA', 'DEC', 'TILING', 'PRIORITY'))
        #data2 = np.recfromtxt(outfile, names=True)
        #print data2

    #return data, data2
    #return data

############################################################

if __name__ == '__main__':

    nights = [['2016/2/10', 'second'],
              ['2016/2/11', 'second'],
              ['2016/2/12', 'second'],
              ['2016/2/13', 'second'],
              ['2016/2/14', 'second'],
              ['2016/2/15', 'second'],
              ['2016/6/27', 'full'],
              ['2016/6/28', 'full'],
              ['2016/6/29', 'full']]

    observation_windows = prepareObservationWindows(nights, outfile='observation_windows.txt')

    #data, data2 = prepareTargetList('smash_fields_alltiles.txt', outfile='list.txt')
    prepareTargetList('%s/maglites/data/smash_fields_alltiles.txt'%(os.environ['MAGLITESDIR']), outfile='target_fields.txt')

############################################################
