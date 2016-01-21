"""
Create diagnostic plots for a particular survey strategy.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import ephem

import maglites.utils.ortho

plt.ion()

############################################################

def movie(infile_accomplished_fields, infile_target_fields=None, outdir=None, chunk=1):
    plt.ioff()

    accomplished_fields = np.recfromtxt(infile_accomplished_fields, delimiter=',', names=True)
    print len(accomplished_fields)

    if infile_target_fields is not None:
        target_fields = np.recfromtxt(infile_target_fields, names=True)
    
    for ii in range(0, 1000):
        print ii, accomplished_fields['DATE'][ii * chunk]
        fig, basemap = maglites.utils.ortho.makePlot(accomplished_fields['DATE'][ii * chunk], figsize=(10.5 / 2., 8.5 / 2.), s=25, dpi=160, moon=False)
        
        if infile_target_fields is not None:
            cut_accomplished = np.in1d(target_fields['ID'], accomplished_fields['ID'][0:ii * chunk])
            #cut_accomplished[ii] = True
            proj = maglites.utils.ortho.safeProj(basemap, 
                                                 target_fields['RA'][np.logical_not(cut_accomplished)],
                                                 target_fields['DEC'][np.logical_not(cut_accomplished)])
            basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=25, vmin=0, vmax=4, cmap='summer_r')
            
        proj = maglites.utils.ortho.safeProj(basemap, accomplished_fields['RA'][0:(ii * chunk) + 1], accomplished_fields['DEC'][0:(ii * chunk) + 1])
        basemap.scatter(*proj, c=accomplished_fields['TILING'][0:(ii * chunk) + 1], edgecolor='none', s=25, vmin=0, vmax=4, cmap='summer_r')
        colorbar = plt.colorbar(label='Tiling')

        # Show the selected field
        if chunk == 1:
            proj = maglites.utils.ortho.safeProj(basemap, [accomplished_fields['RA'][ii]], [accomplished_fields['DEC'][ii]])
        else:
            proj = maglites.utils.ortho.safeProj(basemap, 
                                                 accomplished_fields['RA'][ii * chunk:(ii + 1) * chunk], 
                                                 accomplished_fields['DEC'][ii * chunk:(ii + 1) * chunk])
        basemap.scatter(*proj, c='magenta', edgecolor='none', s=25)

        if outdir is not None:
            plt.savefig('%s/movie%05i.gif'%(outdir, ii))

        #raw_input('WAIT')
        fig.clf()

        if ephem.Date(accomplished_fields['DATE'][(ii + 1) * chunk]) > ephem.Date('2017/1/1 00:00'):
            break
    
    if outdir is not None:
        print 'Generating animated gif...'
        os.system('convert -set delay 10 -loop 0 %s/*.gif %s/output.gif'%(outdir, outdir))

############################################################

def slew(infile_accomplished_fields):
    accomplished_fields = np.recfromtxt(infile_accomplished_fields, delimiter=',', names=True)

    plt.figure()
    plt.hist(accomplished_fields['SLEW'], bins=np.linspace(0., 20., 21), color='green')
    plt.xlabel('Slew Angle (deg)')
    plt.ylabel('Number of Fields')
    plt.xlim(0., 20.)

    # Inset panel
    if np.any(accomplished_fields['SLEW'] > 20.):
        a = plt.axes([.4, .4, .45, .45], axisbg='w')
        max_slew = np.max(accomplished_fields['SLEW'])
        plt.hist(accomplished_fields['SLEW'], bins=np.arange(20., max_slew + 2., 1.), color='green')
        plt.xlabel('Slew Angle (deg)')
        plt.ylabel('Number of Fields')
        plt.xlim(20., max(40., max_slew + 2.))
    if save:
        plt.savefig('slew_hist.pdf')

    plt.figure()
    #plt.scatter(np.arange(len(accomplished_fields['SLEW'])), accomplished_fields['SLEW'], edgecolor='none', alpha=0.33)
    plt.scatter(np.arange(len(accomplished_fields['SLEW'])), accomplished_fields['SLEW'], marker='x')
    plt.xlabel('Sequential Field Observed')
    plt.ylabel('Slew Angle (deg)')
    plt.xlim(0., len(accomplished_fields['SLEW']) + 1)
    plt.ylim(0., 30.)
    if save:
        plt.savefig('slew_sequential.pdf')

    cut = accomplished_fields['SLEW'] > 10.
    for index in np.nonzero(cut)[0]:
        print accomplished_fields['SLEW'][index], accomplished_fields['RA'][index], accomplished_fields['DEC'][index], accomplished_fields['DATE'][index] 

############################################################

def airmass(infile_accomplished_fields):
    accomplished_fields = np.recfromtxt(infile_accomplished_fields, delimiter=',', names=True)

    plt.figure()
    plt.hist(accomplished_fields['AIRMASS'], bins=np.linspace(1., 2., 21), color='red')
    plt.xlabel('Airmass')
    plt.ylabel('Number of Fields')
    plt.xlim(1., 2.)
    if save:
        plt.savefig('airmass_hist.pdf')

    plt.figure()
    plt.scatter(np.arange(len(accomplished_fields['AIRMASS'])), accomplished_fields['AIRMASS'], marker='x', color='red')
    plt.xlabel('Sequential Field Observed')
    plt.ylabel('Airmass')
    plt.xlim(0., len(accomplished_fields['AIRMASS']) + 1)
    plt.ylim(1., 2.)
    if save:
        plt.savefig('airmass_sequential.pdf')

    fig, basemap = maglites.utils.ortho.makePlot(accomplished_fields['DATE'][0], figsize=(10.5, 8.5), s=50, dpi=80, center=(0., -90.), airmass=False, moon=False)
    proj = maglites.utils.ortho.safeProj(basemap, accomplished_fields['RA'], accomplished_fields['DEC'])
    basemap.scatter(*proj, c=accomplished_fields['AIRMASS'], edgecolor='none', s=50, alpha=0.5, vmin=np.min(accomplished_fields['AIRMASS']), vmax=2, cmap='RdYlGn_r')
    colorbar = plt.colorbar(label='Airmass')
    if save:
        plt.savefig('airmass_map.pdf')

############################################################

def progress(infile_accomplished_fields, date, infile_target_fields=None, save=False):
    if type(date) != ephem.Date:
        date = ephem.Date(date)
        
    accomplished_fields = np.recfromtxt(infile_accomplished_fields, delimiter=',', names=True)
    
    date_array = np.tile(0., len(accomplished_fields['DATE']))
    for ii in range(0, len(accomplished_fields['DATE'])):
        date_array[ii] = ephem.Date(accomplished_fields['DATE'][ii]).real
    cut_accomplished = date_array <= date.real 

    fig, basemap = maglites.utils.ortho.makePlot(date, figsize=(10.5, 8.5), s=50, dpi=80, airmass=False, moon=False, center=(0., -90.))
    
    if infile_target_fields is not None:
        target_fields = np.recfromtxt(infile_target_fields, names=True)
        proj = maglites.utils.ortho.safeProj(basemap, 
                                             target_fields['RA'][np.logical_not(cut_accomplished)],
                                             target_fields['DEC'][np.logical_not(cut_accomplished)])
        basemap.scatter(*proj, c=np.tile(0, np.sum(np.logical_not(cut_accomplished))), edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
            
    proj = maglites.utils.ortho.safeProj(basemap, accomplished_fields['RA'][cut_accomplished], accomplished_fields['DEC'][cut_accomplished])
    basemap.scatter(*proj, c=accomplished_fields['TILING'][cut_accomplished], edgecolor='none', s=50, vmin=0, vmax=4, cmap='summer_r')
    colorbar = plt.colorbar(label='Tiling')

    if save:
        plt.savefig('progress_map_%s.pdf'%(maglites.utils.ortho.datestring(date).replace('/', '_').replace(':', '_').replace(' ', '_')))

############################################################

if __name__ == '__main__':
    #movie('accomplished_fields.txt', infile_target_fields='target_fields.txt', outdir='movie')
    #movie('accomplished_fields.txt', infile_target_fields='target_fields.txt', outdir='movie2', chunk=2)
    #movie('accomplished_fields.txt', infile_target_fields='target_fields.txt', outdir='movie3', chunk=2)
    #slew('accomplished_fields.txt')
    #airmass('accomplished_fields.txt')
    progress('accomplished_fields.txt', '2016/6/30 10:32:50', infile_target_fields='target_fields.txt')
    #progress('accomplished_fields.txt', '2017/6/30 10:32:51', infile_target_fields='target_fields.txt')

############################################################
