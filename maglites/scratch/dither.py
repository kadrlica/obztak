import os
import numpy as np
import pylab
import matplotlib.path
from matplotlib.collections import PolyCollection

import maglites.utils.projector
import maglites.utils.fileio as fileio
import maglites.utils.constants

pylab.ion()

############################################################

params = {
    #'backend': 'eps',
    'axes.labelsize': 16,
    #'text.fontsize': 12,           
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'xtick.major.size': 3,      # major tick size in points
    'xtick.minor.size': 1.5,    # minor tick size in points
    'text.usetex': True,
    #'figure.figsize': fig_size,
    'font.family':'serif',
    'font.serif':'Computer Modern Roman',
    'font.size': 12
    }

matplotlib.rcParams.update(params)

############################################################

def rotateFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field):
    proj_center = maglites.utils.projector.Projector(ra_center, dec_center)
    proj_field = maglites.utils.projector.Projector(ra_field, dec_field)
    
    ccd_array_new = []
    for ii in range(0, len(ccd_array)):
        ra, dec = proj_field.imageToSphere(np.transpose(ccd_array[ii])[0], 
                                           np.transpose(ccd_array[ii])[1]) 
        x, y = proj_center.sphereToImage(ra, dec)
        ccd_array_new.append(zip(x, y))

    return ccd_array_new

############################################################

def plotFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field, ax):
    ccd_array_new = rotateFocalPlane(ccd_array, ra_center, dec_center, ra_field, dec_field)
    coll = PolyCollection(ccd_array_new, alpha=0.2, color='red', edgecolors='none')
    ax.add_collection(coll)

############################################################

def applyDither(ra, dec, x, y):
    proj = maglites.utils.projector.Projector(ra, dec)
    ra_dither, dec_dither = proj.imageToSphere(x, y)
    return ra_dither, dec_dither

############################################################

def makeDither():

    X_CCD = 0.29878 # This is the FITS y-axis
    Y_CCD = 0.14939 # This is the FITS x-axis

    #ra_center, dec_center = 182., -88.0
    ra_center, dec_center = 182., -68.0
    ra_center, dec_center = 178., -80.0
    #ra_center, dec_center = 351.6670, -72.0863
    #ra_center, dec_center = 351.6670, -89.

    pattern = 'alex'
    if pattern == 'none':
        dither_array = []
    elif pattern == 'large':
        dither_array = [[4 * X_CCD / 3., 4. * Y_CCD / 3.],
                        [8. * X_CCD / 3., -11 * Y_CCD / 3.]]
    elif pattern == 'small':
        dither_array = [[1 * X_CCD / 3., 1. * Y_CCD / 3.],
                        [2. * X_CCD / 3., -1 * Y_CCD / 3.]]
    elif pattern == 'alex':
        ### ADW: The pattern suggested is actually in SMASH coordinates not celestial.
        dither_array = [[0.75, 0.75],
                        [-0.75, 0.75]]
    #dither_array = [[4 * X_CCD / 3., 4. * Y_CCD / 3.],
    #                [2. * X_CCD / 3., -4 * Y_CCD / 3.]]
    #dither_array = [[4 * X_CCD / 3., 4. * Y_CCD / 3.]]
    #dither_array = [[5 * X_CCD / 3., 5. * Y_CCD / 3.]]
    #mode = 'single' 
    mode = 'fill' 
    if mode == 'single':
        angsep_max = 0.
    if mode == 'fill':
        angsep_max = 3.

    # This should use the environment variable MAGLITESDIR to define the path
    datadir = fileio.get_datadir()
    filename  = os.path.join(datadir,'smash_fields_alltiles.txt')
    data_alltiles = np.recfromtxt(filename, names=True)

    filename = os.path.join(datadir,'ccd_corners_xy_fill.dat')
    data = eval(''.join(open(filename).readlines()))
    ccd_array = []
    for key in data.keys():
        #ccd_array.append(matplotlib.path.Path(data[key]))
        ccd_array.append(data[key])

    """
    n = 400
    x_mesh, y_mesh = np.meshgrid(np.linspace(-1.1, 1.1, n), np.linspace(-1.1, 1.1, n))

    count = np.zeros([n, n])

    for ii in range(0, len(ccd_array)):
    count += ccd_array[ii].contains_points(zip(x_mesh.flatten(), y_mesh.flatten())).reshape([n, n])

    pylab.figure()
    pylab.pcolor(x_mesh, y_mesh, count)
    """

    fig, ax = pylab.subplots(figsize=(8, 8))

    # Make the collection and add it to the plot.
    #coll = PolyCollection(ccd_array, alpha=0.3, color='red', edgecolors='none')
    #ax.add_collection(coll)

    #plotFocalPlane(ccd_array, ra_center, dec_center, ra_center, dec_center, ax)
    #plotFocalPlane(ccd_array, ra_center, dec_center, ra_center, dec_center + 0.1, ax)

    angsep = maglites.utils.projector.angsep(ra_center, dec_center, data_alltiles['RA'], data_alltiles['DEC'])
    for ii in np.nonzero(angsep < (np.min(angsep) + 0.01 + angsep_max))[0]:
        plotFocalPlane(ccd_array, ra_center, dec_center, data_alltiles['RA'][ii], data_alltiles['DEC'][ii], ax)
    
        for x_dither, y_dither in dither_array:
            ra_dither, dec_dither = applyDither(data_alltiles['RA'][ii], data_alltiles['DEC'][ii], 
                                                x_dither, y_dither)
        plotFocalPlane(ccd_array, ra_center, dec_center, ra_dither, dec_dither, ax)

    pylab.xlim(-1.5, 1.5)
    pylab.ylim(-1.5, 1.5)
    pylab.xlabel('x (deg)', labelpad=20)
    pylab.ylabel('y (deg)')
    pylab.title('(RA, Dec) = (%.3f, %.3f)'%(ra_center, dec_center))

    #for ii in range(0, len(ccd_array)):
    #pylab.savefig('dither_ra_%.2f_dec_%.2f_%s_%s.pdf'%(ra_center, dec_center, mode, pattern))

############################################################

def testDither(ra_center, dec_center, infile='target_fields.csv', save=False):

    filename = os.path.join(fileio.get_datadir(),'ccd_corners_xy_fill.dat')
    data = eval(''.join(open(filename).readlines()))
    ccd_array = []
    for key in data.keys():
        #ccd_array.append(matplotlib.path.Path(data[key]))
        ccd_array.append(data[key])

    data_targets = fileio.csv2rec(infile)

    fig, ax = pylab.subplots(figsize=(8, 8))
    
    angsep = maglites.utils.projector.angsep(ra_center, dec_center, data_targets['RA'], data_targets['DEC'])
    cut = (angsep < 3.) & (data_targets['FILTER'] == maglites.utils.constants.BANDS[0]) & (data_targets['TILING'] <= 3)

    print np.sum(angsep < 3.)
    print np.sum(data_targets['FILTER'] == maglites.utils.constants.BANDS[0])
    print np.sum(cut)

    for ii in np.nonzero(cut)[0]:
        plotFocalPlane(ccd_array, ra_center, dec_center, data_targets['RA'][ii], data_targets['DEC'][ii], ax)
    
    pylab.xlim(-1.5, 1.5)
    pylab.ylim(-1.5, 1.5)
    pylab.xlabel('x (deg)', labelpad=20)
    pylab.ylabel('y (deg)')
    pylab.title('(RA, Dec) = (%.3f, %.3f)'%(ra_center, dec_center))

    if save:
        pattern = infile.split('target_fields_')[-1].split('.csv')[0]
        pylab.savefig('dither_ra_%.2f_dec_%.2f_%s.pdf'%(ra_center, dec_center, pattern))

############################################################

#infile = 'target_fields_decam_dither_1.csv'
#infile = 'target_fields_decam_dither_2.csv'
infile = 'target_fields_smash_dither.csv'
#infile = 'target_fields_smash_rotate.csv'
save = True

testDither(100., -70., infile=infile, save=save) # On edge
testDither(125., -75., infile=infile, save=save)
testDither(200., -88., infile=infile, save=save)

############################################################
