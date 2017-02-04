"""
Decide which fields to observe and time windows to observe.
"""

import os,sys
import logging
import numpy.lib.recfunctions as recfunc
import numpy as np
import ephem

import obztak.utils.projector
import obztak.utils.ortho
import obztak.ctio

from obztak.utils import constants
from obztak.utils import fileio
from obztak.field import FieldArray

from obztak.utils.projector import cel2gal
from obztak.utils.date import datestring
from obztak.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS

class Survey(object):
    """Base class for preparing a survey. Creates a list of observation
    windows (observing dates and times) and target fields. The
    `prepare_windows` and `prepare_fields` methods should be
    subclassed for specific surveys.
    """

    # 2016A ACTUAL
    nights_2016A = [
        ['2016/02/11', 'second'],
        ['2016/02/12', 'second'],
        ['2016/02/13', 'second'],
        ['2016/02/14', 'second'],
        ['2016/02/15', 'second'],
        ['2016/02/16', 'second'],
        ['2016/06/27', 'full'],
        ['2016/06/28', 'full'],
        ['2016/06/29', 'full'],
        ]

    # 2017A ACTUAL
    nights_2017A = [
        ['2017/2/21', 'full'],
        ['2017/2/22', 'full'],
        ['2017/2/23', 'full'],
        ['2017/6/18', 'full'],
        ['2017/6/19', 'full'],
        ['2017/6/20', 'full']
        ]

    nights = nights_2016A + nights_2017A

    def prepare_windows(self, nights, horizon=-14., standards=True, outfile=None):
        """Create a list of observation windows consisting of start-stop
        times for all dates of observation.

        Parameters:
        -----------
        nights   : A list of the nights and modes ('full','first','second').
        horizon  : Sun altitude at the start of observing.
        standards: Add time for taking standards
        outfile  : Output file of observing windows.

        Returns:
        --------
        observation_windows : Record array of observation windows.
        """
        observatory = obztak.ctio.CTIO()
        observatory.horizon = str(horizon)

        sun = ephem.Sun()

        observation_windows = []
        for date, mode in nights:
            observatory.date = '%s 03:00'%(date)
            observatory.date = observatory.date + 24. * ephem.hour

            time_setting = ephem.Date(observatory.previous_setting(sun))
            time_rising = ephem.Date(observatory.next_rising(sun))
            time_midpoint = ephem.Date(0.5*(time_setting + time_rising))

            if mode == 'full':
                window = [time_setting, time_rising]
            elif mode == 'first':
                if standards:
                    time_midpoint = time_midpoint - STANDARDS
                window = [time_setting, time_midpoint]
            elif mode == 'second':
                if standards:
                    time_midpoint = time_midpoint + STANDARDS
                window = [time_midpoint, time_rising]
            else:
                msg = "Unrecognized mode: '%s'"%mode
                raise ValueError(msg)

            window = [datestring(w,0) for w in window]
            observation_windows.append(window)

        #dtype=[('UTC_START','S28'),('UTC_END','S28')]
        #observation_windows = np.rec.fromrecords(observation_windows,dtype=dtype)
        names=['UTC_START','UTC_END']
        observation_windows = np.rec.fromrecords(observation_windows,names=names)

        if outfile:
            fileio.rec2csv(outfile,observation_windows)

        return observation_windows

    def prepare_fields(self, infile=None, outfile=None, mode='smash_dither', plot=True, smcnod=False):
        """ Create the list of fields to be targeted by this survey.

        Parameters:
        -----------
        infile : File containing all possible field locations.
        outfile: Output file of selected fields
        mode   : Mode for dithering: 'smash_dither', 'smash_rotate', 'decam_dither', 'none'
        plot   : Create an output plot of selected fields.

        Returns:
        --------
        fields : A FieldArray of the selected fields.
        """
        # Import the dither function here...
        #def dither(ra,dec,dx,dy):
        #    return ra,dec

        if mode is None or mode.lower() == 'none':
            def dither(ra,dec,dx,dy):
                return ra,dec
            TILINGS = [(0,0),(0,0),(0,0),(0,0)]
        elif mode.lower() == 'smash_dither':
            TILINGS = [(0,0), (1.0,0.0), (-1.0,0.0), (0.0,-0.75)]
            dither = self.smash_dither
        elif mode.lower() == 'smash_rotate':
            TILINGS = [(0,0), (0.75,0.75), (-0.75,0.75), (0.0,-0.75)]
            dither = self.smash_rotate
        elif mode.lower() == 'decam_dither':
            TILINGS = [(0., 0.),(8/3.*CCD_X, -11/3.*CCD_Y),
                       (8/3.*CCD_X, 8/3.*CCD_Y),(-8/3.*CCD_X, 0.)]
            dither = self.decam_dither

        if infile is None:
            infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
        data = np.recfromtxt(infile, names=True)

        # Apply footprint selection after tiling/dither
        #sel = obztak.utils.projector.footprint(data['RA'],data['DEC'])

        # This is currently a non-op
        smash_id = data['ID']
        ra       = data['RA']
        dec      = data['DEC']

        nhexes = len(data)
        #ntilings = len(DECAM_DITHERS)
        ntilings = len(TILINGS)
        nbands = len(BANDS)
        nfields = nhexes*nbands*ntilings

        logging.info("Number of hexes: %d"%nhexes)
        logging.info("Number of tilings: %d"%ntilings)
        logging.info("Number of filters: %d"%nbands)

        fields = FieldArray(nfields)
        fields['HEX'] = np.tile(np.repeat(smash_id,nbands),ntilings)
        fields['PRIORITY'].fill(1)
        fields['TILING'] = np.repeat(np.arange(1,ntilings+1),nhexes*nbands)
        fields['FILTER'] = np.tile(BANDS,nhexes*ntilings)

        #for i in range(ntilings):
        for i,tiling in enumerate(TILINGS):
            idx0 = i*nhexes*nbands
            idx1 = idx0+nhexes*nbands
            ra_dither,dec_dither = dither(ra,dec,tiling[0],tiling[1])
            fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
            fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

        # Apply footprint selection after tiling/dither
        sel = self.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION
        if smcnod:
            # Include SMC northern overdensity fields
            sel_smcnod = self.footprintSMCNOD(fields) # SMCNOD OPERATION
            sel = sel | sel_smcnod
            #sel = sel_smcnod
            fields['PRIORITY'][sel_smcnod] = 99
        #if True:
        #    # Include 'bridge' region between Magellanic Clouds
        #    sel_bridge = self.footprintBridge(fields['RA'],fields['DEC'])
        #    sel = sel | sel_bridge
        sel = sel & (fields['DEC'] > constants.SOUTHERN_REACH)
        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            import obztak.utils.ortho

            plt.ion()

            fig, basemap = obztak.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

            proj = basemap.proj(fields['RA'],fields['DEC'])
            basemap.scatter(*proj, c=fields['TILING'], edgecolor='none', s=50, cmap='viridis_r',vmin=0,vmax=len(TILINGS))
            colorbar = plt.colorbar(label='Tiling')

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'.png'
                fig.savefig(outfig,bbox_inches='tight')
            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile: fields.write(outfile)

        return fields

    def survey_prepare(self, args):
        windows = self.prepare_windows(self.nights, outfile=args.windows, standards=args.standards)
        #infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
        infile = None # Use survey default
        return self.prepare_fields(infile=infile,outfile=args.fields,plot=args.plot,smcnod=args.smcnod)


    @staticmethod
    def footprint(ra,dec):
        """" Dummy footprint selection.

        Parameters:
        -----------
        ra : Right ascension (deg)
        dec: Declination (deg)

        Returns:
        --------
        sel : Selection of fields within the footprint

        """
        sel = np.ones(len(ra),dtype=bool)
        return sel

    @staticmethod
    def smash_dither(ra,dec,dx,dy):
        """Convert to SMASH coordinates, then dither, and convert back.

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Offset in x-dimension/ra (decimal degrees)
        dy : Offset in y-dimension/dec (decimal degrees)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple

        """
        ra0,dec0 = SMASH_POLE
        # Rotate the pole at SMASH_POLE to (0,-90)
        R1 = obztak.utils.projector.SphericalRotator(ra0,90+dec0);
        ra1,dec1 = R1.rotate(ra,dec)
        # Dither the offset
        ra2 = ra1+dx/np.cos(np.radians(dec1))
        dec2 = dec1+dy
        # Rotate back to the original frame (keeping R2)
        return R1.rotate(ra2,dec2,invert=True)

    @staticmethod
    def decam_dither(ra,dec,dx,dy):
        """Conventional dither in DECam coordinates. Consider alternative
        DECam dither in "image" coordinates (see scratch/dither.py).

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Offset in x-dimension/ra (decimal degrees)
        dy : Offset in y-dimension/dec (decimal degrees)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple
        """
        out = []
        for _ra,_dec in zip(ra,dec):
            out.append(obztak.utils.projector.SphericalRotator(_ra,_dec).rotate(dx,dy,invert=True))
        return np.array(out).T

    @staticmethod
    def smash_rotate(ra,dec,dx,dy):
        """Rotate to SMASH coordinates, then rotate by desired offset, and
        convert back.  dx, dy specify offset in decimal degrees.

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Offset in x-dimension/ra (decimal degrees)
        dy : Offset in y-dimension/dec (decimal degrees)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple
        """
        ra0,dec0 = SMASH_POLE
        # Rotate the pole at SMASH_POLE to (0,-90)
        R1 = obztak.utils.projector.SphericalRotator(ra0,90+dec0);
        ra1,dec1 = R1.rotate(ra,dec)
        # Rotate around the SMASH pole
        R2 = obztak.utils.projector.SphericalRotator(dx,dy)
        ra2,dec2 = R2.rotate(ra1,dec1)
        # Rotate back to the original frame (keeping the R2 shift)
        return R1.rotate(ra2,dec2,invert=True)

    @staticmethod
    def rotate(ra,dec,dx,dy):
        """Rotate the celesitial coordinate system by desired offset and
        convert back.  dx, dy specify offset in decimal degrees.

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Offset in x-dimension/ra (decimal degrees)
        dy : Offset in y-dimension/dec (decimal degrees)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple
        """
        R = obztak.utils.projector.SphericalRotator(dx,dy)
        return R.rotate(ra,dec)

    @staticmethod
    def coord_rotate(ra,dec,dx,dy):
        """Rotate the celesitial coordinate system by desired offset and
        convert back.  dx, dy specify offset in decimal degrees.

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Offset in x-dimension/ra (decimal degrees)
        dy : Offset in y-dimension/dec (decimal degrees)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple
        """
        ra0,dec0 = (0,-90)
        # Rotate the pole at SMASH_POLE to (0,-90)
        R1 = obztak.utils.projector.SphericalRotator(ra0,90+dec0);
        ra1,dec1 = R1.rotate(ra,dec)
        # Rotate around the SMASH pole
        R2 = obztak.utils.projector.SphericalRotator(dx,dy)
        ra2,dec2 = R2.rotate(ra1,dec1)
        # Rotate back to the original frame (keeping the R2 shift)
        return R1.rotate(ra2,dec2,invert=True)

    @staticmethod
    def decals_rotate(ra,dec,dx,dy):
        """Perform a Euler angle rotation of the celesitial coordinates and
        return the rotated position of ra,dec in the original
        coordinate system. dx,dy specify the Z and Y Euler rotation
        angles in decimal degrees, respectively.

        Parameters:
        -----------
        ra : Input right ascension
        dec: Input declination
        dx : Rotation in x-dimension/Euler Z/ra (deg)
        dy : Rotation in y-dimension/Euler Y/dec (deg)

        Returns:
        --------
        ra, dec : Dithered ra,dec tuple

        """
        from astropy.modeling.rotations import EulerAngleRotation
        R = EulerAngleRotation(dx,dy,0,'zyx')
        ra1,dec1 = R(ra,dec)
        ra1 += 360 * (ra1 < 0)
        return ra1,dec1


def parser():
    import argparse
    from obztak.utils.parser import Parser, DatetimeAction
    description = __doc__
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = Parser(description=description,formatter_class=formatter)
    parser.add_argument('--survey',default=None,
                        help='Type of survey to schedule.')
    parser.add_argument('-p','--plot',action='store_true',
                        help='Plot output.')
    parser.add_argument('-f','--fields',default='target_fields.csv',
                        help='List of all target fields.')
    parser.add_argument('-w','--windows',default='observation_windows.csv',
                        help='List of observation windows.')
    parser.add_argument('-d','--dither',default='smash_dither',
                        help='Dithering scheme.')
    parser.add_argument('-s','--smcnod',action='store_true',
                        help='Include SMC Northern Overdensity fields.')
    parser.add_argument('--no-standards',action='store_false',dest='standards',
                        help = "Don't include time for standard star observations.")
    return parser

############################################################

if __name__ == '__main__':
    main()

############################################################
