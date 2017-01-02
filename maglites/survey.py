"""
Decide which fields to observe and time windows to observe.
"""

import os,sys
import logging
import numpy.lib.recfunctions as recfunc
import numpy as np
import ephem

import maglites.utils.projector
import maglites.utils.constants
import maglites.utils.ortho

from maglites.utils.projector import cel2gal
from maglites.utils.ortho import datestring
from maglites.field import FieldArray
from maglites.utils.constants import BANDS,SMASH_POLE,CCD_X,CCD_Y,STANDARDS
from maglites.utils import fileio

class Survey(object):
    """Base class for preparing a survey. Creates a list of observation
    windows (observing dates and times) and target fields. The
    `prepare_windows` and `prepare_fields` methods should be
    subclassed for specific surveys.
    """

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

    def prepare_fields(self,polygon, mode=None):
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

        if not infile:
            infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
        data = np.recfromtxt(infile, names=True)

        nhexes = len(data)
        #ntilings = len(DECAM_DITHERS)
        ntilings = len(TILINGS)
        nbands = len(BANDS)
        nfields = nhexes*nbands*ntilings

        logging.info("Number of hexes: %d"%nhexes)
        logging.info("Number of tilings: %d"%ntilings)
        logging.info("Number of filters: %d"%nbands)

        #for i in range(ntilings):
        for i,tiling in enumerate(TILINGS):
            idx0 = i*nhexes*nbands
            idx1 = idx0+nhexes*nbands
            ra_dither,dec_dither = dither(ra,dec,tiling[0],tiling[1])
            fields['RA'][idx0:idx1] = np.repeat(ra_dither,nbands)
            fields['DEC'][idx0:idx1] = np.repeat(dec_dither,nbands)

        # Apply footprint selection after tiling/dither
        sel = maglites.utils.projector.footprint(fields['RA'],fields['DEC']) # NORMAL OPERATION

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
        R1 = maglites.utils.projector.SphericalRotator(ra0,90+dec0);
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
            out.append(maglites.utils.projector.SphericalRotator(_ra,_dec).rotate(dx,dy,invert=True))
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
        R1 = maglites.utils.projector.SphericalRotator(ra0,90+dec0);
        ra1,dec1 = R1.rotate(ra,dec)
        # Rotate around the SMASH pole
        R2 = maglites.utils.projector.SphericalRotator(dx,dy)
        ra2,dec2 = R2.rotate(ra1,dec1)
        # Rotate back to the original frame (keeping the R2 shift)
        return R1.rotate(ra2,dec2,invert=True)


class MagLiteS(Survey):
    """ Survey sublcass for MagLiteS. """

    """
    # One season prediction
    nights = [['2016/2/10', 'second'],
              ['2016/2/11', 'second'],
              ['2016/2/12', 'second'],
              ['2016/2/13', 'second'],
              ['2016/2/14', 'second'],
              ['2016/2/15', 'second'],
              ['2016/6/27', 'full'],
              ['2016/6/28', 'full'],
              ['2016/6/29', 'full']]
    """
    """
    # Two seasons prediction
    nights = [['2016/2/10', 'second'],
              ['2016/2/11', 'second'],
              ['2016/2/12', 'second'],
              ['2016/2/13', 'second'],
              ['2016/2/14', 'second'],
              ['2016/2/15', 'second'],
              ['2016/6/27', 'full'],
              ['2016/6/28', 'full'],
              ['2016/6/29', 'full'],
              ['2017/2/18', 'second'],
              ['2017/2/19', 'second'],
              ['2017/2/20', 'second'],
              ['2017/2/21', 'second'],
              ['2017/2/22', 'second'],
              ['2017/2/23', 'second'],
              ['2017/6/27', 'full'],
              ['2017/6/28', 'full'],
              ['2017/6/29', 'full']]
    """

    """
    In 2017A, we are considering a few different possible observing periods.

    Early run (choose one)
    * 2017/2/01
    * 2017/3/02

    + Late run (choose one)
    * 2017/6/17
    * 2017/7/16
    """
    """
    # 2017A prediction (feb1jun17)
    nights = [['2017/2/01', 'second'],
              ['2017/2/02', 'second'],
              ['2017/2/03', 'second'],
              ['2017/2/04', 'second'],
              ['2017/2/05', 'second'],
              ['2017/2/06', 'second'],
              ['2017/6/17', 'full'],
              ['2017/6/18', 'full'],
              ['2017/6/19', 'full']]

    # 2017A prediction (feb1jul16)
    nights = [['2017/2/01', 'second'],
              ['2017/2/02', 'second'],
              ['2017/2/03', 'second'],
              ['2017/2/04', 'second'],
              ['2017/2/05', 'second'],
              ['2017/2/06', 'second'],
              ['2017/7/16', 'full'],
              ['2017/7/17', 'full'],
              ['2017/7/18', 'full']]

    # 2017A prediction (mar2jun17)
    nights = [['2017/3/02', 'second'],
              ['2017/3/03', 'second'],
              ['2017/3/04', 'second'],
              ['2017/3/05', 'second'],
              ['2017/3/06', 'second'],
              ['2017/3/07', 'second'],
              ['2017/6/17', 'full'],
              ['2017/6/18', 'full'],
              ['2017/6/19', 'full']]
    """
    """
    # 2017A prediction (mar2jul16)
    nights = [['2017/3/02', 'second'],
              ['2017/3/03', 'second'],
              ['2017/3/04', 'second'],
              ['2017/3/05', 'second'],
              ['2017/3/06', 'second'],
              ['2017/3/07', 'second'],
              ['2017/7/16', 'full'],
              ['2017/7/17', 'full'],
              ['2017/7/18', 'full']]
    """

    # 2017A prediction (Moon up during second half of night)
    #nights = [['2017/2/18', 'second'],
    #          ['2017/2/19', 'second'],
    #          ['2017/2/20', 'second'],
    #          ['2017/2/21', 'second'],
    #          ['2017/2/22', 'second'],
    #          ['2017/2/23', 'second'],
    #          ['2017/6/17', 'full'],
    #          ['2017/6/18', 'full'],
    #          ['2017/6/19', 'full']]

    # 2017A ACTUAL
    nights = [['2017/2/21', 'full'],
              ['2017/2/22', 'full'],
              ['2017/2/23', 'full'],
              ['2017/6/18', 'full'],
              ['2017/6/19', 'full'],
              ['2017/6/20', 'full']]

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
        #sel = maglites.utils.projector.footprint(data['RA'],data['DEC'])

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
        sel = sel & (fields['DEC'] > maglites.utils.constants.SOUTHERN_REACH)
        fields = fields[sel]

        logging.info("Number of target fields: %d"%len(fields))

        if plot:
            import pylab as plt
            plt.ion()

            fig, basemap = maglites.utils.ortho.makePlot('2016/2/11 03:00',center=(0,-90),airmass=False,moon=False)

            proj = maglites.utils.ortho.safeProj(basemap,fields['RA'],fields['DEC'])
            basemap.scatter(*proj, c=fields['TILING'], edgecolor='none', s=50, cmap='Spectral',vmin=0,vmax=len(TILINGS))
            colorbar = plt.colorbar(label='Tiling')

            if outfile:
                outfig = os.path.splitext(outfile)[0]+'.png'
                fig.savefig(outfig,bbox_inches='tight')
            if not sys.flags.interactive:
                plt.show(block=True)

        if outfile: fields.write(outfile)

        return fields

    @staticmethod
    def footprint(ra,dec):
        l, b = cel2gal(ra, dec)

        angsep_lmc = maglites.utils.projector.angsep(maglites.utils.constants.RA_LMC, maglites.utils.constants.DEC_LMC, ra, dec)
        angsep_smc = maglites.utils.projector.angsep(maglites.utils.constants.RA_SMC, maglites.utils.constants.DEC_SMC, ra, dec)
        sel = (np.fabs(b) > 10.) \
              & ((angsep_lmc < 30.) | (angsep_smc < 30.)) \
              & (dec < -55.) & (ra > 100.) & (ra < 300.)
        #sel = sel | ((dec < -65.) & (angsep_lmc > 5.) & (angsep_smc > 5.))
        sel = sel | ((dec < -65.) & (ra > 300.) & (ra < 360.)) # SMC
        sel = sel | (dec < -80.)

        return sel

    @staticmethod
    def footprintSMCNOD(fields):
        """
        Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
        """
        sel = np.in1d(fields['HEX'], maglites.utils.constants.HEX_SMCNOD) \
              & np.in1d(fields['TILING'], maglites.utils.constants.TILING_SMCNOD)
        return sel

    @staticmethod
    def footprintBridge(ra, dec):
        """
        Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
        """
        sel = (ra > 30.) & (ra < 60.) & (dec < -65.)
        return sel

def parser():
    import argparse
    from maglites.utils.parser import Parser, DatetimeAction
    description = __doc__
    formatter = argparse.ArgumentDefaultsHelpFormatter
    parser = Parser(description=description,formatter_class=formatter)
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
    parser.add_argument('--no-standards',action='store_false',dest='standards', help = "Don't include time for standard star observations.")
    return parser

def main():
    args = parser().parse_args()
    survey = MagLiteS()
    windows = survey.prepare_windows(survey.nights, outfile=args.windows, standards=args.standards)
    infile = os.path.join(fileio.get_datadir(),'smash_fields_alltiles.txt')
    survey.prepare_fields(infile=infile,outfile=args.fields,plot=args.plot,smcnod=args.smcnod)


############################################################

if __name__ == '__main__':
    main()

############################################################
