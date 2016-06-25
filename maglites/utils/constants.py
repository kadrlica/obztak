"""
Constants.
"""
import ephem
from collections import OrderedDict as odict
import numpy as np

# Plotting DECam
DECAM=1.1 # DECam radius (deg)
# Marker size depends on figsize and DPI
FIGSIZE=(10.5,8.5)
SCALE=np.sqrt((8.0*6.0)/(FIGSIZE[0]*FIGSIZE[1]))
DPI=80;

# LMC and SMC
RA_LMC = 80.8939
DEC_LMC = -69.7561
RADIUS_LMC = 5.3667 # semi-major axis (deg)
RA_SMC = 13.1867
DEC_SMC = -72.8286
RADIUS_SMC = 2.667 # semi_major axis (deg)

# SMC Northern Oversdensity
RA_SMCNOD = 12.0
DEC_SMCNOD = -64.8
HEX_SMCNOD = [1651, 1575, 1576, 1804] # baseline
#HEX_SMCNOD = [1575, 1576, 1804, 1960] # alternative
TILING_SMCNOD = [1, 3]

# http://www.ctio.noao.edu/noao/content/coordinates-observatories-cerro-tololo-and-cerro-pachon
LON_CTIO = '-70:48:23.49'
LAT_CTIO = '-30:10:10.78'
ELEVATION_CTIO = 2206.8 # m

# Pole of the SMASH fields (RA,DEC)
SMASH_POLE = (10., -30.)

# Characteristics of the survey
# 90 sec exposures with 30 sec between exposures
EXPTIME   = 90*ephem.second # Exposure time
DOWNTIME  = 30*ephem.second # Time between exposures from readout/slew
NEXP      = 2 # Number of exposures taken in a row
FIELDTIME = EXPTIME+DOWNTIME
BANDS     = ('g','r')

# Time for taking standards
STANDARDS = 10*ephem.minute

# Characteristics of DECam
ARCSEC_TO_DEGREE = 1. / (60. * 60.)
PIXEL_SCALE = 0.2626 * ARCSEC_TO_DEGREE
NPIX_X = 4096
NPIX_Y = 2048
CCD_X = NPIX_X * PIXEL_SCALE # degree
CCD_Y = NPIX_Y * PIXEL_SCALE # degree

# Blanco characteritics
SOUTHERN_REACH = -89.

# SISPI json template formatting
PROPID = '2016A-0366'
OBJECT_FMT = "MAGLITES field - %(ID)d.%(TILING)d.%(PRIORITY)d"
SEQID_FMT = "MAGLITES scheduled - %(DATE)s"
FLOAT_FMT = '%.4f'
SISPI_DICT = odict([
    ("object",  None),
    ("seqnum",  None), # 1-indexed
    ("seqtot",  2),
    ("seqid",   None),
    ("expTime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
    ("propid",  PROPID),
    ("comment", ""),
])

def FIELD2OBJECT(field):
    return OBJECT_FMT%(field)

def OBJECT2FIELD(object_str):
    ID,TILING,PRIORITY = map(int,object_str.split(' - ')[-1].split('.'))
    return dict(ID=ID,TILING=TILING,PRIORITY=PRIORITY)

def FIELD2SEQID(field):
    return SEQID_FMT%(field)

def SEQID2FIELD(seqid_str):
    DATE = str(seqid_str.split(' - ')[-1].strip())
    return dict(DATE=DATE)

