"""
Constants.
"""
import ephem
from collections import OrderedDict as odict

RA_LMC = 80.8939
DEC_LMC = -69.7561
RA_SMC = 13.1867
DEC_SMC = -72.8286

# http://www.ctio.noao.edu/noao/content/coordinates-observatories-cerro-tololo-and-cerro-pachon

LON_CTIO = '-70:48:23.49'
LAT_CTIO = '-30:10:10.78'
ELEVATION_CTIO = 2206.8 # m

# Characteristics of the survey
# 4 minutes = 90 sec exposures in g and r with 30 sec between exposures
EXPTIME   = 90*ephem.second # Exposure time
DOWNTIME  = 30*ephem.second # Time between exposures from readout/slew
NEXP      = 2 # Number of exposures taken in a row
FIELDTIME = NEXP*(EXPTIME+DOWNTIME)
BANDS     = ('g','r')

# SISPI json template formatting
OBJECT_FMT = "MAGLITES field - %(ID)d.%(TILING)d.%(PRIORITY)d"
SEQID_FMT = "MAGLITES scheduled - %(DATE)s"
FLOAT_FMT = '%.4f'
SISPI_DICT = odict([
    ("seqtot",  2),
    ("seqnum",  None), # 1-indexed
    ("seqid",   None),
    ("object",  None),
    ("exptime", 90),
    ("RA",      None),
    ("dec",     None),
    ("filter",  None),
    ("count",   1),
    ("expType", "object"),
    ("program", "maglites"),
    ("wait",    "False"),
])
