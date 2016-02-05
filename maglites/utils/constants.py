"""
Constants.
"""
import ephem

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

# Object formatting
OBJECT_FMT = "MAGLITES field - %d.%d.%d"

