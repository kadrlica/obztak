"""
Observation planning for the Magellanic Satellites Survey.
"""

try:
    from .database import Database
except ImportError:
    pass


############################################################

def datestring(date): 
    import ephem
    date = ephem.Date(date)
    datetuple = date.tuple()
    seconds = round(datetuple[-1],4)
    minutes = datetuple[-2]
    minutes += seconds//60
    seconds = seconds%60.
    return '%s/%s/%s %02i:%02i:%07.4f'%(datetuple[:-2]+(minutes,seconds))

############################################################


