import re
import numpy as np

############################################################

def angsep(lon1,lat1,lon2,lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1,lat1 = np.radians([lon1,lat1])
    lon2,lat2 = np.radians([lon2,lat2])
    
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1,num2), denominator))

############################################################

def airmass(lon_zenith, lat_zenith, lon, lat):
    """
    Safety handling when angular separation to zenith is more than 90 deg
    """
    airmass = 1. / np.cos(np.radians(angsep(lon, lat, lon_zenith, lat_zenith)))
    airmass[airmass < 1.] = 999.
    return airmass

############################################################

def galToCel(ll, bb):
    """
    Converts Galactic (deg) to Celestial J2000 (deg) coordinates
    """
    bb = np.radians(bb)
    sin_bb = np.sin(bb)
    cos_bb = np.cos(bb)

    ll = np.radians(ll)
    ra_gp = np.radians(192.85948)
    de_gp = np.radians(27.12825)
    lcp = np.radians(122.932)

    sin_lcp_ll = np.sin(lcp - ll)
    cos_lcp_ll = np.cos(lcp - ll)

    sin_d = (np.sin(de_gp) * sin_bb) \
            + (np.cos(de_gp) * cos_bb * cos_lcp_ll)
    ramragp = np.arctan2(cos_bb * sin_lcp_ll,
                         (np.cos(de_gp) * sin_bb) \
                         - (np.sin(de_gp) * cos_bb * cos_lcp_ll))
    dec = np.arcsin(sin_d)
    ra = (ramragp + ra_gp + (2. * np.pi)) % (2. * np.pi)
    return np.degrees(ra), np.degrees(dec)

############################################################

def celToGal(ra, dec):
    """
    Converts Celestial J2000 (deg) to Calactic (deg) coordinates
    """
    dec = np.radians(dec)
    sin_dec = np.sin(dec)
    cos_dec = np.cos(dec)

    ra = np.radians(ra)    
    ra_gp = np.radians(192.85948)
    de_gp = np.radians(27.12825)

    sin_ra_gp = np.sin(ra - ra_gp)
    cos_ra_gp = np.cos(ra - ra_gp)

    lcp = np.radians(122.932)    
    sin_b = (np.sin(de_gp) * sin_dec) \
            + (np.cos(de_gp) * cos_dec * cos_ra_gp)
    lcpml = np.arctan2(cos_dec * sin_ra_gp,
                       (np.cos(de_gp) * sin_dec) \
                       - (np.sin(de_gp) * cos_dec * cos_ra_gp))
    bb = np.arcsin(sin_b)
    ll = (lcp - lcpml + (2. * np.pi)) % (2. * np.pi)
    return np.degrees(ll), np.degrees(bb)

############################################################

def hms2dec(hms):
    """
    Convert longitude from hours,minutes,seconds in string or 3-array
    format to decimal degrees.
    """
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.

    if isinstance(hms,basestring):
        hour,minute,second = np.array(hms.split(':')).astype(float)
        #hour,minute,second = np.array(re.split('[hms]',hms))[:3].astype(float)
    else:
        hour,minute,second = hms.T

    decimal = (hour + minute * 1./MINUTE + second * 1./SECOND)*(DEGREE/HOUR)
    return decimal

############################################################
