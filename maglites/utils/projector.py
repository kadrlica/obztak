import re
import numpy as np

import maglites.utils.constants
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

gal2cel = galToCel

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

cel2gal = celToGal

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

class SphericalRotator:
    """
    Base class for rotating points on a sphere.

    The input is a fiducial point (deg) which becomes (0, 0) in rotated coordinates.
    """

    def __init__(self, lon_ref, lat_ref, zenithal=False):
        self.setReference(lon_ref, lat_ref, zenithal)

    def setReference(self, lon_ref, lat_ref, zenithal=False):

        if zenithal:
            phi = (np.pi / 2.) + np.radians(lon_ref)
            theta = (np.pi / 2.) - np.radians(lat_ref)
            psi = 0.
        if not zenithal:
            phi = (-np.pi / 2.) + np.radians(lon_ref)
            theta = np.radians(lat_ref)
            psi = np.radians(90.) # psi = 90 corresponds to (0, 0) psi = -90 corresponds to (180, 0)
        

        cos_psi,sin_psi = np.cos(psi),np.sin(psi)
        cos_phi,sin_phi = np.cos(phi),np.sin(phi)
        cos_theta,sin_theta = np.cos(theta),np.sin(theta)

        self.rotation_matrix = np.matrix([
            [cos_psi * cos_phi - cos_theta * sin_phi * sin_psi,
             cos_psi * sin_phi + cos_theta * cos_phi * sin_psi,
             sin_psi * sin_theta],
            [-sin_psi * cos_phi - cos_theta * sin_phi * cos_psi,
             -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi,
             cos_psi * sin_theta],
            [sin_theta * sin_phi,
             -sin_theta * cos_phi,
             cos_theta]
        ])
        
        self.inverted_rotation_matrix = np.linalg.inv(self.rotation_matrix)

    def cartesian(self,lon,lat):
        lon = np.radians(lon)
        lat = np.radians(lat) 
        
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z =  np.sin(lat)
        return np.array([x,y,z])
        

    def rotate(self, lon, lat, invert=False):
        vec = self.cartesian(lon,lat)

        if invert:
            vec_prime = np.dot(np.array(self.inverted_rotation_matrix), vec)
        else:        
            vec_prime = np.dot(np.array(self.rotation_matrix), vec)

        lon_prime = np.arctan2(vec_prime[1], vec_prime[0])
        lat_prime = np.arcsin(vec_prime[2])

        return (np.degrees(lon_prime) % 360.), np.degrees(lat_prime)

############################################################

class Projector:
    """
    Class for performing two-dimensional map projections from the celestial sphere.
    """

    def __init__(self, lon_ref, lat_ref, proj_type = 'ait'):
        self.lon_ref = lon_ref
        self.lat_ref = lat_ref
        self.proj_type = proj_type

        if proj_type.lower() == 'ait':
            self.rotator = SphericalRotator(lon_ref, lat_ref, zenithal=False)
            self.sphere_to_image_func = aitoffSphereToImage
            self.image_to_sphere_func = aitoffImageToSphere
        elif proj_type.lower() == 'tan':
            self.rotator = SphericalRotator(lon_ref, lat_ref, zenithal=True)
            self.sphere_to_image_func = gnomonicSphereToImage
            self.image_to_sphere_func = gnomonicImageToSphere
        elif proj_type.lower() == 'car':
            def rotate(lon,lat,invert=False):
                if invert:
                    return lon + np.array([lon_ref]), lat + np.array([lat_ref])
                else:
                    return lon - np.array([lon_ref]), lat - np.array([lat_ref])
            self.rotator = SphericalRotator(lon_ref, lat_ref, zenithal=False)
            # Monkey patch the rotate function
            self.rotator.rotate = rotate 
            self.sphere_to_image_func = cartesianSphereToImage
            self.image_to_sphere_func = cartesianImageToSphere
        else:
            print 'WARNING: %s not recognized'%(proj_type)

    def sphereToImage(self, lon, lat):
        lon_rotated, lat_rotated = self.rotator.rotate(lon, lat)
        return self.sphere_to_image_func(lon_rotated, lat_rotated)

    sphere2image = sphereToImage
        
    def imageToSphere(self, x, y):
        lon_rotated, lat_rotated = self.image_to_sphere_func(x, y)
        return self.rotator.rotate(lon_rotated, lat_rotated, invert = True)

    image2sphere = imageToSphere

def aitoffSphereToImage(lon, lat):
    """
    Hammer-Aitoff projection (deg).
    """
    lon = lon - 360.*(lon>180)
    lon = np.radians(lon)
    lat = np.radians(lat)

    half_lon = lon/2.
    cos_lat = np.cos(lat)
     
    gamma = (180. / np.pi) * np.sqrt(2. / (1. + (cos_lat * np.cos(half_lon))))
    x = 2. * gamma * cos_lat * np.sin(half_lon)
    y = gamma * np.sin(lat)
    return x, y

def aitoffImageToSphere(x, y):
    """
    Inverse Hammer-Aitoff projection (deg).
    """
    x = x - 360.*(x>180)
    x = np.asarray(np.radians(x))
    y = np.asarray(np.radians(y))
    z = np.sqrt(1. - (x / 4.)**2 - (y / 2.)**2) # rad
    lon = 2. * np.arctan2((2. * z**2) - 1, (z / 2.) * x)
    lat = np.arcsin( y * z)
    return ((180. - np.degrees(lon)) % 360.), np.degrees(lat)

############################################################

def match(lon1, lat1, lon2, lat2, tol=None, nnearest=1):
    """
    Adapted from Eric Tollerud.
    Finds matches in one catalog to another.
 
    Parameters
    lon1 : array-like
        Longitude of the first catalog (degrees)
    lat1 : array-like
        Latitude of the first catalog (shape of array must match `lon1`)
    lon2 : array-like
        Longitude of the second catalog
    lat2 : array-like
        Latitude of the second catalog (shape of array must match `lon2`)
    tol : float or None, optional
        Proximity (degrees) of a match to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned.
    nnearest : int, optional
        The nth neighbor to find.  E.g., 1 for the nearest nearby, 2 for the
        second nearest neighbor, etc.  Particularly useful if you want to get
        the nearest *non-self* neighbor of a catalog.  To do this, use:
        ``spherematch(lon, lat, lon, lat, nnearest=2)``
 
    Returns
    -------
    idx1 : int array
        Indices into the first catalog of the matches. Will never be
        larger than `lon1`/`lat1`.
    idx2 : int array
        Indices into the second catalog of the matches. Will never be
        larger than `lon2`/`lat2`.
    ds : float array
        Distance (in degrees) between the matches
    """
    from scipy.spatial import cKDTree
 
    lon1 = np.asarray(lon1)
    lat1 = np.asarray(lat1)
    lon2 = np.asarray(lon2)
    lat2 = np.asarray(lat2)
 
    if lon1.shape != lat1.shape:
        raise ValueError('lon1 and lat1 do not match!')
    if lon2.shape != lat2.shape:
        raise ValueError('lon2 and lat2 do not match!')

    rotator = SphericalRotator(0,0)

 
    # This is equivalent, but faster than just doing np.array([x1, y1, z1]).T
    x1, y1, z1 = rotator.cartesian(lon1.ravel(),lat1.ravel())
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1
 
    x2, y2, z2 = rotator.cartesian(lon2.ravel(),lat2.ravel())
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0] = x2
    coords2[:, 1] = y2
    coords2[:, 2] = z2
 
    tree = cKDTree(coords2)
    if nnearest == 1:
        idxs2 = tree.query(coords1)[1]
    elif nnearest > 1:
        idxs2 = tree.query(coords1, nnearest)[1][:, -1]
    else:
        raise ValueError('invalid nnearest ' + str(nnearest))
 
    ds = angsep(lon1, lat1, lon2[idxs2], lat2[idxs2])
 
    idxs1 = np.arange(lon1.size)
 
    if tol is not None:
        msk = ds < tol
        idxs1 = idxs1[msk]
        idxs2 = idxs2[msk]
        ds = ds[msk]
 
    return idxs1, idxs2, ds

def footprint(ra,dec):
    l, b = celToGal(ra, dec)

    angsep_lmc = maglites.utils.projector.angsep(maglites.utils.constants.RA_LMC, maglites.utils.constants.DEC_LMC, ra, dec)
    angsep_smc = maglites.utils.projector.angsep(maglites.utils.constants.RA_SMC, maglites.utils.constants.DEC_SMC, ra, dec)
    cut = (np.fabs(b) > 10.) \
          & ((angsep_lmc < 30.) | (angsep_smc < 30.)) \
          & (dec < -55.) & (ra > 100.) & (ra < 300.)
    #cut = cut | ((dec < -65.) & (angsep_lmc > 5.) & (angsep_smc > 5.))
    cut = cut | ((dec < -65.) & (ra > 300.) & (ra < 360.)) # SMC
    cut = cut | (dec < -80.)

    return cut

def footprintSMCNOD(fields):
    """
    Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
    """
    cut = np.in1d(fields['HEX'], maglites.utils.constants.HEX_SMCNOD) \
          & np.in1d(fields['TILING'], maglites.utils.constants.TILING_SMCNOD)
    return cut

def footprintBridge(ra, dec):
    """
    Special selection for pointings near the SMC Northern Overdensity (SMCNOD)
    """
    cut = (ra > 30.) & (ra < 60.) & (dec < -65.) 
    return cut

