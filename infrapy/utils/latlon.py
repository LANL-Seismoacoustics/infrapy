# latlon.py
#
# Lightweight NumPy-dependent python library of
# functions for latitude & longitude data.
#
# Function Usage Glossary:
#  d=azdiff(az1,az2)
#  tf=azinrng(azrng,az)
#  avg=azmean(az,dim)
#  latlon=fixlatlon(latlon)
#  latloni0,latloni1=gc_intersect(latlona0,latlona1,latlonb0,latlonb1)
#  latloni=gcarc_intersect(latlona0,latlona1,latlonb0,latlonb1)
#  lat=geocentric2geographiclat(lat,ecc)
#  xyz=geocentric2xyz(latlon,radius)
#  lat=geographic2geocentriclat(lat,ecc)
#  radius=geographiclat2radius(lat,ellipsoid)
#  gcdist=haversine(latlon0,latlon1)
#  tf=inlatlonbox(latrng,lonrng,latlon)
#  tf=inlonrng(rng,lon)
#  wlat,px=latmod(lat,wrap)
#  wlon=lonmod(lon,wrap)
#  latlon=randlatlon(n)
#  x=randsphere(n,dim)
#  latlon1,baz=sphericalfwd(latlon0,gcdist,az)
#  gcdist,az,baz=sphericalinv(latlon0,latlon1)
#  latlon,radius=xyz2geocentric(xyz)
#
# To test:
#  python latlon.py -v
#
# Author            Garrett Euler (ggeuler at lanl dot gov)
# Created           Feb. 22, 2016
# Last Modified     Mar.  1, 2016

import numpy as np
from numpy.core.umath_tests import inner1d

# functions to port immediately:
# sph_poly_in
# sph_poly_intersect
# vincentyfwd
# vincentyinv
# xyz2geographic
# geographic2xyz
# enu2geographic
# geographic2enu

# functions to port eventually:
# arrayaperture
# arraycenter
# arrayradius
# gridconv
# geographic2geocentric
# geocentric2geographic
# closest_point_on_gc
# degdist_from_gc
# gc2latlon
# gcarc2latlon
# mean_ellipsoid_radius
# gcarc_count

# todo:
# - more functions
# - more doctests
# - docstring stuff
# - should i do input validation?
# - pole handling for inlatlonbox
# - edge_is_in option for *in*

def azdiff(az1,az2):
    """
    AZDIFF    Returns the angle between azimuths

        Usage:    d=azdiff(az1,az2)

        Description:
         D=AZDIFF(AZ1,AZ2) calculates the angle D from AZ1 to AZ2.  D is
         always within +/-180.

        Notes:

        Examples:
         # 340 --> 20 should give 40:
         >>> azdiff(340,20)
         array([[ 40.]])
         
         # 20 --> 340 should give -40:
         >>> azdiff(20,340)
         array([[-40.]])
         
         # Several examples with equivalent azimuths:
         >>> azdiff([180,0,540],[-180,720,-180])
         array([[ 0.,  0.,  0.]])

        See also: AZMEAN, AZINRNG, INLONRNG, LONMOD, INLATLONBOX, LATMOD
    """

    # convert to numpy array of floats
    try:
        az1=np.array(az1, dtype=np.float)
        az2=np.array(az2, dtype=np.float)
    except (TypeError, ValueError):
        msg="AZ1 & AZ2 must be real-valued inputs!"
        raise TypeError(msg)

    # want angle d between 2 azimuths
    try:
        d=np.atleast_2d(np.mod(az2,360)-np.mod(az1,360))
    except ValueError:
        msg="AZ1 & AZ2 must be equal shaped or scalar!"
        raise ValueError(msg)
    d[d<-180]+=360
    d[d>180]-=360

    return d


def azinrng(azrng,az):
    """
    AZINRNG    Returns TRUE for azimuths within azimuth range

        Usage:    tf=azinrng(azrng,az)

        Description:
         TF=AZINRNG(AZRNG,AZ) returns a logical matrix TF of equal shape to AZ
         with each element set to TRUE for the elements of AZ that are within
         the range AZRNG as [MIN MAX] (e.g., including MIN & MAX).  Note that
         wrapping is done on all azimuths.  So if the range is 10 to 20 deg then
         values within that range and 370 to 380, -350 to -340, etc all are
         within the range.  Also note that the range always extends clockwise
         from AZMIN to AZMAX.  To get the full azimuth range set the limits such
         that AZMAX-AZMIN equals 360.  If you set AZMIN=0 and AZMAX=361 then you
         will actually have a range from 0 to 1.
         
        Notes:
         - This is just a wrapper for INLONRNG.

        Examples:
         % "Dateline" check:
         >>> azinrng([178,-178],np.arange(177,184))
         array([False,  True,  True,  True,  True,  True, False], dtype=bool)

        See also: INLONRNG, AZDIFF, AZMEAN, LONMOD, LATMOD, FIXLATLON,
                  INLATLONBOX
    """

    return inlonrng(azrng,az)


def azmean(az, dim=None):
    """
    AZMEAN    Returns the mean azimuth of a set of azimuths

        Usage:    avg=azmean(az)
                  avg=azmean(az,dim)

        Description:
         AVG=AZMEAN(AZ) returns the average azimuth of the azimuths in AZ.  AZ
         is expected in degrees.  Operates down the first non-singleton
         dimension.

         AVG=AZMEAN(AZ,DIM) takes the mean across along dimension DIM.

        Notes:
         - NaNs are allowed and are ignored.

        Examples:
         # Average North scattered azimuths:
         >>> azmean([-5,-15,5])
         array([-5.])
         
         # Test dimension input:
         >>> azmean([-5,-15,5],0)
         array([ -5., -15.,   5.])
         >>> azmean([-5,-15,5],1)
         array([-5.])
         >>> azmean([[-5,-15,5],[20,40,0]],0)
         array([  7.5,  12.5,   2.5])
         >>> azmean([[-5,-15,5],[20,40,0]],1)
         array([ -5.,  20.])

        See also: AZDIFF, AZINRNG, ARRAYCENTER
    """

    # convert to 2D numpy array of floats
    try:
        az=np.atleast_2d(np.array(az, dtype=np.float))
    except (TypeError, ValueError):
        msg="AZ must be a real-valued input!"
        raise TypeError(msg)

    # default dimension
    if dim==None:
        # find non-singleton dimensions
        nsdim=np.where(np.array(az.shape) != 1)[0]
        if np.size(nsdim)!=0:
            # use first non-singleton dimension
            dim=nsdim[0]
    else:
        try:
            # make sure input is integer
            dim=np.array(dim, dtype=np.int)
        except (TypeError, ValueError):
            msg="DIM must be an integer!"
            raise TypeError(msg)
        
        # require dim to be an appropriately ranged scalar
        if np.size(dim)!=1 or dim>az.ndim:
            msg="DIM must be a scalar integer <= AZ.ndim"
            raise ValueError(msg)

    # azimuths to unit vector, mean, unit vector to azimuth
    avg=180/np.pi*np.arctan2(np.nanmean(np.sin(np.deg2rad(az)), axis=dim),
                             np.nanmean(np.cos(np.deg2rad(az)), axis=dim));

    return avg


#def enu2geographic():


def fixlatlon(latlon):
    """
    FIXLATLON    Returns latitudes & longitudes in reasonable ranges

        Usage:    latlon=fixlatlon(latlon)

        Description:
         LATLON=FIXLATLON(LATLON) returns latitudes within the range +/-90
         and longitudes within +/-180 taking care to preserve the actual
         corresponding location.  The input LATLON must be Nx2 real array.

        Notes:

        Examples:
         # Some dumb programs may go "over the pole" in terms of latitude.
         # FIXLATLON can fix this while handling the accompanying shift in
         # longitude:
         >>> fixlatlon([100, 0])
         array([[  80.,  180.]])
         
         # Bring some longitudes from 0-360 to -180-180:
         >>> fixlatlon([[80, 210],[-45,10],[-5,315]])
         array([[  80., -150.],
                [ -45.,   10.],
                [  -5.,  -45.]])

        See also: LATMOD, LONMOD, INLONRNG, INLATLONBOX
    """

    # check that latlon is numeric & real
    # check that latlon is Nx2
    
    # convert to 2D numpy array of floats
    latlon=np.atleast_2d(np.array(latlon, dtype=np.float))

    # wrap position
    latlon[:,[0]],px=latmod(latlon[:,[0]]);
    latlon[:,[1]]=lonmod(latlon[:,[1]]+np.mod(px,2)*180);

    return latlon


def gc_intersect(latlona0,latlona1,latlonb0,latlonb1):
    """
    GC_INTERSECT    Return intersection points between great circles

        Usage:    latloni0,latloni1=gc_intersect(latlona0,latlona1,
                                                 latlonb0,latlonb1)

        Description:
         LATLONI0,LATLONI1=GC_INTERSECT(LATLONA0,LATLONA1,LATLONB0,LATLONB1)
         finds the intersection points of great circles given by points
         LATLONA0/1 with great circles given by points LATLONB0/1.  Great
         circles either intersect twice or are equal.  LATLONI0 has 1
         intersection point and LATLONI1 gives the other (antipodal to the
         first).  When two great circles are equal both intersection points are
         set to NaNs. All LATLON must either be scalar points (1x2) or arrays
         (Nx2) with the same shape.  This allows finding intersections between
         one great circle and several others or to find intersections between
         distinct pairs.  All inputs must be in degrees! Outputs are in degrees.

        Notes:
         - Assumes positions are given in geocentric coordinates.

        Examples:
         #

        See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
                  GCARC2LATLON, GCARC_INTERSECT
    """
    
    # check that latlona0/1 & latlonb0/1 are (N,1)x2 in shape and real
    
    # convert to 2D numpy arrays of floats
    latlona0=np.atleast_2d(np.array(latlona0, dtype=np.float))
    latlona1=np.atleast_2d(np.array(latlona1, dtype=np.float))
    latlonb0=np.atleast_2d(np.array(latlonb0, dtype=np.float))
    latlonb1=np.atleast_2d(np.array(latlonb1, dtype=np.float))
    
    # convert to xyz
    latlona0=geocentric2xyz(latlona0)
    latlona1=geocentric2xyz(latlona1)
    latlonb0=geocentric2xyz(latlonb0)
    latlonb1=geocentric2xyz(latlonb1)

    # 2. get perpendicular M to plane formed by gc A
    M=np.cross(latlona0,latlona1, axis=1)

    # 3. get perpendicular N to plane formed by gc B
    N=np.cross(latlonb0,latlonb1, axis=1)

    # 4. get perpendicular E to plane formed by gc MN
    # - note that this is fully broadcasted from a0/1 & b0/1
    E=np.cross(M,N, axis=1)

    # 5. get E intersections with sphere
    latloni0=xyz2geocentric(E)
    latloni1=xyz2geocentric(-E)

    # 6. nans if E is basically zero
    # - this occurs when A & B are equivalent
    bad=np.linalg.norm(E, axis=1)
    bad=bad<10.*np.sqrt(np.finfo(np.float).eps)
    if any(bad):
        latloni0[bad,:]=np.nan
        latloni1[bad,:]=np.nan

    return latloni0, latloni1


def gcarc_intersect(latlona0,latlona1,latlonb0,latlonb1):
    """
    GCARC_INTERSECT    Return intersection points between great circle arcs
        
        Usage:    latloni=gcarc_intersect(latlona0,latlona1,latlonb0,latlonb1)
        
        Description:
         LATLONI=GCARC_INTERSECT(LATLONA0,LATLONA1,LATLONB0,LATLONB1) finds the
         intersection points of great circle arcs given by points LATLONA0/1
         with great circle arcs given by points LATLONB0/1.  Arcs are the
         shortest path between points (no arcs are >180deg in length)!  Great
         circle arcs either intersect once or not at all (equal arcs are always
         treated as non-intersecting).  LATLONI gives the intersection points.
         When two great circle arcs do not intersect (or are equal) the
         corresponding intersection is set to NaNs.  All LATLON must either be
         scalar points (1x2) or arrays (Nx2) with the same shape.  This allows
         finding intersections between one great circle arc and several others
         or to find intersections between distinct pairs.  All inputs must be in
         degrees! Outputs are in degrees.
        
        Notes:
         - Arcs are always <=180deg!
         - Assumes positions are given in geocentric coordinates.
        
        Examples:
         # Simple case:
         >>> gcarc_intersect([0,0],[0,10],[-5,5],[5,5])
         array([[-0.,  5.]])
         
         # Do 2 great circle arcs on the same great circle intersect?
         >>> gcarc_intersect([0,0],[0,10],[0,5],[0,15])
         array([[ nan,  nan]])
         
         # How about endpoints?
         >>> gcarc_intersect([0,0],[0,10],[0,0],[10,0])
         array([[ nan,  nan]])
         >>> gcarc_intersect([0,0],[0,10],[-10,0],[10,0])
         array([[ nan,  nan]])
         
         # You can also do several intersections:
         >>> gcarc_intersect([0,0],[0,10],[[0,5],[-5,5]],[[0,15],[5,5]])
         array([[ nan,  nan],
                [ -0.,   5.]])
        
        See also: DEGDIST_FROM_GC, CLOSEST_POINT_ON_GC, GC2LATLON,
                  GCARC2LATLON, GC_INTERSECT
    """
    
    # check that latlona0/1 & latlonb0/1 are (N,1)x2 in shape and real
    
    # convert to 2D numpy arrays of floats
    latlona0=np.atleast_2d(np.array(latlona0, dtype=np.float))
    latlona1=np.atleast_2d(np.array(latlona1, dtype=np.float))
    latlonb0=np.atleast_2d(np.array(latlonb0, dtype=np.float))
    latlonb1=np.atleast_2d(np.array(latlonb1, dtype=np.float))
    
    # convert to xyz
    latlona0=geocentric2xyz(latlona0)
    latlona1=geocentric2xyz(latlona1)
    latlonb0=geocentric2xyz(latlonb0)
    latlonb1=geocentric2xyz(latlonb1)
    
    # broadcast them
    latlona0,latlona1,latlonb0,latlonb1=np.broadcast_arrays(latlona0,latlona1,latlonb0,latlonb1)
    
    # 2. get perpendicular M to plane formed by gc A
    M=np.cross(latlona0,latlona1, axis=1)
    
    # 3. get perpendicular N to plane formed by gc B
    N=np.cross(latlonb0,latlonb1, axis=1)
    
    # 4. get perpendicular E to plane formed by gc MN
    # - note that this is fully broadcasted from a0/1 & b0/1
    E=np.cross(M,N, axis=1)
    
    # 5. get intersection point directions relative to gc end points
    # FIXME - endpoints don't count towards intersection
    d=np.sign(inner1d(np.cross(M,latlona0, axis=1),E))
    d+=np.sign(inner1d(np.cross(latlona1,M, axis=1),E))
    d+=np.sign(inner1d(np.cross(N,latlonb0, axis=1),E))
    d+=np.sign(inner1d(np.cross(latlonb1,N, axis=1),E))
    d=d[:,np.newaxis]
    
    # 6. intersects if directions all agree
    # FIXME - endpoints don't count towards intersection
    latloni,_=xyz2geocentric(np.where(d==-4,-E,np.where(d==4,E,np.nan)))
    
    # 7. nans if E is basically zero
    # - this occurs when A & B are equivalent
    bad=np.linalg.norm(E, axis=1)
    bad=bad<10.*np.sqrt(np.finfo(np.float).eps)
    if any(bad):
        latloni[bad,:]=np.nan
    
    return latloni


def geocentric2geographiclat(lat, ecc=8.181919084262149e-02):
    """
    GEOCENTRIC2GEOGRAPHICLAT    Convert latitude from geocentric to geographic
    
        Usage:    lat=geocentric2geographiclat(lat)
                  lat=geocentric2geographiclat(lat,ecc)
    
        Description:
         LAT=GEOCENTRIC2GEOGRAPHICLAT(LAT) converts geocentric latitudes LAT to
         geographic latitudes.  LAT is in degrees.  Assumes the WGS-84 reference
         ellipsoid.
    
         LAT=GEOCENTRIC2GEOGRAPHICLAT(LAT,ECC) specifies the eccentricity for
         the ellipsoid to use in the conversion.
    
        Notes:
         - If the location is not on the surface use GEOCENTRIC2GEOGRAPHIC.
    
        Examples:
         # At what geographic latitude is 45N & 60N?:
         >>> geocentric2geographiclat((45,60))
         array([ 45.19242322,  60.16636419])
        
         # Show the difference in latitudes (geographic pushes to the poles):
         import matplotlib.pyplot as plt
         x=np.arange(-90,91)
         plt.plot(x,geocentric2geographiclat(x)-x)
         plt.xlabel('geocentric latitude (^o)')
         plt.ylabel('geographic adjustment (^o)')
         plt.show()
    
        See also: GEOGRAPHIC2GEOCENTRICLAT, GEOGRAPHICLAT2RADIUS
    """
    
    # check lat & ecc are equal shaped or scalar and real
    # check ecc is in range (0,1]
        
    # convert to numpy array of floats
    lat=np.array(lat, dtype=np.float)
    ecc=np.array(ecc, dtype=np.float)

    # convert to geographic
    lat=np.deg2rad(lat)
    lat=np.rad2deg(np.arctan2(np.sin(lat),(1-ecc**2)*np.cos(lat)))
    
    return lat


def geocentric2xyz(latlon,radius=1.):
    """
    GEOCENTRIC2XYZ    Converts coordinates from geocentric to cartesian
    
        Usage:    xyz=geocentric2xyz(latlon)
                  xyz=geocentric2xyz(latlon,radius)
    
        Description:
         XYZ=GEOCENTRIC2XYZ(LATLON) converts coordinates in geocentric latitude,
         longitude, radius to Earth-centered, Earth-Fixed (ECEF).  LATLON is in
         degrees and should be a Nx2 array.  XYZ will contain unit vectors.
         
         XYZ=GEOCENTRIC2XYZ(LATLON,RADIUS) converts coordinates in geocentric
         latitude, longitude, radius to Earth-centered, Earth-Fixed (ECEF).
         LATLON is in degrees and should be a Nx2 array.  XYZ will match the
         units of RADIUS.
    
        Notes:
         - The ECEF coordinate system has the X axis passing through the
           equator at the prime meridian, the Z axis through the north pole
           and the Y axis through the equator at 90 degrees longitude.
    
        Examples:
         # Test the ECEF definition:
         >>> geocentric2xyz([[0,0],[0,90],[90,0]],[[1],[1],[1]])
         array([[  1.00000000e+00,   0.00000000e+00,   0.00000000e+00],
                [  6.12323400e-17,   1.00000000e+00,   0.00000000e+00],
                [  6.12323400e-17,   0.00000000e+00,   1.00000000e+00]])
                
        See also: XYZ2GEOCENTRIC
    """
    
    # check that latlon is (N,1)x2 and real
    # check that radius is (N,1) vector and real
    
    # convert to 2D numpy array of floats
    latlon=np.atleast_2d(np.array(latlon, dtype=np.float))
    radius=np.atleast_2d(np.array(radius, dtype=np.float))
    
    # force radius as a column vector
    radius=np.reshape(radius,(-1,1))
    
    # convert latlon from degrees to radians
    latlon=np.deg2rad(latlon)
    
    # preallocate
    xyz=np.empty((np.shape(latlon)[0],3))
    
    # convert to cartesian
    xyz[:,[0]]=radius*np.cos(latlon[:,[1]])*np.cos(latlon[:,[0]]);
    xyz[:,[1]]=radius*np.sin(latlon[:,[1]])*np.cos(latlon[:,[0]]);
    xyz[:,[2]]=radius*np.sin(latlon[:,[0]]);
    
    return xyz


#def geographic2enu():


def geographic2geocentriclat(lat, ecc=8.181919084262149e-02):
    """
    GEOGRAPHIC2GEOCENTRICLAT    Convert latitude from geographic to geocentric
    
        Usage:    lat=geographic2geocentriclat(lat)
                  lat=geographic2geocentriclat(lat,ecc)
    
        Description:
         LAT=GEOGRAPHIC2GEOCENTRICLAT(LAT) converts geographic latitudes LAT to
         geocentric latitudes.  LAT is in degrees.  Assumes the WGS-84 reference
         ellipsoid.
    
         LAT=GEOGRAPHIC2GEOCENTRICLAT(LAT,ECC) specifies the eccentricity for
         the ellipsoid to use in the conversion.
    
        Notes:
         - If the location is not on the surface use GEOGRAPHIC2GEOCENTRIC.
    
        Examples:
         # Get the geocentric latitude for St. Louis, MO & Los Alamos, NM:
         >>> geographic2geocentriclat((38.649,35.891))
         array([ 38.46142446,  35.70841385])
    
         % Show the difference in latitudes (geocentric pushes to the equator):
         import matplotlib.pyplot as plt
         x=np.arange(-90,91)
         plt.plot(x,geographic2geocentriclat(x)-x)
         plt.xlabel('geographic latitude (^o)')
         plt.ylabel('geocentric adjustment (^o)')
         plt.show()
    
        See also: GEOCENTRIC2GEOGRAPHICLAT, GEOGRAPHICLAT2RADIUS
    """

    # check lat & ecc are equal shaped or scalar and real
    # check ecc is in range (0,1]

    # convert to numpy array of floats
    lat=np.array(lat, dtype=np.float)
    ecc=np.array(ecc, dtype=np.float)

    # convert to geocentric
    lat=np.deg2rad(lat)
    lat=np.rad2deg(np.arctan2((1-ecc**2)*np.sin(lat),np.cos(lat)))

    return lat


def geographiclat2radius(lat, ellipsoid=None):
    """
    GEOGRAPHICLAT2RADIUS    Returns the radius at a geographic latitude

        Usage:    radius=geographiclat2radius(lat)
                  radius=geographiclat2radius(lat,ellipsoid)

        Description:
         RADIUS=GEOGRAPHICLAT2RADIUS(LAT) returns the radii at geographic
         latitudes in LAT.  LAT must be in degrees.  Assumes the WGS-84
         reference ellipsoid.

         RADIUS=GEOGRAPHICLAT2RADIUS(LAT,ELLIPSOID) allows specifying the
         ellipsoid as [A, F] where the parameters A (equatorial radius in
         kilometers) and F (flattening) must be scalar and float.

        Notes:

        Examples:
         # Get the radius for St. Louis, MO & Los Alamos, NM:
         >>> geographiclat2radius((38.649,35.891))
         array([ 6369.83836784,  6370.82788438])
         
         # What is the radius at the pole & equator?:
         >>> geographiclat2radius((90,0))
         array([ 6356.75231425,  6378.137     ])

        See also: GEOCENTRIC2GEOGRAPHICLAT, GEOGRAPHIC2GEOCENTRICLAT
    """

    # check that lat is real
    # check that ellipsoid is 2 element vector of reals

    # default ellipsoid
    if ellipsoid==None:
        ellipsoid=[6378.137, 1/298.257223563]
    
    # convert to numpy array of floats
    lat=np.array(lat, dtype=np.float)

    # convert to radians
    lat=np.deg2rad(lat)

    # get ellipsoid parameters
    a2=ellipsoid[0]**2
    b2=(ellipsoid[0]*(1-ellipsoid[1]))**2

    # optimization
    a2cos2lat=a2*np.cos(lat)**2;
    b2sin2lat=b2*np.sin(lat)**2;

    # get radius
    radius=np.sqrt((a2*a2cos2lat+b2*b2sin2lat)/(a2cos2lat+b2sin2lat));

    return radius


#def geographic2xyz():


def haversine(latlon0,latlon1):
    """
    HAVERSINE    Returns distance between 2 points using the Haversine formula

        Usage:    gcdist=haversine(latlon0,latlon1)

        Description:
         GCDIST=HAVERSINE(LATLON0,LATLON1) finds the spherical great-circle-arc
         degree distance between two points GCDIST.  LATLON0 & LATLON1 must all
         be in degrees and the latitudes must be geocentric.

        Notes:
         - 'half versed sine' is better suited for accuracy at small distances
           compared to SPHERICALINV as it uses a half-versed sine function
           rather than a cosine which becomes inefficient at small distances.

        Examples:
         # Plot distance discrepancy for SPHERICALINV and HAVERSINE:
         d2m=1000*6371*np.pi/180
         d0=10**(np.linspace(-1,7.2,num=83))/d2m
         d0=d0[:,np.newaxis]
         pos=np.zeros((83,2))
         pos[:,[1]]=d0
         dc,_,_=sphericalinv((0,0),pos)
         dh=haversine((0,0),pos)
         hsph,=plt.plot(d2m*d0,d2m*abs(d0-dc),'r',label='SPHERICALINV')
         hhav,=plt.plot(d2m*d0,d2m*abs(d0-dh),'g',label='HAVERSINE')
         plt.xscale('log')
         plt.yscale('log')
         plt.xlabel('distance (m)')
         plt.ylabel('discrepancy (m)')
         plt.legend(handles=[hsph,hhav])
         plt.show()
         # demonstrates convincingly this function is more accurate!
         
         # St. Louis, MO to Yaounde, Cameroon:
         >>> haversine((38.649,-90.305),(3.861,11.521))
         array([[ 96.75578437]])
         
         # St. Louis, MO to Isla Isabella, Galapagos:
         >>> haversine((38.649,-90.305),(-0.823,-91.097))
         array([[ 39.47872366]])
         
         # St. Louis, MO to Los Alamos, NM:
         >>> haversine((38.649,-90.305),(35.891,-106.298))
         array([[ 13.00431661]])

        See also: SPHERICALINV, VINCENTYINV, SPHERICALFWD, VINCENTYFWD
    """
        
    # check that latlon0/latlon1 are (N,1)x2/(N,1)x2 in shape and real
    
    # convert to 2D numpy arrays of floats
    latlon0=np.atleast_2d(np.array(latlon0, dtype=np.float))
    latlon1=np.atleast_2d(np.array(latlon1, dtype=np.float))
    
    # convert from degrees to radians
    latlon0=np.deg2rad(latlon0)
    latlon1=np.deg2rad(latlon1)

    # get haversine distance
    a=np.sin((latlon1[:,[0]]-latlon0[:,[0]])/2)**2  \
     +np.cos(latlon0[:,[0]])*np.cos(latlon1[:,[0]]) \
     *np.sin((latlon1[:,[1]]-latlon0[:,[1]])/2)**2
    gcdist=np.rad2deg(2.*np.arctan2(np.sqrt(a),np.sqrt(1-a)))

    return gcdist


def inlatlonbox(latrng, lonrng, latlon):
    """
        INLATLONBOX    Returns TRUE for positions within the specified region
        
        Usage:    tf=inlatlonbox(latrng,lonrng,latlon)
        
        Description:
        TF=INLATLONBOX(LATRNG,LONRNG,LATLON) returns TRUE or FALSE depending
        on if the positions given by LATLON are within the region specified
        by LATRNG & LONRNG.  LATRNG & LONRNG should be specified as [MIN MAX]
        and LATLON should be given as [LAT LON].  TF is a NROWSx1 logical
        array where NROWS is the number of rows in LATLON.
        
        Notes:
        - Does not call FIXLATLON!
        
        Examples:
        # Make a region covering 1/4th the Earth and see how close to 1/4th
        # of a random set of positions is in the region:
        np.sum(inlatlonbox([0 90],[0 180],randlatlon(1000)))/1000
        
        # Longitude wrapping just works:
        >>> inlatlonbox([0,45],[160, 200],[[30,-170],[0,-180],[20,190]])
        array([ True,  True,  True], dtype=bool)
        >>> inlatlonbox([0,45],[160, 200],[[30,-160],[45,180],[-20,-170]])
        array([ True,  True, False], dtype=bool)
        
        # Test pole handling:
        >>> inlatlonbox([90,90],[180,180],[90,0])
        array([False], dtype=bool)
        >>> inlatlonbox([90,90],[180,180],[90,180])
        array([ True], dtype=bool)
        >>> inlatlonbox([90,90],[180,180],[90,-180])
        array([ True], dtype=bool)
        
        See also: INLONRNG, FIXLATLON, LATMOD, LONMOD, AZINRNG
        """
    
    # check that latrng, lonrng, latlon are numeric & real
    # check that latrng, lonrng, latlon are Nx2 or 1x2
    
    # convert to 2D numpy array of floats
    latrng=np.atleast_2d(np.array(latrng, dtype=np.float))
    lonrng=np.atleast_2d(np.array(lonrng, dtype=np.float))
    latlon=np.atleast_2d(np.array(latlon, dtype=np.float))
    
    # find those in box
    tf=np.logical_and(latlon[:,0]>=latrng[:,0],latlon[:,0]<=latrng[:,1])
    tf=np.logical_and(tf,inlonrng(lonrng,latlon[:,1]))
    
    return tf


def inlonrng(rng, lon):
    """
    INLONRNG    Returns TRUE for longitudes within the specified range
    
        Usage:    tf=inlonrng(rng,lon)
    
        Description:
         TF=INLONRNG(RNG,LON) returns TRUE or FALSE depending on if LON is in
         the longitude range specified by RNG.  RNG must be [MIN MAX] and
         handles wraparound of longitudes.
    
        Notes:
    
        Examples:
         % A few tests that should return TRUE:
         >>> inlonrng([-10, 10],360)
         array([ True], dtype=bool)
         >>> inlonrng([170, 190],-180)
         array([ True], dtype=bool)
         >>> inlonrng([170, 190],180)
         array([ True], dtype=bool)
         >>> inlonrng([-190, -170],-180)
         array([ True], dtype=bool)
         >>> inlonrng([-190, -170],180)
         array([ True], dtype=bool)
         >>> inlonrng([350, 370],0)
         array([ True], dtype=bool)
    
         % A few tougher cases that should return TRUE:
         >>> inlonrng([0, 0],360)
         array([ True], dtype=bool)
         >>> inlonrng([180, 180],-180)
         array([ True], dtype=bool)
         >>> inlonrng([180, 180],180)
         array([ True], dtype=bool)
         >>> inlonrng([-180, -180],-180)
         array([ True], dtype=bool)
         >>> inlonrng([-180, -180],180) # FIXME IT FAILS
         array([ True], dtype=bool)
         >>> inlonrng([360, 360],0)
         array([ True], dtype=bool)
    
         % Yet more cases that should return TRUE:
         >>> inlonrng([0, 360],0)
         array([ True], dtype=bool)
         >>> inlonrng([0, 360],360)
         array([ True], dtype=bool)
         >>> inlonrng([-180, 180],-180)
         array([ True], dtype=bool)
         >>> inlonrng([-180, 180],180)
         array([ True], dtype=bool)
    
        See also: AZINRNG, LONMOD, LATMOD, FIXLATLON, INLATLONBOX
    """

    # check that rng, lon are (N,1)x2 and (N,1)x1 and real
    
    # convert to numpy array of floats
    rng=np.array(rng, dtype=np.float)
    lon=np.array(lon, dtype=np.float)
    
    # convert rng to 2D
    rng=np.atleast_2d(rng)

    # get the effective range
    drng=np.diff(rng,n=1,axis=1)
    drng[np.abs(drng)>360]=360.
    drng[drng<0]+=360
    
    # use lonmod to get things close enough for a simple test
    rng[:,[0]]=lonmod(rng[:,[0]])
    rng[:,[1]]=rng[:,[0]]+drng
    lon=lonmod(lon)
    tf=np.logical_and(lon>=rng[:,0],lon<=rng[:,1]);
    tf=np.logical_or(tf,np.logical_and(lon+360>=rng[:,0],lon+360<=rng[:,1]))
    tf=np.logical_or(tf,np.logical_and(lon-360>=rng[:,0],lon-360<=rng[:,1]))

    return tf


def latmod(lat, wrap=90):
    """
    LATMOD    Returns a latitude modulus (i.e., wraps latitudes into range)
    
        Usage:    wlat,px=latmod(lat)
                  wlat,px=latmod(lat,wrap)
    
        Description:
         WLAT,PX=LATMOD(LAT) returns latitudes LAT to be within the range +/-90
         and  also returns the number of pole-crossings each latitude in LAT
         made in PX.  For example, if LAT=100 & WRAP=90 then WLAT=80 & PX=1.
         This is useful for coupling LATMOD with LONMOD to preserve the
         positions while reducing the values to reasonable ranges.  See the
         Examples section below for instructions on how to use this output.
         Note that LATMOD expects floats and so integers are converted to
         floats such that WLAT and PX are always float.  LAT may be a numpy
         array.
    
         WLAT,PX=LATMOD(LAT,WRAP) is S.*(LAT-N.*WRAP) where N=round(LAT./WRAP)
         if WRAP~=0 and S=1-2.*MOD(N,2).  Thus LATMOD(LAT,WRAP) is always
         within the range +/-WRAP and forms a continuous function (but is
         discontinuous in the 1st derivative) that looks like a "triangle wave".
         The function is primarily intended for wrapping latitude values back
         to a valid range.  The inputs LAT and WRAP must be equal in shape or
         scalar.  WRAP is optional and defaults to 90.
    
        Notes:
         By convention:
          C,PX=latmod(A,0) returns C=A, PX=inf
          C,PX=latmod(A,A) returns C=A, PX=0
          C,PX=latmod(0,0) returns C=0, PX=nan
    
        Examples:
         # Modifying the latitude should also take into account the longitude
         # shift necessary to preserve the actual position.  This may be done
         # by utilizing the second output to shift the longitude by 180 degrees
         # if there are an odd number of pole-crossings:
         lat,px=latmod(lat)
         lon=lonmod(lon+np.mod(px,2)*180)
         
         # Try a few cases:
         >>> latmod((-89,-90,-91))
         (array([-89., -90., -89.]), array([-0.,  0.,  1.]))
         >>> latmod((-269,-270,-271))
         (array([ 89.,  90.,  89.]), array([ 1.,  1.,  2.]))
         >>> latmod((89,90,91))
         (array([ 89.,  90.,  89.]), array([-0.,  0.,  1.]))
         >>> latmod((269,270,271))
         (array([-89., -90., -89.]), array([ 1.,  1.,  2.]))
         
         # Now playing with wrap=0:
         >>> latmod((180, 0, 54.3),0)
         (array([ 180. ,    0. ,   54.3]), array([ inf,  nan,  inf]))
         
        See also: LONMOD, FIXLATLON, INLONRNG, INLATLONBOX
    """
    
    # check that lat & wrap are scalar or equal size and real

    # convert to float
    lat=np.array(lat, dtype=np.float)
    wrap=np.array(wrap, dtype=np.float)

    # modulus without discontinuities
    # i.e., ramp up & ramp down rather than all ramp up
    # i.e., /\/\/\ rather than /|/|/|
    n=np.round(0.5*lat/wrap)
    s=1-2*np.mod(n,2)
    wlat=s*(lat-2*n*wrap)

    # get the number of pole crossings
    px=np.ceil((np.abs(lat)-wrap)/(2*wrap))

    # handle zero modulus case
    d=wrap==0
    if np.any(d):
        if np.size(d)==1:
            wlat=lat
        else:
            wlat[d]=lat[d]

    return wlat, px


def lonmod(lon, wrap=360):
    """
    LONMOD    Returns a longitude modulus (i.e., wraps longitudes into range)

        Usage:    wlon=lonmod(lon)
                  wlon=lonmod(lon,wrap)

        Description:
         WLON=LONMOD(LON) returns the longitudes LON within the range +/-180.
         For example, LON=181 is output as -179.  See the other usage form for
         more algorithm details.  LON must be array_like and composed of floats
         or integers.  WLON is an ndarray of floats.

         WLON=LONMOD(LON,WRAP) is LON-N.*WRAP where N=round(LON./WRAP) if
         WRAP~=0.  Thus LONMOD(LON,WRAP) is always within the range +/-(WRAP/2).
         The inputs LON and WRAP must be equal in shape or scalar.  WRAP is
         optional and defaults to 360.

        Notes:
         - By convention:
            C=lonmod(A,0) returns C=A
            C=lonmod(A,A) returns C=0

        Examples:
         # Try a few cases:
         >>> lonmod((181, 180, -180, -181, 360, -360, 540, -540))
         array([-179.,  180., -180.,  179.,    0.,    0.,  180., -180.])
         
         # Now playing with wrap=0:
         >>> lonmod((180, 0, 54.3),0)
         array([ 180. ,    0. ,   54.3])
         
         # wrap=lon:
         >>> lonmod((180, 0, 54.3),(180, 0, 54.3))
         array([ 0.,  0.,  0.])
         
        See also: LATMOD, FIXLATLON, INLONRNG, INLATLONBOX
    """
    
    # check that lon & wrap are scalar or equal size and real
    
    # convert to float
    lon=np.array(lon, dtype=np.float)
    wrap=np.array(wrap, dtype=np.float)

    # this always flips the sign of -180/180 and so we don't use it
    #wlon=lon-np.round(lon/wrap)*wrap

    # this preserves the sign of -180/180
    wlon=lon-np.sign(lon)*np.ceil((np.abs(lon)-wrap/2)/wrap)*wrap

    # handle zero modulus case
    d=wrap==0
    if np.any(d):
        if np.size(d)==1:
            wlon=lon
        else:
            wlon[d]=lon[d]

    return wlon


def randlatlon(n):
    """
    RANDLATLON    Returns lat/lon points uniformly distributed on a sphere
    
        Usage:    latlon=randlatlon(n)
    
        Description:
         LATLON=RANDLATLON(N) returns the random position of N points uniformly
         distributed on a sphere (not within).  The output LATLON is formatted
         as a Nx2 array of [LAT LON] where the latitudes and longitudes are in
         degrees.  Longitudes are in the range -180 to 180 degrees.
    
        Notes:
    
        Examples:
         # Get 5 random points in latitude & longitude:
         randlatlon(5)
         
         # Map 1000 random lat/lon points:
         p=randlatlon(1000)
         from mpl_toolkits.basemap import Basemap
         import matplotlib.pyplot as plt
         map=Basemap(projection='hammer',lon_0=0)
         map.drawcoastlines()
         x,y=map(p[:,1],p[:,0])
         map.plot(x,y,'ro')
         plt.show()
         
         # Output is Nx2:
         >>> np.shape(randlatlon(1000))
         (1000, 2)
    
        See also: RANDSPHERE
    """
    
    # check that n>=1 & is scalar integer

    # get 3D normal distribution
    xyz=np.random.normal(size=(n, 3))

    # normalize so all lie on the sphere
    xyz/=np.linalg.norm(xyz, axis=1)[:, np.newaxis]
    
    # convert to latlon
    latlon,_=xyz2geocentric(xyz)

    return latlon


def randsphere(n, dim=3):
    """
    RANDSPHERE    Returns points uniformly distributed within a sphere
    
        Usage:    x=randsphere(n)
                  x=randsphere(n,dim)
        
        Description:
         X=RANDSPHERE(N) returns the random positions of N points uniformly
         distributed within a unit sphere.  The output X is a Nx3 array of
         cartesian coordinates.
         
         X=RANDSPHERE(N,DIM) allows getting points in an DIM dimensional sphere.
         DIM should be a scalar integer >=2.
        
        Notes:
        
        Examples:
         # Return 1000 points and plot in 3D:
         x=randsphere(1000)
         import matplotlib.pyplot as plt
         from mpl_toolkits.mplot3d import Axes3D
         fig = plt.figure()
         ax = Axes3D(fig)
         ax.scatter(x[:,0],x[:,1],x[:,2])
         plt.show()
         
         # output is NxM
         >>> np.shape(randsphere(1000))
         (1000, 3)
         >>> np.shape(randsphere(1000,5))
         (1000, 5)
        
        See also: RANDLATLON
    """
    
    # check that n >=1 & dim >=2 and are scalar integers

    # get N-D normal distribution
    x=np.random.normal(size=(n, dim))

    # normalize so all lie on the hypersphere
    x/=np.linalg.norm(x, axis=1)[:, np.newaxis]

    # multiply by a cube-root distribution from 0-1 in radius
    x*=np.power(np.random.uniform(size=(n, 1)), 1./dim)

    return x


#def sph_poly_in():


#def sph_poly_intersect():


def sphericalfwd(latlon0,gcdist,az):
    """
    SPHERICALFWD    Finds a point on a sphere relative to another point

        Usage:    latlon1,baz=sphericalfwd(latlon0,gcdist,az)

        Description:
         LATLON1,BAZ=SPHERICALFWD(LATLON0,GCDIST,AZ) finds geocentric latitudes
         and longitudes in LATLON1 of destination point(s), as well as the back
         azimuths BAZ, given the great circle distances GCDIST and forward
         azimuths AZ from initial point(s) with geocentric latitudes and
         longitudes in LATLON0.  Inputs must all be in degrees.  Outputs are
         also all in degrees.

        Notes:
         - Latitudes are geocentric (0 deg lat == equator, range +/-90)
         - Longitudes are returned in the range +/-180
         - Backazimuth is returned in the range 0-360

        Examples:
         # Heading 45deg NW of St. Louis, MO to ???:
         >>> sphericalfwd((38.649,-90.305),45,-30)
         (array([[  66.90805234, -154.65352366]]), array([[ 95.35931627]]))
         
         # Go 10deg South of Los Alamos, NM to ???:
         >>> sphericalfwd((35.891,106.298),10,180)
         (array([[  25.891,  106.298]]), array([[ 0.]]))

        See also: SPHERICALINV
    """
    
    # check that latlon0/gcdist/az are (N,1)x2/(N,1)/(N,1) in shape and real
    
    # convert to 2D numpy arrays of floats
    latlon0=np.atleast_2d(np.array(latlon0, dtype=np.float))
    gcdist=np.atleast_2d(np.array(gcdist, dtype=np.float))
    az=np.atleast_2d(np.array(az, dtype=np.float))
    
    # force gcdist & az to column vector
    gcdist=np.reshape(gcdist,[-1,1])
    az=np.reshape(az,[-1,1])
    
    # convert from degrees to radians
    latlon0=np.deg2rad(latlon0)
    gcdist=np.deg2rad(gcdist)
    az=np.deg2rad(az)
    
    # preallocate
    latlon1=np.empty(np.shape(latlon0))

    # optimize computation over memory
    sinlat0=np.sin(latlon0[:,[0]])
    singc=np.sin(gcdist)
    coslat0=np.cos(latlon0[:,[0]])
    cosgc=np.cos(gcdist)

    # get destination point
    # - use real on latitudes to avoid occasional
    #   complex value when points coincide
    latlon1[:,[0]]=np.real(np.arcsin(sinlat0*cosgc+coslat0*singc*np.cos(az)))
    sinlat1=np.sin(latlon1[:,[0]])
    latlon1[:,[1]]=latlon0[:,[1]]+np.arctan2(np.sin(az)*singc*coslat0,
                                             cosgc-sinlat0*sinlat1)

    # get back azimuth
    baz=np.arctan2(np.sin(latlon0[:,[1]]-latlon1[:,[1]])*coslat0,
                   np.cos(latlon1[:,[0]])*sinlat0
                   -sinlat1*coslat0*np.cos(latlon0[:,[1]]-latlon1[:,[1]]))
    
    # convert back to degrees and get in proper range
    latlon1=np.rad2deg(latlon1)
    latlon1[:,[1]]=lonmod(latlon1[:,[1]])
    baz=np.mod(np.rad2deg(baz),360)

    return latlon1, baz


def sphericalinv(latlon0,latlon1):
    """
    SPHERICALINV    Return distance and azimuths between 2 locations on sphere

        Usage:    gcdist,az,baz=sphericalinv(latlon0,latlon1)

        Description:
         GCDIST,AZ,BAZ=SPHERICALINV(LATLON0,LATLON1) returns the
         great-circle-arc degree distances GCDIST, forward azimuths AZ and
         backward azimuths BAZ between initial point(s) with geocentric
         latitudes and longitudes LATLON0 and final point(s) with geocentric
         latitudes and longitudes LATLON1 on a sphere.  Inputs must be in
         degrees.  Outputs are also all in degrees.

        Notes:
         - Will always return the shorter distance (GCDIST<=180)
         - GCDIST accuracy degrades when < 3000km (see HAVERSINE example!)
         - Azimuths are returned in the range 0<=az<=360

        Examples:
         # St. Louis, MO to Yaounde, Cameroon:
         >>> sphericalinv((38.649,-90.305),(3.861,11.521))
         (array([[ 96.75578437]]), array([[ 79.53972827]]), array([[ 309.66814964]]))

         # St. Louis, MO to Isla Isabella, Galapagos:
         >>> sphericalinv((38.649,-90.305),(-0.823,-91.097))
         (array([[ 39.47872366]]), array([[ 181.24562107]]), array([[ 0.97288389]]))
         
         # St. Louis, MO to Los Alamos, NM:
         >>> sphericalinv((38.649,-90.305),(35.891,-106.298))
         (array([[ 13.00431661]]), array([[ 262.71497244]]), array([[ 72.98729659]]))

        See also: SPHERICALFWD
    """
    
    # check that latlon0/latlon1 are (N,1)x2/(N,1)x2 in shape and real
    
    # convert to 2D numpy arrays of floats
    latlon0=np.atleast_2d(np.array(latlon0, dtype=np.float))
    latlon1=np.atleast_2d(np.array(latlon1, dtype=np.float))
    
    # convert from degrees to radians
    latlon0=np.deg2rad(latlon0)
    latlon1=np.deg2rad(latlon1)
    
    # optimize for minimum computation at the expense of more memory usage
    sinlat0=np.sin(latlon0[:,[0]])
    sinlat1=np.sin(latlon1[:,[0]])
    coslat0=np.cos(latlon0[:,[0]])
    coslat1=np.cos(latlon1[:,[0]])
    coslo=np.cos(latlon1[:,[1]]-latlon0[:,[1]])
    sinlo=np.sin(latlon1[:,[1]]-latlon0[:,[1]])

    # get law-of-cosines distance
    # - use real to avoid occasional complex values when points coincide
    gcdist=np.real(np.arccos(sinlat0*sinlat1+coslat0*coslat1*coslo))

    # azimuths
    az=np.arctan2(sinlo*coslat1,coslat0*sinlat1-sinlat0*coslat1*coslo)
    baz=np.arctan2(-sinlo*coslat0,coslat1*sinlat0-sinlat1*coslat0*coslo)

    # convert back to degrees and get in proper range
    gcdist=np.rad2deg(gcdist)
    az=np.mod(np.rad2deg(az),360)
    baz=np.mod(np.rad2deg(baz),360)

    # force equal points to be 0,0,0
    eqpo=np.logical_and(latlon0[:,[0]]==latlon1[:,[0]],
                        latlon0[:,[1]]==latlon1[:,[1]])
    if any(eqpo):
        gcdist[eqpo]=0.
        az[eqpo]=0.
        baz[eqpo]=0.

    return gcdist, az, baz


#def vincentyfwd():


#def vincentyinv():


def xyz2geocentric(xyz):
    """
    XYZ2GEOCENTRIC    Converts coordinates from cartesian to geocentric

        Usage:    latlon,radius=xyz2geocentric(xyz)

        Description:
         LATLON,RADIUS=XYZ2GEOCENTRIC(XYZ) converts arrays of coordinates in
         Earth-centered, Earth-Fixed (ECEF) to geocentric latitude, longitude, &
         radius.  LAT and LON are in degrees.  X, Y and Z of XYZ (a Nx3 array of
         floats) must have the same units (RADIUS will will be in those units).

        Notes:
         - The ECEF coordinate system has the X axis passing through the equator
           at the prime meridian, the Z axis through the north pole and the Y
           axis through the equator at 90 degrees longitude.

        Examples:
         # Test the ECEF definition:
         >>> xyz2geocentric([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
         (array([[  0.,   0.],
                [  0.,  90.],
                [ 90.,   0.]]), array([[ 1.],
                [ 1.],
                [ 1.]]))

        See also: GEOCENTRIC2XYZ
    """
    
    # check that xyz is Nx3 and real

    # convert to 2D numpy array of floats
    xyz=np.atleast_2d(np.array(xyz, dtype=np.float))
    
    # preallocate
    latlon=np.empty((np.shape(xyz)[0],2))

    # convert to geocentric
    radius=np.sqrt(xyz[:,[0]]**2+xyz[:,[1]]**2+xyz[:,[2]]**2)
    latlon[:,[0]]=np.rad2deg(np.arcsin(xyz[:,[2]]/radius))
    latlon[:,[1]]=np.rad2deg(np.arctan2(xyz[:,[1]],xyz[:,[0]]))

    return latlon, radius


#def xyz2geographic(xyz):


# for testing functions using "python latlon.py -v"
if __name__ == "__main__":
    import doctest
    doctest.testmod()


