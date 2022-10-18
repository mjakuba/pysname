import numpy as np

ELLIPSOID_WGS84 = (6.378137e6, 298.257223563) # (equatorial axis, inverse flattening)

# State is a 13x1 vector:
# [u v w r p q xe ye ze q0 q1 q2 q3]'
# where:
# [u v w r p q]' [m/s, rad/s] describe the body-fixed linear and angular velocities in std SNAME notation
# [xe ye ze]' [m] are ECEF coordinates (earth-centered-earth-fixed).  We assume these can be treated as 
#          inertial for the purposes of vehicle control.  These will be rarely accessed for the purposes
#          of navigation.
# [q0 q1 q2 q3] a unit quaternion describing the orientiation of the vehicle relative to ECEF coordinates.
#
# methods are defined to permit forward and inverse computation of various derived quantities:
#
# latitude,longitude [rad]: geodetic lat/lon (using the WGS84 ellipsoid as the datum)
# north,east,down [m]: north,east,down coordinates (locally tangent to the WGS84 ellipsoid)  This is an orthographic projection, _not_ a mercator or transverse mercator projection.
#
#
# for now depth is equivalent to down, that is, we assume the earth's surface == the ellipsoid 
# (i.e. we use no geoid approximation to the earth's surface)
#
# The methods and notation here are that of Fossen 2002, mostly.
#
# Unfortunately this is very slow, but I don't understand why.  In ipython try 
#  import state
#  s = state.State()
#  %prun s.set_latitude(0.45) 
# This results in a large number of function calls (as it should from the math) but it's not clear if the get/set stuff implemented here is super 
# inefficient or what.  Unfortunately the interface is really the only nice part of this module.  If speeding it up destroys the interface then it has
# no value.  
# Also try
# %lprun -m state s.set_latitude(0.45) 
# %lprun -f state.State.get_latitude s.latitude = 0.45
#
# Only thought is this might be faster if instead of e.g. self.xe we called self[6].  If that doesn't work then there's probably 
# something more fundamental.  Tried that on one line where it was easy and code was actually slower.

# geographic distance between two points
def dist(s1,s2):

    return(np.linalg.norm(s1.pe-s2.pe,'fro'))

# 2022-10-17  MVJ  Never used these and seem useless.
# # ECEF azimuth of s relative to s0
# def azimuth(s0,s):
#
#     return(np.arctan2(s.ye-s0.ye,s.xe-s0.xe))
#
# # ECEF elevation of s relative to s0
# def elevation(s0,s):
#
#     dr = np.linalg.norm(s.pe[0:2]-s0.pe[:,2],'fro')
#     dz = s.ze-s0.ze
#     return(np.arctan(dz/dr))

def ull(s0,s):

    '''unit vector to s in local-level (NED) coordinates with origin at s0'''
    
    lat0 = s0.latitude
    lon0 = s0.longitude
    dn = s.get_north(lat0,lon0)
    de = s.get_east(lat0,lon0)
    dd = s.get_down(lat0,lon0)

    return(np.mat([dn,de,dd]).T/dist(s0,s))

def azimuth(s0,s):

    '''celestial azimuth of s1 with origin at s0'''

    u = ull(s0,s)
    return(np.arctan2(u[1],u[0]))

def elevation(s0,s):

    '''celestial elevation of s1 with origin at s0; positive for dz<0'''    
    
    u = ull(s0,s)
    dr = np.linalg.norm(u[:2])
    dz = u[2]
    return(np.arctan(-dz/dr))
    
class State(np.matrixlib.defmatrix.matrix):

    def __new__(self,*args):
        obj = super(State,self).__new__(self,np.empty((13,1,)),dtype=float)

        return obj

    def __init__(self,x0=None):

        # Initialize.  Avoid singularity at center of earth 
        # (singularity applies only to derived geodetic 
        # coordinates, not underlying representation)
        if (type(x0) == type(self)):
            self[:] = x0
        elif x0 is None:
            self.fill(0)
            self.xe=ELLIPSOID_WGS84[0]
            self.q[0,0]=1
        else:
            raise(BaseException)

    # Print in human-readable units
    def __str__(self):
        # @@@ do a better job with this later.
        
        # pick a "reasonable" origin for the purposes of displaying northing/easting.
        lat0 = np.fix(self.latitude*180/np.pi)*np.pi/180
        lon0 = np.fix(self.longitude*180/np.pi)*np.pi/180

        latdeg=int(np.fix(self.latitude*180/np.pi))
        latmin=(self.latitude-latdeg*np.pi/180)*180/np.pi*60
        londeg=int(np.fix(self.longitude*180/np.pi))
        lonmin=(self.longitude-londeg*np.pi/180)*180/np.pi*60

        latdec = self.latitude*180/np.pi
        londec = self.longitude*180/np.pi

        return( 

            "latitudedec={latdec: .5f} "
            "longitudedec={londec: .5f}  "
            "latitude={latdeg: 02d}{latmin: .6f} "
            "longitude={londeg: 03d}{lonmin: .6f}  "
            "northing={northing: .2f} "
            "easting={easting: .2f} "
            "depth={self.depth: .2f} "
            "heading={heading: .2f}".format(self=self,
                                            latdec=latdec,londec=londec,
                                            latdeg=latdeg,londeg=londeg,
                                            latmin=latmin,lonmin=lonmin,
                                            northing=self.get_north(lat0,lon0),
                                            easting=self.get_east(lat0,lon0),
                                            heading=self.heading*180.0/np.pi)
            )

    def _asmatrix(self,sub=None):

        """ Create a regular np.mat from all or a portion of a state.State """

        if sub is None:
            return(np.mat(np.asarray(self)))
        else:
            return(np.mat(np.asarray(sub)))
            

    # nu: vector of body-fixed velocities
    def get_nu(self):
        return(self._asmatrix(self[0:6,0]))
    def set_nu(self,nu):
        self[0:6,0]=nu
    nu=property(get_nu,set_nu)

    # vb: vector of body-fixed linear velocities
    def get_vb(self):
        return(self._asmatrix(self[0:3,0]))
    def set_vb(self,vb):
        self[0:3,0]=vb
    vb=property(get_vb,set_vb)

    # omegab: vector of body-fixed angular rates
    def get_omegab(self):
        return(self._asmatrix(self[3:6,0]))
    def set_omegab(self,omegab):
        self[3:6,0]=omegab
    omegab=property(get_omegab,set_omegab)


    # ECEF coordinates
    def get_pe(self):
        return(self._asmatrix(self[6:9,0]))
    def set_pe(self,pe):
        self[6:9,0]=pe
    pe=property(get_pe,set_pe)

    def get_xe(self):
        return(self.pe[0,0])
    def set_xe(self,val):
        self.pe[0,0]=val
    xe=property(get_xe,set_xe)

    def get_ye(self):
        return(self.pe[1,0])
    def set_ye(self,val):
        self.pe[1,0]=val
    ye=property(get_ye,set_ye)

    def get_ze(self):
        return(self.pe[2,0])
    def set_ze(self,val):
        self.pe[2,0]=val
    ze=property(get_ze,set_ze)

    # Geodetic latitude and longitude and ellipsoidal hieght
    def get_latitude(self):

        f = 1/ELLIPSOID_WGS84[1]
        re = ELLIPSOID_WGS84[0]
        rp = ELLIPSOID_WGS84[0]*(1-f)

        e2 = 2*f-f**2
        ep2 = f*(2-f)/(1-f)**2

        # Bowring's one-step method (cm accuracy for |heights| < 1000 km)
        p = np.sqrt(self.xe**2 + self.ye**2)
        theta = np.arctan((self.ze*re)/(p*rp))
        mu = np.arctan((self.ze + ep2*rp*np.sin(theta)**3)/(p - e2*re*np.cos(theta)**3))

        return(mu)

    def set_latitude(self,latitude):

        f = 1/ELLIPSOID_WGS84[1]
        re = ELLIPSOID_WGS84[0]
        rp = ELLIPSOID_WGS84[0]*(1-f)

        mu = latitude
        h = self.height  # calls get_latitude and so uses latitude value prior to it being reset.
        N = re**2/np.sqrt(re**2*np.cos(mu)**2+rp**2*np.sin(mu)**2)

        l = self.longitude

        self.xe = (N+h)*np.cos(mu)*np.cos(l)
        self.ye = (N+h)*np.cos(mu)*np.sin(l)
        self.ze = (rp**2/re**2*N+h)*np.sin(mu)

    latitude=property(get_latitude,set_latitude)

    def get_longitude(self):
        l = np.arctan2(self.ye,self.xe)
        return(l)
    def set_longitude(self,longitude):
        
        f = 1/ELLIPSOID_WGS84[1]
        re = ELLIPSOID_WGS84[0]
        rp = ELLIPSOID_WGS84[0]*(1-f)

        mu = self.latitude
        h = self.height
        N = re**2/np.sqrt(re**2*np.cos(mu)**2+rp**2*np.sin(mu)**2)

        l = longitude

        self.xe = (N+h)*np.cos(mu)*np.cos(l)
        self.ye = (N+h)*np.cos(mu)*np.sin(l)

    longitude=property(get_longitude,set_longitude)

    def get_height(self):

        f = 1/ELLIPSOID_WGS84[1]
        re = ELLIPSOID_WGS84[0]
        rp = ELLIPSOID_WGS84[0]*(1-f)

        mu = self.latitude
        N = re**2/np.sqrt(re**2*np.cos(mu)**2+rp**2*np.sin(mu)**2)
        p = np.sqrt(self.xe**2+self.ye**2)

        h = p/np.cos(mu) - N

        return(h)

    def set_height(self,height):
        
        f = 1/ELLIPSOID_WGS84[1]
        re = ELLIPSOID_WGS84[0]
        rp = ELLIPSOID_WGS84[0]*(1-f)

        mu = self.latitude
        l = self.longitude
        N = re**2/np.sqrt(re**2*np.cos(mu)**2+rp**2*np.sin(mu)**2)

        h = height

        self.xe = (N+h)*np.cos(mu)*np.cos(l)
        self.ye = (N+h)*np.cos(mu)*np.sin(l)
        self.ze = (rp**2/re**2*N+h)*np.sin(mu)

    height=property(get_height,set_height)

    # depth: Assumes earth's surface is equivalent to the ellipsoid, and no geoid approximation is used.
    def get_depth(self):
        return(-self.height)
    def set_depth(self,depth):
        self.height = -depth
    depth=property(get_depth,set_depth)

    # NED coordinates (orthographic coordinates, with origin on WGS84 ellipsoid)
    # note north,east,down properties are useful only to perform a motion,
    # e.g. s.north=45 will move the vehicle 45 m north in the NED frame.  
    # getting the values, as in s.north, is useless as these are always 0, except for s.down.
    # To view orthographic coordinates relative to a particular origin use, e.g., s.get_north(lat,lon)
    def _get_Rne(self):
        return(self.get_Rne(self.latitude,self.longitude))
    def get_Rne(self,lat,lon):
        return(np.transpose(self.get_Ren(lat,lon)))
    def _get_Ren(self):
        return(self.get_Ren(self.latitude,self.longitude))
    def get_Ren(self,lat,lon):
        
        cm = np.cos(lat)
        sm = np.sin(lat)
        cl = np.cos(lon)
        sl = np.sin(lon)

        R = np.mat([[-cl*sm, -sl, -cl*cm],
                    [-sl*sm, cl, -sl*cm],
                    [cm, 0, -sm]])
        
        return(R)

    Ren = property(_get_Ren,None)
    Rne = property(_get_Rne,None)

    def _get_pn(self):
        return(self.get_pn(self.latitude,self.longitude))
    def get_pn(self,lat,lon):

        R = self.get_Rne(lat,lon)

        s0=State()
        s0.latitude=lat
        s0.longitude=lon


        pn = R*(self.pe-s0.pe)
        return(pn)
    
    def _set_pn(self,pn):
        self.set_pn(pn,self.latitude,self.longitude)
    def set_pn(self,pn,lat,lon):

        s0=State()
        s0.latitude=lat
        s0.longitude=lon

        R = self.get_Ren(lat,lon)
        self.pe = R*pn + s0.pe

    pn=property(_get_pn,_set_pn)

    def _get_north(self):
        return(self.get_north(self.latitude,self.longitude))
    def get_north(self,lat,lon):
        return(self.get_pn(lat,lon)[0,0])
    def _set_north(self,north):
        self.set_north(north,self.latitude,self.longitude)
    def set_north(self,north,lat,lon):
        pn = np.mat(self.get_pn(lat,lon))
        pn[0,0]=north
        self.set_pn(pn,lat,lon)
    north=property(_get_north,_set_north)

    def _get_east(self):
        return(self.get_east(self.latitude,self.longitude))
    def get_east(self,lat,lon):
        return(self.get_pn(lat,lon)[1,0])
    def _set_east(self,east):
        self.set_east(east,self.latitude,self.longitude)
    def set_east(self,east,lat,lon):
        pn = np.mat(self.get_pn(lat,lon))
        pn[1,0]=east
        self.set_pn(pn,lat,lon)
    east=property(_get_east,_set_east)
    
    def _get_down(self):
        return(self.get_down(self.latitude,self.longitude))
    def get_down(self,lat,lon):
        return(self.get_pn(lat,lon)[2,0])
    def _set_down(self,down):
        self.set_down(down,self.latitude,self.longitude)
    def set_down(self,down,lat,lon):
        pn = np.mat(self.get_pn(lat,lon))
        pn[2,0]=down
        self.set_pn(pn,lat,lon)
    down=property(_get_down,_set_down)
    
    # eta: vector composed of Position (NED) and attitude (Euler angles)
    # not implemented.
    def get_eta(self):
        pass
    def set_eta(self,eta):
        pass
    eta=property(get_eta,set_eta)    

    # Elements of nu
    def get_u(self):
        return(self.nu[0,0])
    def set_u(self,val):
        self.nu[0,0]=val
    u=property(get_u,set_u)
    speed=property(get_u,set_u)

    # v,w,p,q,r not implemented...

    # Quaternions.
    def get_q(self):
        return(self._asmatrix(self[9:13,0]))
    def set_q(self,q):
        self[9:13,0] = q
    q=property(get_q,set_q)

    def get_Rbn(self):
        return(np.transpose(self.get_Rnb()))
    def get_Rnb(self):

        h = self.q[0,0]
        e1 = self.q[1,0]
        e2 = self.q[2,0]
        e3 = self.q[3,0]
        Se = np.mat([[0, -e3, e2],
                     [e3, 0, -e1],
                     [-e2, e1, 0]])
        R = np.eye(3) + 2*h*Se + 2*Se**2

        return(R)
    Rbn = property(get_Rbn,None)
    Rnb = property(get_Rnb,None)

    # Euler angles.
    def get_euler(self):
        
        R = self.get_Rnb()

        phi = np.arctan2(R[2,1],R[2,2])
        theta = -np.arcsin(R[2,0])
        psi = np.arctan2(R[1,0],R[0,0])

        euler = np.mat([[phi,],[theta,],[psi,]])
        return(euler)

    def set_euler(self,euler):

        phi = euler[0,0]
        theta = euler[1,0]
        psi = euler[2,0]

        cf=np.cos(phi)
        sf=np.sin(phi)
        ct=np.cos(theta)
        st=np.sin(theta)
        cs=np.cos(psi)
        ss=np.sin(psi)

        R = np.mat([[cs*ct, -ss*cf+cs*st*sf, ss*sf+cs*cf*st],
                    [ss*ct, cs*cf+sf*st*ss, -cs*sf+st*ss*cf],
                    [-st, ct*sf, ct*cf]])
        
        R44 = np.trace(R)
        i = np.argmax([R[0,0],R[1,1],R[2,2],R44])
        if i==3:
            pi = np.abs(np.sqrt(1+R44))
        else:
            pi = np.abs(np.sqrt(1+2*R[i,i]-R44))
        
        p4p1 = R[2,1]-R[1,2]
        p4p2 = R[0,2]-R[2,0]
        p4p3 = R[1,0]-R[0,1]
        p2p3 = R[2,1]+R[1,2]
        p3p1 = R[0,2]+R[2,0]
        p1p2 = R[1,0]+R[0,1]

        if i==0:
            p1=pi
            p2=p1p2/p1
            p3=p3p1/p1
            p4=p4p1/p1
        elif i==1:
            p2=pi
            p1=p1p2/p2
            p3=p2p3/p2
            p4=p4p2/p2
        elif i==2:
            p3=pi
            p1=p3p1/p3
            p2=p2p3/p3
            p4=p4p3/p3
        elif i==3:
            p4=pi
            p1=p4p1/p4
            p2=p4p2/p4
            p3=p4p3/p4

        self.q[0,0] = p4/2
        self.q[1,0] = p1/2
        self.q[2,0] = p2/2
        self.q[3,0] = p3/2

    euler=property(get_euler,set_euler)

    def get_psi(self):
        return(self.euler[2,0])
    def set_psi(self,psi):
        euler = np.mat(self.euler)
        euler[2,0]=psi
        self.set_euler(euler)
    psi=property(get_psi,set_psi)
    heading=property(get_psi,set_psi)

    def get_theta(self):
        return(self.euler[1,0])
    def set_theta(self,theta):
        euler = np.mat(self.euler)
        euler[1,0]=theta
        self.set_euler(euler)
    theta=property(get_theta,set_theta)
    pitch=property(get_theta,set_theta)

    def get_phi(self):
        return(self.euler[0,0])
    def set_phi(self,phi):
        euler = np.mat(self.euler)
        euler[0,0]=phi
        self.set_euler(euler)
    phi=property(get_phi,set_phi)
    roll=property(get_phi,set_phi)
