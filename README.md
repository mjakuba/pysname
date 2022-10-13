# sname-state

Vehicle/vessel state and conversions with SNAME notation for naval architecture and marine robotics applications.

## Installation

    pip3 install -e .

## Definitions and notation
    # State is a 13x1 vector:
    # [u v w r p q xe ye ze q0 q1 q2 q3]'
    # where:
    # [u v w r p q]' [m/s, rad/s] describe the body-fixed linear and angular velocities in std SNAME notation
    # [xe ye ze]' [m] are ECEF coordinates (earth-centered-earth-fixed).  We assume these can be treated as 
    #          inertial for the purposes of vehicle control.  These will be rarely accessed for the purposes
    #          of navigation.
    # [q0 q1 q2 q3] a unit quaternion describing the orientiation of the vehicle relative to ECEF coordinates.
    
## Conversions between coordinate frames
    
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
   
## Performance

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
